#!/usr/bin/env python3
"""
segment_cnv_changepoint_mp.py
-----------------------------
Parallel CNV caller (ruptures + PELT) with asymmetric post-processing
Comai Laboratory, UC Davis Genome Center, 2025


Stages:
1) PELT segmentation (global optimum, RBF model)
2) Label segments (medians vs del-cut / ins-cut)
3) Merge same-type segments across short gaps
   • deletions use --gap-bins-del (default: --gap-bins)
   • insertions use --gap-bins-ins (default: 10)
4) Promote core runs inside Neutral
   • deletions: threshold = --post-core-cut, run >= --post-run-bins
   • insertions: threshold = --post-ins-core-cut, run >= --post-ins-run-bins
5) Extend edges while local core density stays high
   • deletions: frac = --extend-del-frac (default: --extend-frac)
   • insertions: frac = --extend-ins-frac (default: 0.75)
6) Insertion quality gate (only): require both
   • length >= --min-ins-length-bins
   • core_frac >= --min-ins-core-frac  (within the segment)
"""

import argparse, zipfile, numpy as np, pandas as pd, ruptures as rpt
from scipy import stats, linalg
from multiprocessing import Pool, cpu_count
from functools import partial


# ───────────── I/O helpers ────────────────────────────────────────────
def load_bins(path):
    if path.lower().endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            name = next(f for f in z.namelist()
                        if f.lower().endswith(('.tsv', '.txt')))
            with z.open(name) as fh:
                df = pd.read_csv(fh, sep='\t', dtype=str)
    else:
        df = pd.read_csv(path, sep='\t', dtype=str)
    df.rename(columns={df.columns[0]: 'chrom',
                       df.columns[1]: 'start',
                       df.columns[2]: 'end'}, inplace=True)
    df[['start','end']] = df[['start','end']].apply(pd.to_numeric)
    return df


def neutral_mode(ratios):
    band = ratios[(0.8 <= ratios) & (ratios <= 1.2)]
    if len(band) == 0:
        return 1.0
    try:
        kde = stats.gaussian_kde(np.log(band))
        xs = np.linspace(np.log(0.8), np.log(1.2), 200)
        return float(np.exp(xs[np.argmax(kde(xs))]))
    except (linalg.LinAlgError, ValueError):
        hist, edges = np.histogram(band, bins=200, range=(0.8, 1.2))
        i = int(np.argmax(hist))
        return float((edges[i] + edges[i+1]) / 2.0)


# ───────────── segmentation & postproc ────────────────────────────────
def initial_segments(ratios, pen, min_size, neut, dcut, icut):
    algo = rpt.Pelt(model="rbf", min_size=min_size).fit(ratios)
    ends = algo.predict(pen=pen)  # 1-based ends
    segs = []
    s = 0
    for e in ends:
        med = float(np.median(ratios[s:e]))
        if med < dcut * neut:
            typ = "Deletion"
        elif med > icut * neut:
            typ = "Insertion"
        else:
            typ = "Neutral"
        segs.append([s, e-1, typ, med])
        s = e
    return segs


def merge_gaps_asym(segs, gap_del, gap_ins):
    """Merge same-type segments across small gaps.
       Deletions use gap_del; Insertions use gap_ins.
    """
    if not segs:
        return segs
    out = [segs[0]]
    i = 1
    while i < len(segs):
        prev = out[-1]
        cur = segs[i]
        if cur[2] == prev[2]:
            # direct adjacency or tiny neutral gap
            direct_gap = cur[0] - prev[1] - 1
            if prev[2] == "Deletion":
                allow = gap_del
            elif prev[2] == "Insertion":
                allow = gap_ins
            else:
                allow = 0  # Neutral shouldn't merge with Neutral here

            if direct_gap <= allow:
                prev[1] = cur[1]
                i += 1
                continue

            # Neutral bridge case: prev , Neutral , cur (same type)
            if (i + 1 < len(segs) and
                segs[i][2] == "Neutral" and
                segs[i+1][2] == prev[2]):
                bridge_len = segs[i][1] - segs[i][0] + 1
                if prev[2] == "Deletion" and bridge_len <= gap_del:
                    prev[1] = segs[i+1][1]
                    i += 2
                    continue
                if prev[2] == "Insertion" and bridge_len <= gap_ins:
                    prev[1] = segs[i+1][1]
                    i += 2
                    continue

        out.append(cur)
        i += 1
    return out


def scan_neutral_promote(ratios, segs,
                         del_core, del_run,
                         ins_core, ins_run):
    """Promote runs inside Neutral segments to CNV calls."""
    adds = []
    for s, e, typ, _ in segs:
        if typ != "Neutral":
            continue
        # Deletion-like runs (≤ del_core)
        run = 0; rs = None
        for idx in range(s, e+1):
            if ratios[idx] <= del_core:
                if run == 0: rs = idx
                run += 1
            else:
                if run >= del_run:
                    adds.append([rs, idx-1, "Deletion",
                                 float(np.median(ratios[rs:idx]))])
                run = 0
        if run >= del_run:
            adds.append([rs, e, "Deletion",
                         float(np.median(ratios[rs:e+1]))])

        # Insertion-like runs (≥ ins_core)
        run = 0; rs = None
        for idx in range(s, e+1):
            if ratios[idx] >= ins_core:
                if run == 0: rs = idx
                run += 1
            else:
                if run >= ins_run:
                    adds.append([rs, idx-1, "Insertion",
                                 float(np.median(ratios[rs:idx]))])
                run = 0
        if run >= ins_run:
            adds.append([rs, e, "Insertion",
                         float(np.median(ratios[rs:e+1]))])

    segs.extend(adds)
    segs.sort(key=lambda x: x[0])
    return segs


def extend_edges_asym(ratios, segs,
                      del_core, del_win, del_frac,
                      ins_core, ins_win, ins_frac):
    """Extend both edges while local core density stays high.
       Deletions use ≤ del_core with frac del_frac over window del_win.
       Insertions use ≥ ins_core with frac ins_frac over window ins_win.
    """
    n = len(ratios)
    for seg in segs:
        typ = seg[2]
        if typ == "Deletion":
            core_fn = lambda w: (w <= del_core).mean()
            need_win, need_frac = del_win, del_frac
        elif typ == "Insertion":
            core_fn = lambda w: (w >= ins_core).mean()
            need_win, need_frac = ins_win, ins_frac
        else:
            continue

        # Left extension
        s = seg[0]
        while s - need_win >= 0:
            if core_fn(ratios[s-need_win:s]) >= need_frac:
                s -= 1
            else:
                break
        seg[0] = s

        # Right extension
        e = seg[1]
        while e + need_win < n:
            if core_fn(ratios[e+1:e+1+need_win]) >= need_frac:
                e += 1
            else:
                break
        seg[1] = e

    segs.sort(key=lambda x: x[0])
    return segs


# ───────────── worker ─────────────────────────────────────────────────
def process_sample(pair, bins_df, a):
    sample, col = pair
    out_rows = []

    ratios_all = pd.to_numeric(bins_df[col], errors='coerce').fillna(1.0).values
    neut = neutral_mode(ratios_all)

    for chrom, cdf in bins_df.groupby('chrom'):
        cdf = cdf.reset_index(drop=True)
        ratios = pd.to_numeric(cdf[col], errors='coerce').fillna(neut).values

        # Stage 1: segmentation + label
        segs = initial_segments(ratios, a.pen, a.min_size, neut, a.del_cut, a.ins_cut)

        # Stage 2: asymmetric gap merge
        gap_del = a.gap_bins_del if a.gap_bins_del is not None else a.gap_bins
        gap_ins = a.gap_bins_ins if a.gap_bins_ins is not None else 10  # stricter default
        segs = merge_gaps_asym(segs, gap_del, gap_ins)

        # Stage 3a: neutral scan promote (both types)
        ins_core = a.post_ins_core_cut
        ins_run  = a.post_ins_run_bins
        segs = scan_neutral_promote(ratios, segs,
                                    a.post_core_cut, a.post_run_bins,
                                    ins_core, ins_run)

        # Stage 3b: asymmetric edge extension
        del_frac = a.extend_del_frac if a.extend_del_frac is not None else a.extend_frac
        ins_frac = a.extend_ins_frac if a.extend_ins_frac is not None else 0.75  # stricter
        segs = extend_edges_asym(ratios, segs,
                                 a.post_core_cut, a.extend_win, del_frac,
                                 ins_core, a.extend_win, ins_frac)

        # Stage 4: insertion quality gate (length + core fraction)
        for s, e, typ, med in segs:
            if typ == "Neutral":
                continue
            if typ == "Insertion":
                seg = ratios[s:e+1]
                core_frac = float((seg >= ins_core).mean())
                length_bins = e - s + 1
                if (length_bins < a.min_ins_length_bins) or (core_frac < a.min_ins_core_frac):
                    # discard weak/loose insertions
                    continue

            out_rows.append({
                'chrom': chrom,
                'sample': sample,
                'start': int(cdf.loc[s, 'start']),
                'end':   int(cdf.loc[e, 'end']),
                'length_bins': e - s + 1,
                'median_ratio': round(med, 3),
                'type': typ
            })

    return out_rows


# ───────────── CLI / main ────────────────────────────────────────────
def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--bins', required=True)
    p.add_argument('--out',  default='cnv_changepoint.tsv')
    p.add_argument('--pen', type=float, default=25.0)
    p.add_argument('--min-size', type=int, default=5)
    p.add_argument('--del-cut', type=float, default=0.80)
    p.add_argument('--ins-cut', type=float, default=1.15)
    p.add_argument('--gap-bins', type=int, default=30)
    # New: asymmetric gap sizes
    p.add_argument('--gap-bins-del', type=int, default=None,
                   help='Override deletion gap merge (default: --gap-bins)')
    p.add_argument('--gap-bins-ins', type=int, default=None,
                   help='Override insertion gap merge (default: 10)')

    # Second pass thresholds
    p.add_argument('--post-core-cut', type=float, default=0.85)   # deletions
    p.add_argument('--post-run-bins', type=int,   default=10)     # deletions
    p.add_argument('--post-ins-core-cut', type=float, default=1.18,
                   help='Insertion core threshold for promotion/extension')
    p.add_argument('--post-ins-run-bins', type=int, default=12,
                   help='Min consecutive high bins to promote insertion')

    # Edge extension
    p.add_argument('--extend-win',  type=int,   default=10)
    p.add_argument('--extend-frac', type=float, default=0.60)     # deletions default
    p.add_argument('--extend-del-frac', type=float, default=None,
                   help='Override deletion extension frac (default: --extend-frac)')
    p.add_argument('--extend-ins-frac', type=float, default=None,
                   help='Override insertion extension frac (default: 0.75)')

    # Insertion quality gate
    p.add_argument('--min-ins-length-bins', type=int, default=10,
                   help='Drop insertions shorter than this many bins')
    p.add_argument('--min-ins-core-frac', type=float, default=0.65,
                   help='Drop insertions with < this fraction of bins ≥ post-ins-core-cut')

    p.add_argument('--sample-suffix', default='/CONTROL')
    p.add_argument('--nproc', type=int, default=cpu_count())
    return p


def main():
    a = build_parser().parse_args()
    bins_df = load_bins(a.bins)

    # normalized columns only (suffix)
    sample_cols = {col[:-len(a.sample_suffix)]: col
                   for col in bins_df.columns[3:]
                   if col.endswith(a.sample_suffix)}
    if not sample_cols:
        raise SystemExit(f"No columns end with suffix “{a.sample_suffix}”")

    worker = partial(process_sample, bins_df=bins_df, a=a)
    with Pool(a.nproc) as pool:
        chunks = pool.map(worker, sample_cols.items())

    flat = [row for ch in chunks for row in ch]
    pd.DataFrame(flat).sort_values(['chrom','sample','start']) \
        .to_csv(a.out, sep='\t', index=False)
    print(f"Wrote {len(flat)} CNV segments → {a.out}")


if __name__ == '__main__':
    main()


