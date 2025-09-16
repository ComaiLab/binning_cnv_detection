# binning_cnv_detection - segment_cnv_changepoint_mp.py
Calls CNVs from bin-by-sam.py per-bin normalized coverage. For each sample×chromosome it segments with PELT, labels segments as Deletion/Insertion by median ratio, merges across small neutral gaps, promotes shallow core runs, and extends edges. Output TSV lists chrom, sample, start/end, length_bins, median_ratio, and type. Code was tuned using 100-kb bins and is written in python3.


## Input & I/O

### `--bins` (required)
Path to a binned coverage table (`.tsv/.txt` or `.zip` containing one).  
- First 3 columns: `chrom`, `start`, `end`.  
- Remaining columns: per-sample normalized ratios.  
- If zipped, the script auto-detects the first `*.tsv`/`*.txt` file inside.

### `--out` (default: `cnv_changepoint.tsv`)
Output file (TSV).  
Columns: `chrom`, `sample`, `start`, `end`, `length_bins`, `median_ratio`, `type`.

### `--sample-suffix` (default: `/CONTROL`)
Only columns ending with this control name suffix are considered as sample coverage tracks. This would be the second set of columns in a normal bin-by-sam output containing only the normalized counts. This is typically either "/NA" if no control or "/Control-Lib-Name" if a control is used 

### `--nproc` (default: all CPUs)
Number of worker processes.  
Match this with `--cpus-per-task` in SLURM.

---

## Global Model & Segmentation (Stage 1)

### Neutral mode (auto-estimated)
Per-sample neutral ratio estimated via KDE (or histogram fallback) on bins within 0.8–1.2.  

### `--pen` (float, default: `25.0`)
PELT penalty term.  
- Higher → fewer, larger segments.  
- Lower → more, smaller segments.

### `--min-size` (int, default: `5`)
Minimum bins per segment in PELT.  

---

## Primary Labeling Thresholds

Segments are compared to the neutral ratio:  

### `--del-cut` (default: `0.80`)
If `median < del_cut * neutral` → **Deletion**.

### `--ins-cut` (default: `1.15`)
If `median > ins_cut * neutral` → **Insertion**.  

Anything in between is **Neutral**.

---

## Gap Merging (Stage 2; Asymmetric)

Merge same-type segments across small gaps or Neutral bridges:  

### `--gap-bins` (default: `30`)
Default allowable gap.  

### `--gap-bins-del` (default: None → `--gap-bins`)
Deletion gap size override.  

### `--gap-bins-ins` (default: None → `10`)
Insertion gap size override (stricter default).  

---

## Neutral Scan & Promotion (Stage 3a)

Detect CNV-like **runs inside Neutral segments**.

### Deletion promotion
- `--post-core-cut` (default: `0.85`) → bins ≤ this ratio count as deletion.  
- `--post-run-bins` (default: `10`) → minimum run length to call a deletion.  

### Insertion promotion
- `--post-ins-core-cut` (default: `1.18`) → bins ≥ this ratio count as insertion.  
- `--post-ins-run-bins` (default: `12`) → minimum run length to call an insertion.  

---

## Edge Extension (Stage 3b)

Extend segment edges while local window density remains strong.

### Shared
- `--extend-win` (default: `10`) → sliding window size in bins.  

### Deletions
- `--extend-frac` (default: `0.60`) → fraction of bins ≤ `post-core-cut`.  
- `--extend-del-frac` (default: None → `--extend-frac`).  

### Insertions
- `--extend-ins-frac` (default: None → `0.75`) → fraction of bins ≥ `post-ins-core-cut`.  

---

## Insertion Quality Gate (Stage 4)

Final filter for insertions. Requires both:  

- `--min-ins-length-bins` (default: `10`) → minimum segment length.  
- `--min-ins-core-frac` (default: `0.65`) → fraction of bins ≥ `post-ins-core-cut`.  

---

## Practical Tuning Tips

- **Fragmented deletions:** raise `--gap-bins-del` or lower `--pen`.  
- **Over-merged insertions:** keep `--gap-bins-ins` small, raise `--extend-ins-frac`.  
- **Missing shallow gains:** lower `--ins-cut`, `--post-ins-core-cut`, and `--post-ins-run-bins`.  
- **Excess noisy insertions:** raise `--min-ins-length-bins` and `--min-ins-core-frac`.  
- **Jittery breakpoints:** increase `--pen` or `--min-size`, and raise extension fractions.  

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bins` | **required** | Input binned coverage table (`.tsv/.txt` or `.zip` containing one). First 3 columns: `chrom`, `start`, `end`; remaining: per-sample normalized ratios. |
| `--out` | `cnv_changepoint.tsv` | Output TSV file with CNV calls (`chrom`, `sample`, `start`, `end`, `length_bins`, `median_ratio`, `type`). |
| `--sample-suffix` | `/CONTROL` | Only columns ending with this suffix are processed as sample coverage tracks, "/NA" if no control, "/name-of-control-lib" if a control was used. This would be the second set of columns in a normal bin-by-sam output containing only the normalized counts. |
| `--nproc` | all CPUs | Number of worker processes (per-sample parallelism). Match with `--cpus-per-task` in SLURM. |
| `--pen` | `25.0` | Penalty for PELT segmentation. Higher → fewer, larger segments; lower → more, smaller segments. |
| `--min-size` | `5` | Minimum number of bins per segment in PELT. |
| `--del-cut` | `0.80` | Segments with `median < del_cut * neutral` → **Deletion**. |
| `--ins-cut` | `1.15` | Segments with `median > ins_cut * neutral` → **Insertion**. |
| `--gap-bins` | `30` | Max gap in bins for merging same-type segments (default for deletions). |
| `--gap-bins-del` | `None` (→ `--gap-bins`) | Override max gap size for **Deletion** merging. |
| `--gap-bins-ins` | `None` (→ `10`) | Override max gap size for **Insertion** merging (default stricter). |
| `--post-core-cut` | `0.85` | Deletion promotion: bins ≤ this ratio count toward a deletion run. |
| `--post-run-bins` | `10` | Deletion promotion: minimum consecutive bins ≤ `post-core-cut`. |
| `--post-ins-core-cut` | `1.18` | Insertion promotion: bins ≥ this ratio count toward an insertion run. |
| `--post-ins-run-bins` | `12` | Insertion promotion: minimum consecutive bins ≥ `post-ins-core-cut`. |
| `--extend-win` | `10` | Window size (bins) for extension at segment edges. |
| `--extend-frac` | `0.60` | Fraction threshold for extending **deletions** (bins ≤ `post-core-cut`). |
| `--extend-del-frac` | `None` (→ `--extend-frac`) | Override extension fraction for deletions only. |
| `--extend-ins-frac` | `None` (→ `0.75`) | Extension fraction for insertions (bins ≥ `post-ins-core-cut`). |
| `--min-ins-length-bins` | `10` | Minimum length (bins) for final insertion calls. |
| `--min-ins-core-frac` | `0.65` | Within insertion segments, required fraction of bins ≥ `post-ins-core-cut`. |

---

## Output Format

| Column         | Description                                                                                                      |
|----------------|------------------------------------------------------------------------------------------------------------------|
| `chrom`        | Chromosome/contig ID.                                                                                            |
| `sample`       | Sample name (suffix removed).                                                                                    |
| `start`, `end` | Genomic coordinates in bp (inclusive).                                                                           |
| `length_bins`  | Number of bins spanned by the CNV segment.                                                                       |
| `median_ratio` | Median normalized coverage within the segment (Neutral ≈ 1.0; Deletions ~0.5–0.9; Insertions ≥1.15).             |
| `type`         | CNV call type: either **Deletion** or **Insertion**.                                                             |

