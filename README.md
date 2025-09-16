# binning_cnv_detection
Calls CNVs from bin-by-sam.py per-bin normalized coverage. For each sample×chromosome it segments with PELT, labels segments as Deletion/Insertion by median ratio, merges across small neutral gaps, promotes shallow core runs, and extends edges. Output TSV lists chrom, sample, start/end, length_bins, median_ratio, and type.
