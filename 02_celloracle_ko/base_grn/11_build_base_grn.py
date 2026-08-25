#!/usr/bin/env python
# ============================================================================
# SKIN base-GRN build, STEP 2 (Python/CellOracle): TSS annotation + motif scan
# Follows CellOracle official tutorials:
#   notebooks/01_ATAC-seq_data_processing/.../02_preprocess_peak_data
#   notebooks/02_motif_scan/02_atac_peaks_to_TFinfo_with_celloracle_20200801
# Output = TF_info_matrix (peak_id, gene_short_name, one column per TF) -> parquet
# ============================================================================
import os, sys
import pandas as pd
import numpy as np
import celloracle as co
from celloracle import motif_analysis as ma

OUT   = "/data/ssc/basegrn/cicero_out"
PARQ  = "/data/ssc/basegrn/skin_base_GRN_dataframe.parquet"
N_CPUS = 10                     # leave ~6 cores for concurrent SCENIC+
REF   = "hg38"

# ---- load Cicero outputs ---------------------------------------------------
peaks = pd.read_csv(os.path.join(OUT, "all_peaks.csv"), index_col=0)
peaks = peaks.x.values                                   # R vector -> column 'x'
cicero_connections = pd.read_csv(os.path.join(OUT, "cicero_connections.csv"), index_col=0)
print("[py] n_peaks:", len(peaks), " n_cicero_conns:", len(cicero_connections), flush=True)

# ---- 1. annotate TSS + integrate co-accessibility --------------------------
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=REF)
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                               cicero_connections=cicero_connections)
print("[py] integrated peaks (pre-filter):", integrated.shape, flush=True)

# ---- 2. coaccess >= 0.8 (tutorial threshold) -------------------------------
peak = integrated[integrated.coaccess >= 0.8]
peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)
peak.to_csv(os.path.join(OUT, "processed_peak_file.csv"))
print("[py] peaks passing coaccess>=0.8:", peak.shape[0],
      " unique genes:", peak.gene_short_name.nunique(), flush=True)

# ---- 3. peak-format QC (drops peaks outside chrom bounds) -------------------
peak = ma.check_peak_format(peak, ref_genome=REF, genomes_dir=None)
print("[py] peaks after check_peak_format:", peak.shape[0], flush=True)

# ---- 4. gimmemotifs motif scan (tutorial params) ---------------------------
tfi = ma.TFinfo(peak_data_frame=peak, ref_genome=REF, genomes_dir=None)
tfi.scan(fpr=0.02, motifs=None, n_cpus=N_CPUS, verbose=True)   # motifs=None -> CellOracle default (gimme CisBP)
tfi.to_hdf5(file_path=os.path.join(OUT, "skin_TFinfo.celloracle.tfinfo"))  # checkpoint
tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=10)                       # tutorial threshold
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

# ---- 5. base-GRN dataframe -> parquet --------------------------------------
df = tfi.to_dataframe()
df.to_parquet(PARQ)
n_tf = df.shape[1] - 2
n_edges = int((df.iloc[:, 2:].values != 0).sum())
print("[py] BASE-GRN rows(peaks-genes):", df.shape[0], " TF cols:", n_tf,
      " TF-peak edges:", n_edges, flush=True)
# gate (b): HES1/SMAD3 etc must be source TFs with a non-trivial target set
for g in ["HES1","SMAD3","FOSB","MEF2C","HEYL","CEBPD"]:
    if g in df.columns:
        nz = df[g].values != 0
        n_rows = int(nz.sum())
        n_genes = int(df.loc[nz, "gene_short_name"].nunique())
        print(f"   TF {g}: present=True  n_target_rows={n_rows}  n_target_genes={n_genes}", flush=True)
    else:
        print(f"   TF {g}: present=False (no motif matched -> NO regulator edges)", flush=True)
print("SAVED", PARQ, flush=True)
