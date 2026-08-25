# -*- coding: utf-8 -*-
"""Build the RCTD single-cell reference (fibro_immune_ref.h5ad) that 06_spatial expects.
RCTD deconvolves spatial spots against a labelled scRNA reference; 06's _export_reference()
reads layers['counts'] (raw counts) + obs['celltype']. We build a STRATIFIED subsample of the
full annotated atlas (~1500 cells per celltype -> compact + representative; RCTD needs raw
counts, and a per-type subsample is standard practice, not the 424K full atlas). npc env.
Usage: python build_rctd_ref.py"""
import os, sys, io
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings; warnings.filterwarnings('ignore')
import numpy as np, scipy.sparse as sp, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

ATLAS = '/data/ssc/powered/full_integrated_annotated.h5ad'
OUT = '/data/ssc/raw/atlas/fibro_immune_ref.h5ad'
N_PER = 1500          # cells per celltype (cap; types with fewer keep all)
SEED = 0

print('[rctd-ref] loading', ATLAS, flush=True)
A = sc.read_h5ad(ATLAS)
assert 'counts' in A.layers, 'atlas lacks layers["counts"] (raw counts) required by RCTD'
assert 'celltype' in A.obs, 'atlas lacks obs["celltype"]'

ct = A.obs['celltype'].astype(str).values
rng = np.random.default_rng(SEED)
keep = []
for c in np.unique(ct):
    pos = np.where(ct == c)[0]
    take = min(len(pos), N_PER)
    keep.append(rng.choice(pos, size=take, replace=False))
keep = np.sort(np.concatenate(keep))
print('[rctd-ref] stratified subsample: %d / %d cells, %d celltypes' %
      (len(keep), A.n_obs, len(np.unique(ct))), flush=True)

ref = A[keep].copy()
# RCTD needs RAW COUNTS: put counts in BOTH X and layers['counts'] (06 reads layers['counts']).
cnt = ref.layers['counts']
cnt = cnt.copy() if sp.issparse(cnt) else sp.csr_matrix(cnt)
# ★ spacexr/RCTD check_counts(require_int=T) REJECTS non-integer counts; the atlas layers['counts']
# are decontX-CORRECTED counts (fractional) -> round to nearest integer (RCTD's Poisson model
# requires integer counts). Without this the reference fails check_counts and every slice's RCTD dies.
cnt = cnt.astype('float32')
cnt.data = np.rint(cnt.data)
cnt.eliminate_zeros()
ref.X = cnt
ref.layers.clear(); ref.layers['counts'] = cnt.copy()
ref.obs = ref.obs[['celltype']].copy()   # barcodes (index) + celltype label only
ref.raw = None
os.makedirs(os.path.dirname(OUT), exist_ok=True)
ref.write(OUT)
print('[rctd-ref] wrote %s: %d cells, %d types (%s), %d genes' %
      (OUT, ref.n_obs, ref.obs['celltype'].nunique(),
       ', '.join(sorted(ref.obs['celltype'].astype(str).unique())), ref.n_vars), flush=True)
print('RCTD_REF_DONE', flush=True)
