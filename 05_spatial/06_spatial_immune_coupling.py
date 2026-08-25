# -*- coding: utf-8 -*-
"""
======================================================================================
06_spatial_immune_coupling.py   ·   POWERED immune-layer core (Fig B)
SSc virtual-perturbation paper — spatial immune-stromal coupling at REGULON resolution
--------------------------------------------------------------------------------------
WHAT (main line, honest framing):
  The engine (CellOracle in-silico TF-KO + decoupler CollecTRI regulon activity)
  nominated HES1/Notch (-> nirogacestat) and FOSB/AP-1 as drivers of the SSc
  myofibroblast fibrotic program. This script asks, PER SPATIAL SLICE, whether the
  macrophage micro-environment is spatially COUPLED to HES1/FOSB regulon activity in
  the myofibroblast front, and WHICH compartment sources the Notch ligand.

  HONEST FRAMING (baked into comments and outputs — do NOT overclaim):
    * Regulon activity + virtual-KO = HYPOTHESIS + CORRELATION, never causal.
    * Immune coupling = spatial ASSOCIATION / co-context. NEVER "immune -> HES1".
      Ligand sourcing shows proximity, not a proven signalling axis.
    * TREM2-hi / SAM macrophages are PROTECTIVE in SSc skin (Xu et al.), NOT
      pro-fibrotic. They are scored and reported as a protective-axis contrast.
    * nirogacestat has NO immune endpoints — any immune read-out here is a spatial
      hypothesis for the discussion, not a drug-effect claim.

ANALYSES per slice (each macrophage sub-state x each lead regulon):
  1. decoupler ULM regulon activity per spot: HES1, FOSB (+ SMAD3 pos ctrl, MEF2C neg ctrl).
  2. macrophage sub-state scores per spot: SPP1-hi, TREM2-hi/SAM (PROTECTIVE),
     IL1B/FCN1-inflammatory, LYVE1-perivascular  (score_genes; panels reuse the
     figA_macrophage_regulon_coupling.py definitions, scaled up).
  3. RCTD deconvolution (module 505 via rctd_deconvolve.R subprocess) -> myofibroblast
     & macrophage spot fractions; marker-score proxy fallback if spacexr absent.
  4. BANKSY spatial domains (module 541 via banksy_domains.R subprocess, lambda=0
     baseline plus a recorded platform/SSc-skin setting) -> which domain is
     the fibrotic/interface domain.
  5. CORE coupling test, per (substate, regulon):
       (a) bivariate Moran's I (esda.Moran_BV) with squidpy graph;
       (b) squidpy neighborhood-enrichment z-score on a discrete niche split;
       (c) coordinate-permutation null (999 perms) -> empirical p (specificity vs geometry);
       (d) size-matched RANDOM-regulon null (999 draws) -> specificity vs any regulon.
  6. Interface / CellDegree score (heterotypic KNN-6 fraction, module 505 logic).
  7. Notch-ligand SOURCING: JAG1/JAG2/DLL1/DLL4 vs HES1-regulon-high spots; which
     compartment (endothelial / myeloid / autocrine-fibroblast) is proximal source.

REUSE (no hand-written analysis from scratch — real APIs only, read from source):
  * decoupler dc.mt.ulm(tmin=5)                 <- regulon_activity.py / spatial.py
  * squidpy sq.gr.spatial_neighbors / nhood_enrichment / co_occurrence
                                                <- module 543 (verified squidpy 1.8.2 API)
  * esda.Moran / Moran_BV / Moran_Local_BV      <- spatial_enhanced/01,03,05
  * niche-class median-split + regulon-by-class <- spatial_enhanced/05_immune_niche.py
  * macrophage sub-state panels                 <- figA_macrophage_regulon_coupling.py
  * RCTD (spacexr) + CellDegree                 <- module 505 (via R subprocess)
  * BANKSY (lambda=0 baseline)                  <- module 541 (via R subprocess)
  * multi-format slice loading (Visium/Stereo-seq mtx/Xenium csv)
                                                <- spatial_enhanced/03_extra_cohorts.py

REPRODUCIBLE: fixed seeds; checkpoint + idempotent (skip slice if its outputs exist);
  ALL paths from the CFG block (no hardcoding); memory-friendly (backed/sparse reads,
  per-slice processing, light pickles only); progress logged with flush.
  squidpy uses multiprocessing -> heavy code under `if __name__ == '__main__':`.

ENV: npc (scanpy / squidpy / decoupler / esda / libpysal). RCTD + BANKSY steps shell
  out to R (CFG['r_env'] / CFG['rscript']); if R or the package is missing, the step
  degrades gracefully (marker-proxy for RCTD, expression-Leiden for domains) and the
  coupling core still runs (it does not depend on RCTD/BANKSY).

STAGE-BY-STAGE DEBUG: run with --slices <name,name> to test one/few slices first,
  --skip-rctd / --skip-banksy to isolate the pure-Python coupling core, and
  --nperm 99 for a fast smoke test before the full 999-perm run.
======================================================================================
"""
import os
import sys
import json
import glob
import pickle
import argparse
import warnings
import subprocess
import traceback
import zlib

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

import scanpy as sc
import anndata as ad
import decoupler as dc
import squidpy as sq
from scipy.stats import spearmanr
from libpysal.weights import WSP
import esda


def _stable_seed(base, *keys):
    """Deterministic per-test seed for esda.Moran_BV (F11). esda.Moran_BV has NO seed
    arg and draws permutations from numpy's GLOBAL RNG, so p_sim depends on loop/slice
    order unless we re-pin np.random immediately before each call. Python's built-in
    hash() is PYTHONHASHSEED-salted (NOT reproducible run-to-run), so we use a STABLE
    digest (zlib.crc32) of the test-identity string instead."""
    key = '|'.join(str(k) for k in keys)
    return (int(base) + zlib.crc32(key.encode('utf-8'))) % (2**31)

# ======================================================================================
# CONFIG  (edit ONLY here on the server; nothing below hardcodes a path)
# ======================================================================================
CFG = {
    # --- inputs (assumed already downloaded to the server; scripts DO NOT download) ---
    # env-var overrides match the sibling powered scripts (SSC_RAW / SSC_NET / SSC_OUT).
    'raw_root':   os.environ.get('SSC_RAW', '/data/ssc/raw'),   # root of all downloaded raw data
    'net_csv':    os.environ.get('SSC_NET', '/data/ssc/data/collectri_net.csv'),  # SAME net as discovery
    # single-cell reference for RCTD (fibroblast + immune atlas subset; barcode,celltype)
    'ref_h5ad':   os.environ.get('SSC_REF', '/data/ssc/raw/atlas/fibro_subtype_ref.h5ad'),  # ★v2: subtype-level (Fib_Myofibroblast) so RCTD/SPOTlight give a MYOFIBROBLAST fraction
    'ref_celltype_key': 'celltype',                      # obs column with the RCTD labels

    # --- spatial slices manifest ---------------------------------------------------
    # Each entry: name -> dict(platform, ...paths). Platforms handled:
    #   'visium'    : a 10x spaceranger 'outs' dir readable by sc.read_visium
    #   'visium_mtx': MatrixMarket export (expr.mtx gene x spot + genes/spots/meta) [GSE334710 style]
    #   'stereoseq' : Stereo-seq bin-aggregated .h5ad (obsm['spatial'] present) OR gem->h5ad
    #   'xenium'    : Xenium single-cell; a per-cell h5ad (obsm['spatial']) OR the 10x xenium dir
    # Fill real paths on the server. Names below mirror the MAX-DATA plan.
    'slices': {
        # ---- Visium ----
        'GSE249279_A': dict(platform='visium',     path='{raw}/spatial/GSE249279/sA/A'),
        'GSE249279_B': dict(platform='visium',     path='{raw}/spatial/GSE249279/sB/B'),
        'GSE249279_C': dict(platform='visium',     path='{raw}/spatial/GSE249279/sC/C'),
        'GSE249279_D': dict(platform='visium',     path='{raw}/spatial/GSE249279/sD/D'),
        'GSE334710':   dict(platform='visium_mtx', mtx='{raw}/spatial/GSE334710/expr.mtx',
                            genes='{raw}/spatial/GSE334710/genes.txt',
                            spots='{raw}/spatial/GSE334710/spots.txt',
                            meta='{raw}/spatial/GSE334710/meta_fixed.csv', sample_key='sample'),
        # ---- Zenodo 14577696: 10 SSc + 4 HC Visium + 4 Stereo-seq ----
        'Zenodo_visium':   dict(platform='visium_h5ad_multi',
                                h5ad='{raw}/spatial/zenodo_14577696/visium_all.h5ad',
                                sample_key='sample'),
        'Zenodo_stereoseq': dict(platform='stereoseq_multi',
                                 h5ad='{raw}/spatial/zenodo_14577696/stereoseq_all.h5ad',
                                 sample_key='sample'),
        # ---- Xenium ----
        'GSE312932_Xenium': dict(platform='xenium',
                                 h5ad='{raw}/spatial/GSE312932/xenium_all.h5ad',
                                 sample_key='sample'),
    },

    # --- outputs -------------------------------------------------------------------
    'out_dir':    os.environ.get('SSC_SPATIAL_OUT', '/data/ssc/powered/spatial'),

    # --- R subprocess (RCTD module 505, BANKSY module 541) -------------------------
    'rscript':    os.environ.get('RSCRIPT', '/data/ssc/miniconda3/envs/r44/bin/Rscript'),  # ★r44 abs path (spacexr+SPOTlight live there; npc has no R)
    'r_env':      None,                 # e.g. 'spatialR'; if set, wraps rscript in `conda run`
    'rctd_R':     None,                 # auto -> alongside this file: rctd_deconvolve.R
    'banksy_R':   None,                 # auto -> alongside this file: banksy_domains.R
    'rctd_max_cores': 12,   # RCTD parallelizes cleanly across spots (32C host; leave cores for concurrent 03/04/07)

    # --- knobs ---------------------------------------------------------------------
    'n_neighs':   6,        # KNN spatial graph
    'nperm':      999,      # coordinate-perm + random-regulon nulls (bivariate Moran uses this too)
    'min_spots':  200,      # skip a slice with fewer usable spots
    'min_counts': 200,      # per-spot count filter (Visium); Xenium single-cell -> lower inside
    # BANKSY separates platform geometry from SSc-skin biological calibration.
    # Visium retains the platform recommendation (0.2). High-resolution skin
    # uses the published SSc setting (lambda=0.25, mean k=15, gradient k=30,
    # first 20 PCs; doi:10.1016/j.ard.2025.06.002). lambda=0 is always computed.
    'banksy_lambda_visium':  float(os.environ.get('SSC_BANKSY_LAMBDA_VISIUM', '0.2')),
    'banksy_lambda_imaging': float(os.environ.get('SSC_BANKSY_LAMBDA_IMAGING', '0.25')),
    'banksy_k_geom':         int(os.environ.get('SSC_BANKSY_K_GEOM', '15')),
    'banksy_npcs':           int(os.environ.get('SSC_BANKSY_NPCS', '20')),
    'banksy_lambda': 0.2,   # back-compat default (unused by run_banksy; kept for safety)
    'seed':       0,

    # --- DOCUMENTED tuning knobs (D1: promoted from buried magic numbers; NOT headline
    #     thresholds. Named + defaulted so every choice is explicit/reproducible). ------
    'celldegree_dom_floor_sd': 0.25,  # cell_degree: min max-z (SD) to call a spot myofib/macro-dominant, else 'other'
    'niche_split_z':          0.0,    # nhood niche split: z> this = regulon/substate 'hi' (mean split at z=0)
    'randreg_size_lo_frac':   0.6,    # random-regulon size window lower = this * min(lead regulon size)
    'randreg_size_hi_frac':   1.5,    # random-regulon size window upper = this * max(lead regulon size)
}

# gene panels ---------------------------------------------------------------------------
# ★2026-07-09 MAINLINE PIVOT: the immune arm is dropped; 06 now tests whether HES1 regulon
# activity spatially CO-LOCALISES with the MYOFIBROBLAST / fibrotic program (not macrophage
# niches). The myofibroblast signature is DATA-DRIVEN from THIS atlas: the top-DE genes of the
# two activated pro-fibrotic subtypes (Myofibroblast leiden-7 + SFRP4_proFib leiden-11, argmax
# marker score in fibro_subtype_marker_scores.csv), noise-filtered, which RECOVERS the canonical
# published core (POSTN/COL11A1/CTHRC1/collagens; Tabib 2021 PMC8289865 SFRP4+/PRSS23+, Gur 2022
# PMC7612792 CTHRC1) and ADDS SSc-dermal specifics (CTGF/THY1/ASPN/TNN/EDNRA). ACTA2/TAGLN are
# deliberately DROPPED (not top-DE in this dermal atlas -> would not fit our data).
MYO = ['POSTN', 'COL11A1', 'CTHRC1', 'ASPN', 'TNN', 'COL1A1', 'COL1A2', 'COL3A1',
       'COL5A2', 'PRSS23', 'THY1', 'CTGF', 'CCN2', 'SPARC', 'LUM', 'ELN', 'MFAP5',
       'FBN1', 'BGN', 'MMP2', 'EDNRA', 'F2R']
MAC = ['CD68', 'CD163', 'LYZ', 'AIF1', 'C1QA', 'C1QB', 'MRC1', 'MSR1', 'FCGR3A', 'ITGAX', 'CD14', 'MS4A7']
# COUPLING TARGET(S): mainline = HES1 activity vs the myofibroblast program (the disease front
# that in-silico HES1-KO reverses). Kept as a dict named MAC_SUB so every downstream reference
# (score_slice / coupling_core / nhood_by_substate / figures / output) uses it unchanged. RCTD /
# CARD myofibroblast-fraction targets are appended at runtime (multi-method consensus) when
# deconvolution succeeds.
MAC_SUB = {
    'myofib_program': MYO,      # data-driven signature score (score_genes)
}
LEADS = ['HES1', 'FOSB']                 # engine-nominated leads
CTRL_POS = 'SMAD3'                       # canonical fibrotic-program TF (positive spatial ref)
CTRL_NEG = 'MEF2C'                       # negative control (activity anti-correlated w/ program)
TFS = LEADS + [CTRL_POS, CTRL_NEG]
# Notch ligands + candidate source compartments (for ligand sourcing)
NOTCH_LIG = ['JAG1', 'JAG2', 'DLL1', 'DLL4']
COMPARTMENT = {
    'endothelial': ['PECAM1', 'VWF', 'CDH5', 'CLDN5'],
    'myeloid':     ['CD68', 'LYZ', 'AIF1', 'ITGAX'],
    'fibroblast':  ['DCN', 'LUM', 'PDGFRA', 'COL1A1'],   # autocrine
}
HERE = os.path.dirname(os.path.abspath(__file__))


def P(*a):
    print(*a, flush=True)


# ======================================================================================
# I/O helpers
# ======================================================================================
def _fill(pathlike):
    """Expand '{raw}' in a manifest path against CFG['raw_root']."""
    if pathlike is None:
        return None
    return pathlike.replace('{raw}', CFG['raw_root'])


def _rscript_cmd():
    base = CFG['rscript']
    if CFG['r_env']:
        return ['conda', 'run', '-n', CFG['r_env'], 'Rscript']
    # allow "conda run -n env Rscript" as a single string too
    return base.split() if isinstance(base, str) else base


def load_net():
    """CollecTRI net (source,target,weight). Same file the whole project uses."""
    net = pd.read_csv(CFG['net_csv'])
    keep = [c for c in ['source', 'target', 'weight'] if c in net.columns]
    return net[keep]


# ---- slice loaders (multi-platform; format handling per DECISION POINTS) --------------
def _basic_qc(a, min_counts):
    a.var_names_make_unique()
    if 'counts' not in a.layers:
        a.layers['counts'] = a.X.copy()          # keep raw counts for RCTD export
    sc.pp.filter_cells(a, min_counts=min_counts)
    sc.pp.filter_genes(a, min_cells=3)
    # normalise a COPY view for scoring/ULM; RCTD uses raw counts from layers['counts']
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def load_slice(name, spec):
    """Return a dict {sub_name: AnnData} — a manifest entry may hold several sub-slices.
    Coordinates land in obsm['spatial']; raw counts kept in layers['counts'].
    coord_type for the spatial graph: 'grid' for genuine Visium hex arrays, 'generic'
    (KNN) for mtx exports / Stereo-seq bins / Xenium cells (no regular lattice)."""
    plat = spec['platform']
    out = {}

    if plat == 'visium':
        a = sc.read_visium(_fill(spec['path']))
        a = _basic_qc(a, CFG['min_counts'])
        a.uns['coord_type'] = 'grid'
        out[name] = a

    elif plat == 'visium_mtx':
        # GSE334710-style export (spatial_enhanced/03_extra_cohorts.py contract)
        expr = sio.mmread(_fill(spec['mtx'])).tocsr().T                # spots x genes
        genes = [l.strip() for l in open(_fill(spec['genes']))]
        spots = [l.strip() for l in open(_fill(spec['spots']))]
        meta = pd.read_csv(_fill(spec['meta'])).set_index('spot').loc[spots].reset_index()
        A = ad.AnnData(X=expr, obs=meta.set_index('spot'), var=pd.DataFrame(index=genes))
        A.obsm['spatial'] = meta[['x', 'y']].values
        skey = spec.get('sample_key', 'sample')
        subs = meta[skey].unique() if skey in meta.columns else [name]
        for sl in subs:
            sub = A[A.obs[skey] == sl].copy() if skey in A.obs.columns else A.copy()
            sub = _basic_qc(sub, CFG['min_counts']); sub.uns['coord_type'] = 'generic'
            out[f'{name}/{sl}'] = sub

    elif plat in ('visium_h5ad_multi', 'stereoseq_multi'):
        # Zenodo 14577696 pre-assembled multi-sample .h5ad (obsm['spatial'] present).
        # Stereo-seq: same access pattern once binned to an AnnData; bins are irregular
        # -> coord_type='generic' (KNN graph), NOT 'grid'.
        A = sc.read_h5ad(_fill(spec['h5ad']))
        if plat == 'stereoseq_multi':
            A = _ensembl_to_symbol(A)          # ★Ensembl IDs -> symbols (else HES1/MYO absent)
        skey = spec.get('sample_key', 'sample')
        ctype = 'grid' if plat == 'visium_h5ad_multi' else 'generic'
        subs = A.obs[skey].unique() if skey in A.obs.columns else [name]
        mc = CFG['min_counts'] if plat == 'visium_h5ad_multi' else max(50, CFG['min_counts'] // 4)
        for sl in subs:
            sub = A[A.obs[skey] == sl].copy() if skey in A.obs.columns else A.copy()
            if 'spatial' not in sub.obsm and {'x', 'y'} <= set(sub.obs.columns):
                sub.obsm['spatial'] = sub.obs[['x', 'y']].values
            sub = _basic_qc(sub, mc); sub.uns['coord_type'] = ctype
            out[f'{name}/{sl}'] = sub

    elif plat == 'xenium':
        # Xenium = single-cell resolution (imaging-based, few-hundred-plex panel).
        # Prefer a pre-built per-cell .h5ad (obsm['spatial']); irregular -> 'generic'.
        A = sc.read_h5ad(_fill(spec['h5ad']))
        skey = spec.get('sample_key', 'sample')
        subs = A.obs[skey].unique() if skey in A.obs.columns else [name]
        for sl in subs:
            sub = A[A.obs[skey] == sl].copy() if skey in A.obs.columns else A.copy()
            if 'spatial' not in sub.obsm and {'x', 'y'} <= set(sub.obs.columns):
                sub.obsm['spatial'] = sub.obs[['x', 'y']].values
            # single-cell: lower count floor; panel is targeted so many genes absent
            sub = _basic_qc(sub, max(10, CFG['min_counts'] // 20))
            sub.uns['coord_type'] = 'generic'; sub.uns['is_xenium'] = True
            out[f'{name}/{sl}'] = sub

    else:
        raise ValueError(f'unknown platform {plat!r} for slice {name}')

    for k, a in out.items():
        a.uns['slice_name'] = k
        a.uns['platform'] = plat        # A8/A9: platform drives RCTD doublet_mode + BANKSY lambda
    return out


# ======================================================================================
# scoring: regulon activity + macrophage sub-states + compartments
# ======================================================================================
def score_slice(a, net):
    """decoupler ULM regulon activity (TFS) + macrophage/myofib/substate/compartment
    score_genes, all per spot. Missing genes are simply skipped (Xenium panels are small)."""
    present = set(a.var_names)

    def score(genes, name):
        gg = [g for g in genes if g in present]
        if len(gg) >= 1:
            sc.tl.score_genes(a, gg, score_name=name)
        else:
            a.obs[name] = np.nan

    score(MYO, 'myofib')
    score(MAC, 'macro')
    for nm, gs in MAC_SUB.items():
        score(gs, nm)
    for nm, gs in COMPARTMENT.items():
        score(gs, f'comp_{nm}')

    # ULM regulon activity (skip if too few TF targets are in the panel, e.g. Xenium)
    net_here = net[net.target.isin(present)]
    n_tf_ok = net_here.groupby('source').size()
    have = [tf for tf in TFS if n_tf_ok.get(tf, 0) >= 5]
    acts = None
    if have:
        dc.mt.ulm(data=a, net=net_here, tmin=5)
        acts = a.obsm['score_ulm']
    # ★2026-07-09: per-TF activity with an EXPRESSION-MODE fallback for small panels (Xenium:
    # HES1's CollecTRI regulon has <5 targets in the ~5k-plex panel -> ULM cannot run, but HES1
    # itself IS on the panel). Then the TF's own log-norm expression stands in for regulon
    # activity; recorded in obs[f'{tf}_mode'] so Visium (regulon) vs Xenium (expression) stays
    # transparent and is never silently mixed in the pooled coupling.
    for tf in TFS:
        if acts is not None and tf in acts.columns:
            a.obs[tf] = acts[tf].values
            a.obs[f'{tf}_mode'] = 'regulon_ulm'
        elif tf in present:
            a.obs[tf] = np.asarray(_col(a, tf)).ravel()
            a.obs[f'{tf}_mode'] = 'expression'
        else:
            a.obs[tf] = np.nan
            a.obs[f'{tf}_mode'] = 'absent'
    if acts is None:
        a.uns['ulm_skipped'] = True
    # Notch ligand expression per spot (for sourcing)
    for lg in NOTCH_LIG:
        if lg in present:
            a.obs[f'lig_{lg}'] = np.asarray(_col(a, lg)).ravel()
        else:
            a.obs[f'lig_{lg}'] = np.nan
    return a


def _col(a, gene):
    x = a[:, gene].X
    return x.toarray() if sp.issparse(x) else np.asarray(x)


def _ensembl_to_symbol(A):
    """★2026-07-09: Stereo-seq (Zenodo) exports Ensembl gene IDs (ENSG...); map them to HGNC
    symbols via mygene so score_genes / the CollecTRI ULM find HES1/SMAD3/MYO by symbol. Drops
    unmapped IDs, collapses duplicate symbols (var_names_make_unique). No-op if already symbols."""
    if not any(str(g).startswith('ENSG') for g in list(A.var_names[:50])):
        return A
    import mygene
    ids = [str(g).split('.')[0] for g in A.var_names]           # strip version suffix
    res = mygene.MyGeneInfo().querymany(ids, scopes='ensembl.gene', fields='symbol',
                                        species='human', as_dataframe=True, verbose=False)
    sym = ({} if 'symbol' not in res.columns
           else res[~res.index.duplicated(keep='first')]['symbol'].to_dict())
    new = [sym.get(i) for i in ids]
    keep = [j for j, s in enumerate(new) if isinstance(s, str) and s]
    A = A[:, keep].copy()
    A.var_names = [new[j] for j in keep]
    A.var_names_make_unique()
    P('    [ensembl->symbol] mapped %d/%d genes (Stereo-seq)' % (len(keep), len(ids)))
    return A


# ======================================================================================
# spatial graph + weights (reused from spatial_enhanced/01,03,05)
# ======================================================================================
def build_graph(a):
    ct = a.uns.get('coord_type', 'generic')
    sq.gr.spatial_neighbors(a, coord_type=ct, n_neighs=CFG['n_neighs'])
    w = WSP(a.obsp['spatial_connectivities']).to_W()
    w.transform = 'r'
    return w


def _couple_targets(a):
    """★v2 multi-method: coupling targets = the myofibroblast PROGRAM score (MAC_SUB, data-driven
    signature) PLUS any deconvolution-derived myofibroblast-FRACTION obs columns that are present
    and non-constant ('myofib_RCTD', 'myofib_SPOTlight'). So the pooled coupling reports HES1 vs
    up to THREE independent myofibroblast-abundance measures (signature + RCTD + SPOTlight)."""
    t = list(MAC_SUB.keys())
    for c in ('myofib_RCTD', 'myofib_SPOTlight'):
        if c in a.obs.columns:
            v = np.asarray(a.obs[c], float)
            if np.all(np.isfinite(v)) and np.nanstd(v) > 0:
                t.append(c)
    return t


# ======================================================================================
# CORE coupling test: bivariate Moran + nhood enrichment + TWO nulls
# ======================================================================================
def coupling_core(a, w, net, rng):
    """For each (macrophage sub-state, lead-regulon) pair, quantify spatial coupling with
    (a) bivariate Moran's I, (b) coordinate-permutation empirical p (999),
    (c) size-matched random-regulon null empirical p (specificity control).
    Returns a list[dict]. Honest: coupling = spatial association, not causation."""
    nperm = CFG['nperm']
    coord = np.asarray(a.obsm['spatial'], float)
    rows = []

    # snapshot the REAL lead/ctrl TF activities BEFORE building the random-regulon
    # pool (the pool's dc.mt.ulm overwrites obsm['score_ulm']); use the snapshot so
    # coupling never depends on obs being restored correctly.
    tf_act = {tf: np.asarray(a.obs[tf], float) for tf in TFS}

    # pre-compute a pool of size-matched random regulons for the specificity null.
    # A "random regulon" = ULM activity of a randomly-drawn TF with a target-set size
    # matched to the lead. This asks: is the coupling SPECIFIC to HES1/FOSB, or would
    # any regulon of similar size couple to the macrophage state just as well?
    rand_pool = _random_regulon_pool(a, net, rng, tf_act, n_draws=min(60, nperm))

    # ---- valid (substate, TF) pairs: finite + non-constant on both axes ----
    valid_pairs, y_snap = [], {}
    for sub in _couple_targets(a):
        y = np.asarray(a.obs[sub], float)
        if not np.all(np.isfinite(y)) or np.nanstd(y) == 0:
            continue
        y_snap[sub] = y
        for tf in TFS:
            x = tf_act[tf]
            if not np.all(np.isfinite(x)) or np.nanstd(x) == 0:
                continue
            valid_pairs.append((sub, tf))

    # (b) coordinate-permutation null -- build the permuted spatial graph ONCE PER PERMUTATION
    #     and REUSE it across ALL (sub,tf) pairs (the scrambled geometry is pair-independent).
    #     The old code rebuilt the KNN graph per pair per permutation (~n_pairs*nperm ~= 16k
    #     builds) -> a multi-hour/day hang on 100K-cell Xenium/Stereo slices. This keeps the
    #     FULL nperm (NO resolution loss) at exactly nperm graph builds. Sharing one permutation
    #     set across pairs leaves each pair's null nperm valid draws, so every individual
    #     perm_p_coord stays statistically valid (only cross-pair independence changes).
    null_geom = {p: np.empty(nperm) for p in valid_pairs}
    for i in range(nperm):
        order = rng.permutation(a.n_obs)
        a.obsm['spatial'] = coord[order]
        sq.gr.spatial_neighbors(a, coord_type=a.uns.get('coord_type', 'generic'),
                                n_neighs=CFG['n_neighs'])
        wp = WSP(a.obsp['spatial_connectivities']).to_W(); wp.transform = 'r'
        for (sub, tf) in valid_pairs:
            null_geom[(sub, tf)][i] = esda.Moran_BV(tf_act[tf], y_snap[sub], wp, permutations=0).I
    a.obsm['spatial'] = coord                                        # restore real geometry
    sq.gr.spatial_neighbors(a, coord_type=a.uns.get('coord_type', 'generic'),
                            n_neighs=CFG['n_neighs'])

    # per-pair reporting: observed bivariate Moran + spearman + p_geom + regulon-specificity null
    for (sub, tf) in valid_pairs:
        y, x = y_snap[sub], tf_act[tf]
        # esda.Moran_BV has no seed arg + uses the numpy GLOBAL RNG; re-pin with a STABLE
        # per-test seed so p_sim is reproducible independent of slice/loop order.
        np.random.seed(_stable_seed(CFG['seed'], a.uns.get('slice_name'), sub, tf))
        bv = esda.Moran_BV(x, y, w, permutations=nperm)              # value-permutation p
        rho, prho = spearmanr(x, y)
        p_geom = (np.sum(null_geom[(sub, tf)] >= bv.I) + 1) / (nperm + 1)

        # (c) size-matched random-regulon null (specificity): is the lead in the tail?
        if tf in LEADS and rand_pool:
            null_reg = np.array([esda.Moran_BV(rx, y, w, permutations=0).I
                                 for rx in rand_pool])
            p_reg = (np.sum(null_reg >= bv.I) + 1) / (len(null_reg) + 1)
        else:
            p_reg = np.nan

        rows.append(dict(
            substate=sub, regulon=tf,
            moranBV=float(bv.I), moranBV_p=float(bv.p_sim),
            spearman=float(rho), spearman_p=float(prho),
            perm_p_coord=float(p_geom),        # geometry-specificity
            perm_p_randreg=float(p_reg),       # regulon-specificity
            kind=('LEAD' if tf in LEADS else ('pos_ctrl' if tf == CTRL_POS else 'neg_ctrl')),
        ))
    return rows


def _random_regulon_pool(a, net, rng, tf_act, n_draws=60):
    """ULM activity for n_draws random TFs whose target-set size (in the current panel)
    is matched to the lead regulons. Returns list of arrays. Restores the real TF
    activities in a.obs from the tf_act snapshot (dc.mt.ulm here overwrote score_ulm)."""
    if a.uns.get('ulm_skipped'):
        return []
    present = set(a.var_names)
    net_here = net[net.target.isin(present)]
    sizes = net_here.groupby('source').size()
    lead_sizes = [sizes.get(tf, 0) for tf in LEADS if sizes.get(tf, 0) >= 5]
    if not lead_sizes:
        return []
    # size-match window for the specificity null (D1: documented CFG knobs, was 0.6/1.5)
    lo = CFG['randreg_size_lo_frac'] * min(lead_sizes)
    hi = CFG['randreg_size_hi_frac'] * max(lead_sizes)
    cand = sizes[(sizes >= max(5, lo)) & (sizes <= hi)].index.tolist()
    cand = [c for c in cand if c not in TFS]
    if len(cand) < 5:
        cand = sizes[sizes >= 5].index.tolist()
    if not cand:
        return []
    pick = list(rng.choice(cand, size=min(n_draws, len(cand)), replace=False))
    sub_net = net_here[net_here.source.isin(pick)]
    dc.mt.ulm(data=a, net=sub_net, tmin=5)
    acts = a.obsm['score_ulm']
    pool = []
    for tf in pick:
        if tf in acts.columns:
            v = np.asarray(acts[tf].values, float)
            if np.all(np.isfinite(v)) and np.nanstd(v) > 0:
                pool.append(v)
    # dc.mt.ulm overwrote obsm['score_ulm'] with the random regulons; restore the
    # real HES1/FOSB/SMAD3/MEF2C activities in a.obs from the snapshot (downstream
    # nhood_by_substate / ligand_sourcing read a.obs[tf]).
    for tf in TFS:
        a.obs[tf] = tf_act[tf]
    return pool


# ======================================================================================
# niche split + neighborhood enrichment  (reuse 05_immune_niche.py logic, per substate)
# ======================================================================================
def nhood_by_substate(a, sub):
    """Median-split myofib(=HES1/FOSB front proxy via 'macro' co-context) niche classes as in
    05_immune_niche.py, but here we build a 4-class split of (regulon-high x substate-high)
    for neighborhood enrichment: does the macrophage sub-state sit ADJACENT to the
    regulon-high myofibroblast front? Returns dict per lead regulon."""
    res = {}
    z = lambda v: (v - np.nanmean(v)) / (np.nanstd(v) + 1e-9)
    ss = z(np.asarray(a.obs[sub], float))
    # niche split threshold (D1: documented CFG knob; default 0.0 = mean-split at z=0)
    _zsplit = CFG['niche_split_z']
    for tf in LEADS:
        x = np.asarray(a.obs[tf], float)
        if not np.all(np.isfinite(x)) or np.nanstd(x) == 0:
            continue
        zr = z(x)
        rh, sh = zr > _zsplit, ss > _zsplit
        cls = np.where(rh & sh, 'Interface',
                       np.where(rh & ~sh, f'{tf}-hi',
                                np.where(~rh & sh, f'{sub}-hi', 'Other')))
        cats = ['Interface', f'{tf}-hi', f'{sub}-hi', 'Other']
        a.obs['_ncls'] = pd.Categorical(cls, categories=cats)
        # every niche class must be populated or nhood_enrichment's z-matrix indexing
        # is undefined; assert instead of swallowing a degenerate split.
        vc = a.obs['_ncls'].value_counts()
        if (vc.reindex(cats).fillna(0) == 0).any():
            # A degenerate split (an empty niche class -> nhood z-matrix indexing undefined)
            # is exactly what a strongly-coupled or small slice hits. SKIP this (sub,tf) pair
            # (surfaced via print, NOT silently) instead of raising, so the slice stays
            # otherwise complete and its checkpoint consistent rather than aborting the slice.
            print(f'    [nhood] ({sub},{tf}) degenerate split {vc.to_dict()} -> skipped', flush=True)
            continue
        sq.gr.nhood_enrichment(a, cluster_key='_ncls', seed=CFG['seed'],
                               n_perms=CFG['nperm'], show_progress_bar=False)
        zmat = a.uns['_ncls_nhood_enrichment']['zscore']
        i_r = cats.index(f'{tf}-hi'); i_s = cats.index(f'{sub}-hi')
        res[tf] = dict(z_core_adj=float(zmat[i_r, i_s]), cats=cats,
                       zmat=zmat.tolist(),
                       comp={c: float((a.obs['_ncls'] == c).mean()) for c in cats})
    if '_ncls' in a.obs:
        del a.obs['_ncls']
    return res


# ======================================================================================
# CellDegree interface score (module 505 base KNN-6 logic)
# ======================================================================================
def cell_degree(a):
    """Heterotypic-neighbour fraction per spot on a simple compartment call (myofib vs
    macro vs other by argmax of z-scores). High = mixing / interface zone. module 505."""
    zc = {}
    for k in ['myofib', 'macro']:
        v = np.asarray(a.obs[k], float)
        zc[k] = (v - np.nanmean(v)) / (np.nanstd(v) + 1e-9)
    zmat = np.vstack([zc['myofib'], zc['macro']]).T
    # NOTE (D1: documented judgment-call knob, not a published cutoff): a spot is
    # 'other' (neither myofib- nor macro-dominant) when its max compartment z-score
    # is below this SD floor. Named CFG knob so the choice is explicit/reproducible.
    _dom_floor = CFG['celldegree_dom_floor_sd']
    dom = np.where(zmat.max(1) < _dom_floor, 'other',
                   np.array(['myofib', 'macro'])[zmat.argmax(1)])
    # use the already-built KNN graph connectivities to find neighbours
    conn = a.obsp['spatial_connectivities'].tocsr()
    het = np.zeros(a.n_obs)
    for i in range(a.n_obs):
        nn = conn.indices[conn.indptr[i]:conn.indptr[i + 1]]
        if len(nn):
            het[i] = np.mean(dom[nn] != dom[i])
    a.obs['celldegree'] = het
    a.obs['compartment'] = dom
    return float(np.mean(het))


# ======================================================================================
# Notch-ligand sourcing: which compartment sources the ligand near HES1-hi spots
# ======================================================================================
def ligand_sourcing(a, w):
    """Bivariate spatial coupling (esda.Moran_BV) of each Notch ligand's expression vs
    the CONTINUOUS HES1 regulon activity across ALL spots (no self-created high/low
    binarisation), plus ligand<->compartment bivariate Moran to indicate the co-
    localising compartment. Honest: proximity/association, NOT a proven axis.
    NOTE (F09): bivariate Moran is computed on the CONTINUOUS activity, so there is NO
    'top tertile' threshold; the docstring is corrected to match the code.
    Returns per-slice ligand + compartment sourcing rows."""
    rows = []
    if 'HES1' not in a.obs or not np.isfinite(np.asarray(a.obs['HES1'], float)).all():
        return rows
    hes = np.asarray(a.obs['HES1'], float)
    if np.nanstd(hes) == 0:
        return rows
    for lg in NOTCH_LIG:
        lv = np.asarray(a.obs[f'lig_{lg}'], float)
        if not np.all(np.isfinite(lv)) or np.nanstd(lv) == 0:
            continue
        # F11: stable per-test re-seed before the permutation-bearing Moran_BV.
        np.random.seed(_stable_seed(CFG['seed'], a.uns.get('slice_name'), 'lig', lg))
        bv = esda.Moran_BV(lv, hes, w, permutations=CFG['nperm'])   # ligand <-> HES1 activity
        # which compartment co-localises with THIS ligand (bivariate Moran ligand<->comp)
        comp_bv = {}
        for nm in COMPARTMENT:
            cv = np.asarray(a.obs[f'comp_{nm}'], float)
            if np.all(np.isfinite(cv)) and np.nanstd(cv) > 0:
                comp_bv[nm] = float(esda.Moran_BV(lv, cv, w, permutations=0).I)
        # ★2026-07-08 audit: do NOT hard-argmax three bivariate Moran I's into one confident
        # "source" — require a margin between top-1 and top-2, else call it 'ambiguous'. The full
        # comp_bv vector is reported alongside so the raw co-localisation I's stay inspectable.
        if comp_bv:
            _rk = sorted(comp_bv.items(), key=lambda kv: kv[1], reverse=True)
            _v2 = _rk[1][1] if len(_rk) > 1 else float('-inf')
            src = _rk[0][0] if (_rk[0][1] - _v2) >= float(CFG.get('ligand_source_margin', 0.05)) else 'ambiguous'
        else:
            src = 'NA'
        rows.append(dict(ligand=lg, moranBV_vs_HES1=float(bv.I), p_sim=float(bv.p_sim),
                         source_compartment=src, **{f'bv_{k}': v for k, v in comp_bv.items()}))
    return rows


# ======================================================================================
# RCTD (module 505 via R subprocess) + fallback
# ======================================================================================
def run_rctd(a, sl_key, net):
    """Export slice counts, call rctd_deconvolve.R (spacexr). On any failure -> a
    marker-score proxy (normalised myofib/macro/substate scores) so downstream fractions
    still exist. Returns a DataFrame [spots x fractions]. Idempotent via out csv."""
    out_csv = os.path.join(CFG['out_dir'], 'rctd', f'{_safe(sl_key)}_weights.csv')
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    if os.path.exists(out_csv) and os.path.getsize(out_csv) > 10:
        P(f'    [RCTD] {os.path.basename(out_csv)} exists -> load')
        return pd.read_csv(out_csv, index_col=0)

    rctd_R = CFG['rctd_R'] or os.path.join(HERE, 'rctd_deconvolve.R')
    ref_ok = os.path.exists(_fill(CFG['ref_h5ad']))
    if ref_ok and os.path.exists(rctd_R):
        try:
            tmp = os.path.join(CFG['out_dir'], 'rctd', '_tmp', _safe(sl_key))
            os.makedirs(tmp, exist_ok=True)
            # export spatial slice raw counts (natural spot x gene; _export_mm
            # transposes to the gene x spot MatrixMarket that RCTD expects) + coords
            cnt = a.layers['counts']
            cnt = cnt if sp.issparse(cnt) else sp.csr_matrix(cnt)
            _export_mm(cnt, a.var_names, a.obs_names, os.path.join(tmp, 'sp'))
            pd.DataFrame(a.obsm['spatial'], index=a.obs_names, columns=['x', 'y']).to_csv(
                os.path.join(tmp, 'sp_coords.csv'), index_label='barcode')
            # export reference once (cached across slices)
            ref_dir = os.path.join(CFG['out_dir'], 'rctd', '_ref')
            if not os.path.exists(os.path.join(ref_dir, 'ref.mtx')):
                _export_reference(ref_dir)
            # A8: doublet_mode per platform (spacexr contract). 'doublet' = high-res
            # 1-2 cells/pixel (Xenium single-cell, Stereo-seq cell-bins); 'full' = low-res
            # multi-cell spots (Visium 55um). ★2026-07-07 修: 原映射反了(Visium 曾→doublet)。
            plat = a.uns.get('platform', '')
            dmode = 'doublet' if plat in ('stereoseq_multi', 'xenium') else 'full'
            cmd = _rscript_cmd() + [rctd_R,
                os.path.join(ref_dir, 'ref.mtx'), os.path.join(ref_dir, 'ref_genes.txt'),
                os.path.join(ref_dir, 'ref_barcodes.txt'), os.path.join(ref_dir, 'ref_labels.csv'),
                os.path.join(tmp, 'sp.mtx'), os.path.join(tmp, 'sp_genes.txt'),
                os.path.join(tmp, 'sp_barcodes.txt'), os.path.join(tmp, 'sp_coords.csv'),
                out_csv, str(CFG['rctd_max_cores']), dmode]
            P(f'    [RCTD] {sl_key}: doublet_mode={dmode} (platform={plat}); '
              f'{" ".join(cmd[:3])} ... (R subprocess)')
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            if r.returncode == 0 and os.path.exists(out_csv):
                P(f'    [RCTD] OK -> {os.path.basename(out_csv)}')
                return pd.read_csv(out_csv, index_col=0)
            raise RuntimeError(
                f'[RCTD] {sl_key}: spacexr deconvolution FAILED rc={r.returncode}. '
                f'stderr tail:\n{r.stderr[-1500:]}\n'
                f'No marker-score proxy is substituted (would produce silently-wrong '
                f'"fractions"). Fix the R/spacexr env or pass --skip-rctd to run the '
                f'coupling core without deconvolution.')
        except subprocess.TimeoutExpired:
            raise
        except Exception as e:
            raise RuntimeError(f'[RCTD] {sl_key}: subprocess raised -> {e}') from e
    else:
        # ref or R helper genuinely absent: HARD FAIL, do not fabricate fractions.
        raise FileNotFoundError(
            f'[RCTD] reference h5ad ({_fill(CFG["ref_h5ad"])}, exists={ref_ok}) or '
            f'rctd_deconvolve.R ({rctd_R}, exists={os.path.exists(rctd_R)}) missing. '
            f'RCTD cannot run. Provide them or pass --skip-rctd (coupling core does not '
            f'depend on RCTD).')


def run_spotlight(a, sl_key):
    """★v2 2nd deconvolution: SPOTlight, REUSING the MatrixMarket ref + slice files that run_rctd
    already exported (so no re-export). Returns [spots x celltype] fractions (incl Fib_Myofibroblast)
    or None on any failure (never blocks; the coupling then uses signature + RCTD only). Idempotent."""
    spot_R = os.path.join(HERE, 'spotlight_deconvolve.R')
    ref_dir = os.path.join(CFG['out_dir'], 'rctd', '_ref')
    tmp = os.path.join(CFG['out_dir'], 'rctd', '_tmp', _safe(sl_key))
    out_csv = os.path.join(CFG['out_dir'], 'spotlight', f'{_safe(sl_key)}_spotlight.csv')
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    if os.path.exists(out_csv) and os.path.getsize(out_csv) > 10:
        return pd.read_csv(out_csv, index_col=0)
    need = ([os.path.join(ref_dir, x) for x in ('ref.mtx', 'ref_genes.txt', 'ref_barcodes.txt', 'ref_labels.csv')]
            + [os.path.join(tmp, x) for x in ('sp.mtx', 'sp_genes.txt', 'sp_barcodes.txt', 'sp_coords.csv')])
    if not (os.path.exists(spot_R) and all(os.path.exists(x) for x in need)):
        P(f'    [SPOTlight] {sl_key}: R script or RCTD-export files missing -> skip (uses RCTD+signature)')
        return None
    cmd = _rscript_cmd() + [spot_R] + need + [out_csv, str(CFG['rctd_max_cores'])]
    P(f'    [SPOTlight] {sl_key}: R subprocess (reusing RCTD export) ...')
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except subprocess.TimeoutExpired:
        P(f'    [SPOTlight] {sl_key}: TIMEOUT -> skip'); return None
    if r.returncode == 0 and os.path.exists(out_csv):
        P(f'    [SPOTlight] OK -> {os.path.basename(out_csv)}')
        return pd.read_csv(out_csv, index_col=0)
    P(f'    [SPOTlight] {sl_key}: FAILED rc={r.returncode}. stderr tail:\n{r.stderr[-800:]}')
    return None


def _export_mm(mat_sp_spotxgene, genes, barcodes, prefix):
    """Write MatrixMarket (gene x spot as R expects) + genes/barcodes txt.
    mat_sp_spotxgene is spot x gene -> transpose to gene x spot for R."""
    gxs = mat_sp_spotxgene.T.tocoo()
    sio.mmwrite(prefix + '.mtx', gxs)
    with open(prefix + '_genes.txt', 'w') as f:
        f.write('\n'.join(map(str, genes)))
    with open(prefix + '_barcodes.txt', 'w') as f:
        f.write('\n'.join(map(str, barcodes)))


def _export_reference(ref_dir):
    """Export the single-cell reference (raw counts + labels) once for RCTD."""
    os.makedirs(ref_dir, exist_ok=True)
    R = sc.read_h5ad(_fill(CFG['ref_h5ad']))
    R.var_names_make_unique()
    key = CFG['ref_celltype_key']
    X = R.layers['counts'] if 'counts' in R.layers else R.X
    X = X if sp.issparse(X) else sp.csr_matrix(X)
    _export_mm(X, R.var_names, R.obs_names, os.path.join(ref_dir, 'ref'))
    # rename ref.mtx (prefix helper writes ref.mtx already) -> ensure names
    if not os.path.exists(os.path.join(ref_dir, 'ref.mtx')):
        os.replace(os.path.join(ref_dir, 'ref') + '.mtx', os.path.join(ref_dir, 'ref.mtx'))
    pd.DataFrame({'barcode': R.obs_names, 'celltype': R.obs[key].astype(str).values}).to_csv(
        os.path.join(ref_dir, 'ref_labels.csv'), index=False)
    P(f'    [RCTD] reference exported ({R.n_obs} cells, {R.obs[key].nunique()} types)')


# ======================================================================================
# BANKSY domains (module 541 via R subprocess) + fallback
# ======================================================================================
def run_banksy(a, sl_key):
    """BANKSY domains with a lambda=0 baseline and recorded spatial parameters.
    Fallback: expression-Leiden domain (lambda=0 equivalent). Returns per-spot domain
    label Series. Idempotent."""
    plat = a.uns.get('platform', '')
    lam = CFG['banksy_lambda_imaging'] if plat in ('stereoseq_multi', 'xenium') \
        else CFG['banksy_lambda_visium']
    k_geom = CFG['banksy_k_geom']
    npcs = CFG['banksy_npcs']
    param_tag = f"lambda{lam:g}_k{k_geom}_pc{npcs}".replace('.', 'p')
    out_csv = os.path.join(
        CFG['out_dir'], 'banksy', f'{_safe(sl_key)}_domains_{param_tag}.csv')
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    if os.path.exists(out_csv) and os.path.getsize(out_csv) > 10:
        d = pd.read_csv(out_csv)
        return d.set_index(d.columns[0]).iloc[:, -1]

    banksy_R = CFG['banksy_R'] or os.path.join(HERE, 'banksy_domains.R')
    if os.path.exists(banksy_R):
        try:
            tmp = os.path.join(
                CFG['out_dir'], 'banksy', '_tmp', _safe(sl_key), param_tag)
            os.makedirs(tmp, exist_ok=True)
            long_csv = os.path.join(tmp, 'slice_long.csv')
            _export_banksy_long(a, long_csv)
            cmd = _rscript_cmd() + [banksy_R, '--input', long_csv,
                                    '--outdir', tmp, '--lambda', str(lam),
                                    '--k_geom', str(k_geom), '--npcs', str(npcs)]
            P(f'    [BANKSY] {sl_key}: lambda={lam}, k={k_geom}/{2*k_geom}, '
              f'PCs={npcs} (platform={plat}); R subprocess ...')
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            assign = os.path.join(tmp, 'spot_domain_assignments.csv')
            if r.returncode == 0 and os.path.exists(assign):
                dd = pd.read_csv(assign)
                dd.to_csv(out_csv, index=False)
                col = 'spatial' if 'spatial' in dd.columns else dd.columns[-1]
                return dd.set_index(dd.columns[0])[col]
            raise RuntimeError(
                f'[BANKSY] {sl_key}: FAILED rc={r.returncode}. stderr tail:\n'
                f'{r.stderr[-1500:]}\nNo expression-Leiden proxy is substituted '
                f'(a non-spatial Leiden clustering is NOT a BANKSY spatial domain). '
                f'Fix the R/Banksy env or pass --skip-banksy.')
        except subprocess.TimeoutExpired:
            raise
        except Exception as e:
            raise RuntimeError(f'[BANKSY] {sl_key}: subprocess raised -> {e}') from e
    else:
        raise FileNotFoundError(
            f'[BANKSY] banksy_domains.R ({banksy_R}) not found. Provide it or pass '
            f'--skip-banksy (coupling core does not depend on BANKSY domains).')


def _export_banksy_long(a, path):
    """Long-format CSV for module 541: spot,x,y,<HVGs...>. Use top-HVG raw counts to
    keep the file small (memory-friendly). module 541 auto-handles no 'domain' col."""
    b = a.copy()
    # HVG selection must succeed deterministically; do not silently fall back to
    # exporting all genes (that would change what BANKSY sees between runs/envs).
    sc.pp.highly_variable_genes(b, n_top_genes=min(400, b.n_vars), flavor='seurat')
    b = b[:, b.var['highly_variable']].copy()
    X = b.layers['counts'] if 'counts' in b.layers else b.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    df = pd.DataFrame(X, columns=list(b.var_names), index=b.obs_names)
    df.insert(0, 'y', a.obsm['spatial'][:, 1])
    df.insert(0, 'x', a.obsm['spatial'][:, 0])
    df.insert(0, 'spot', b.obs_names)
    df.to_csv(path, index=False)


# ======================================================================================
# figures (NO plain bars — dot / heatmap / lollipop / spatial-scatter)
# ======================================================================================
def make_figs(pooled_coupling, pooled_ligand, light):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig_dir = os.path.join(CFG['out_dir'], 'figs')
    os.makedirs(fig_dir, exist_ok=True)

    # --- Fig B1: coupling dotplot (substate x regulon), color=mean BV, size=-log10 p ---
    if len(pooled_coupling):
        cp = pooled_coupling.groupby(['substate', 'regulon']).agg(
            bv=('moranBV', 'mean'), p=('perm_p_coord', 'median')).reset_index()
        # ★v2: substates come from the ACTUAL coupling data (myofib_program + myofib_RCTD +
        # myofib_SPOTlight), NOT the MAC_SUB dict, so multi-method targets are all plotted.
        subs = sorted(cp['substate'].unique())
        regs = [t for t in TFS if t in set(cp['regulon'])] or list(TFS)
        fig, axd = plt.subplots(figsize=(6.4, 0.6 * len(subs) + 2.6))
        for _, r in cp.iterrows():
            if r['regulon'] not in regs or r['substate'] not in subs:
                continue
            xi = regs.index(r['regulon']); yi = subs.index(r['substate'])
            size = 30 + 220 * min(3.0, -np.log10(max(r['p'], 1e-3))) / 3.0
            axd.scatter(xi, yi, s=size, c=[r['bv']], cmap='RdBu_r',
                        vmin=-0.3, vmax=0.3, edgecolor='k', linewidth=0.4)
        axd.set_xticks(range(len(regs))); axd.set_xticklabels(regs, rotation=30, ha='right')
        axd.set_yticks(range(len(subs))); axd.set_yticklabels(subs)
        axd.set_title('Spatial coupling: HES1/FOSB activity × myofibroblast measure\n'
                      '(signature / RCTD / SPOTlight; colour = mean bivariate Moran I; '
                      'size = -log10 coord-perm p)')
        sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=plt.Normalize(-0.3, 0.3))
        fig.colorbar(sm, ax=axd, fraction=0.046, pad=0.04, label="bivariate Moran's I")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, 'FigB1_coupling_dotplot.pdf'))
        fig.savefig(os.path.join(fig_dir, 'FigB1_coupling_dotplot.png'), dpi=150)
        plt.close(fig)

    # --- Fig B2: per-slice heatmap of LEAD coupling (substate x slice), HES1 & FOSB ---
    if len(pooled_coupling):
        for tf in LEADS:
            sub = pooled_coupling[pooled_coupling.regulon == tf]
            if not len(sub):
                continue
            piv = sub.pivot_table(index='substate', columns='slice', values='moranBV')
            fig, axh = plt.subplots(figsize=(max(6, 0.5 * piv.shape[1] + 3), 3.2))
            im = axh.imshow(piv.values, cmap='RdBu_r', vmin=-0.3, vmax=0.3, aspect='auto')
            axh.set_xticks(range(piv.shape[1]))
            axh.set_xticklabels(piv.columns, rotation=45, ha='right', fontsize=6)
            axh.set_yticks(range(piv.shape[0])); axh.set_yticklabels(piv.index)
            axh.set_title(f'{tf} regulon × macrophage sub-state coupling (per slice)')
            fig.colorbar(im, ax=axh, fraction=0.046, pad=0.04, label="bivariate Moran's I")
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, f'FigB2_{tf}_coupling_heatmap.pdf'))
            fig.savefig(os.path.join(fig_dir, f'FigB2_{tf}_coupling_heatmap.png'), dpi=150)
            plt.close(fig)

    # --- Fig B3: Notch-ligand sourcing lollipop (ligand vs HES1, colored by source) ---
    if len(pooled_ligand):
        lg = pooled_ligand.groupby('ligand').agg(bv=('moranBV_vs_HES1', 'mean')).reset_index()
        src = pooled_ligand.groupby('ligand')['source_compartment'].agg(
            lambda s: s.value_counts().index[0])
        lg['source'] = lg['ligand'].map(src)
        lg = lg.sort_values('bv')
        cmap = {'endothelial': '#2166ac', 'myeloid': '#b2182b', 'fibroblast': '#1a9850', 'NA': '#888'}
        fig, axl = plt.subplots(figsize=(5.6, 3.2))
        yp = np.arange(len(lg))
        axl.hlines(yp, 0, lg['bv'], color='#bbb', lw=1.2)
        axl.scatter(lg['bv'], yp, s=90,
                    c=[cmap.get(s, '#888') for s in lg['source']], edgecolor='k', zorder=3)
        axl.set_yticks(yp); axl.set_yticklabels(lg['ligand'])
        axl.axvline(0, color='k', lw=0.6)
        axl.set_xlabel("bivariate Moran's I  (ligand × HES1 activity)")
        axl.set_title('Notch-ligand sourcing near HES1-high spots\n(dot color = dominant source compartment)')
        from matplotlib.lines import Line2D
        axl.legend(handles=[Line2D([0], [0], marker='o', color='w', markerfacecolor=c,
                                   markersize=8, label=k) for k, c in cmap.items() if k != 'NA'],
                   fontsize=7, frameon=False, loc='lower right')
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, 'FigB3_notch_ligand_sourcing.pdf'))
        fig.savefig(os.path.join(fig_dir, 'FigB3_notch_ligand_sourcing.png'), dpi=150)
        plt.close(fig)

    # --- Fig B4: spatial scatter of one representative slice (HES1 act + top substate) ---
    if light:
        k0 = next(iter(light))
        d = light[k0]
        coord = d['spatial']; obs = d['obs']
        panels = [('HES1', 'HES1 regulon activity'),
                  ('SPP1_hi', 'SPP1-hi macrophage'),
                  ('TREM2_hi_SAM', 'TREM2-hi/SAM (PROTECTIVE)'),
                  ('celldegree', 'Interface (CellDegree)')]
        panels = [(k, t) for k, t in panels if k in obs.columns]
        fig, axs = plt.subplots(1, len(panels), figsize=(4.3 * len(panels), 4.2))
        if len(panels) == 1:
            axs = [axs]
        for axp, (k, t) in zip(axs, panels):
            sc2 = axp.scatter(coord[:, 0], coord[:, 1], c=np.asarray(obs[k], float),
                              cmap='magma', s=8, edgecolor='none')
            axp.set_title(t, fontsize=9); axp.set_aspect('equal'); axp.axis('off')
            fig.colorbar(sc2, ax=axp, fraction=0.046, pad=0.04)
        fig.suptitle(f'Representative slice: {k0}', fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, 'FigB4_spatial_scatter.pdf'))
        fig.savefig(os.path.join(fig_dir, 'FigB4_spatial_scatter.png'), dpi=150)
        plt.close(fig)

    P(f'    [figs] -> {fig_dir}')


# ======================================================================================
# small utils
# ======================================================================================
def _safe(s):
    return s.replace('/', '__').replace(' ', '_')


def _slice_done(sl_key):
    """Idempotency check: a slice is done if its per-slice coupling csv exists."""
    return os.path.exists(os.path.join(CFG['out_dir'], 'per_slice',
                                       f'{_safe(sl_key)}_coupling.csv'))


# ======================================================================================
# per-slice driver
# ======================================================================================
def process_slice(sl_key, a, net, args):
    rng = np.random.default_rng(CFG['seed'])
    ps_dir = os.path.join(CFG['out_dir'], 'per_slice')
    os.makedirs(ps_dir, exist_ok=True)

    if a.n_obs < CFG['min_spots']:
        P(f'  [{sl_key}] only {a.n_obs} spots (< {CFG["min_spots"]}) -> SKIP')
        return None

    P(f'  [{sl_key}] scoring (ULM regulon + substates + compartments) ...')
    a = score_slice(a, net)
    if a.uns.get('ulm_skipped'):
        P(f'  [{sl_key}] NOTE: ULM skipped (too few TF targets in panel, e.g. Xenium)')

    P(f'  [{sl_key}] spatial graph (coord_type={a.uns.get("coord_type")}) ...')
    w = build_graph(a)

    # RCTD + SPOTlight + BANKSY (optional; do not block the coupling core)
    frac = None
    if not args.skip_rctd:
        frac = run_rctd(a, sl_key, net)
        if frac is not None:
            frac = frac.reindex(a.obs_names)
            for c in frac.columns:
                if c not in ('_proxy',):
                    a.obs[f'rctd_{c}'] = frac[c].values
            # ★v2 multi-method: the RCTD Fib_proFibrotic (myofibroblast-front) fraction is a
            # coupling target. Collapsed class (Myofibroblast+SFRP4_proFib) -> deconvolution
            # separates it cleanly (fine subtypes are transcriptionally inseparable on Visium).
            if 'Fib_proFibrotic' in frac.columns:
                a.obs['myofib_RCTD'] = frac['Fib_proFibrotic'].values
        # SPOTlight = 2nd deconvolution (reuses the RCTD-exported ref+slice mtx); never blocks
        if not getattr(args, 'skip_spotlight', False):
            sfrac = run_spotlight(a, sl_key)
            if sfrac is not None and 'Fib_proFibrotic' in sfrac.columns:
                a.obs['myofib_SPOTlight'] = sfrac.reindex(a.obs_names)['Fib_proFibrotic'].values
    domain = None
    if not args.skip_banksy:
        domain = run_banksy(a, sl_key)
        if domain is not None:
            a.obs['banksy_domain'] = domain.reindex(a.obs_names).astype(str).values

    P(f'  [{sl_key}] CellDegree interface ...')
    mean_het = cell_degree(a)

    P(f'  [{sl_key}] CORE coupling (bivariate Moran + coord-perm + rand-regulon nulls) ...')
    cp_rows = coupling_core(a, w, net, rng)
    for r in cp_rows:
        r['slice'] = sl_key
    cp = pd.DataFrame(cp_rows)
    # NOTE: coupling.csv (the _slice_done() checkpoint key) is written at the END of this
    # function, only AFTER nhood/ligand/light all succeed (see below), so a mid-slice failure
    # never leaves a silently-partial 'done' checkpoint.

    P(f'  [{sl_key}] neighborhood enrichment per substate ...')
    nhood = {sub: nhood_by_substate(a, sub) for sub in MAC_SUB
             if np.isfinite(np.asarray(a.obs[sub], float)).all()
             and np.nanstd(np.asarray(a.obs[sub], float)) > 0}
    with open(os.path.join(ps_dir, f'{_safe(sl_key)}_nhood.pkl'), 'wb') as fh:
        pickle.dump(nhood, fh)

    P(f'  [{sl_key}] Notch-ligand sourcing ...')
    lg_rows = ligand_sourcing(a, w)
    for r in lg_rows:
        r['slice'] = sl_key
    lg = pd.DataFrame(lg_rows)
    lg.to_csv(os.path.join(ps_dir, f'{_safe(sl_key)}_ligand.csv'), index=False)

    # light bundle for figures (small — no full matrices)
    keep_obs = ['myofib', 'macro', 'celldegree', 'compartment'] + list(MAC_SUB) + TFS + \
               [f'lig_{lg_}' for lg_ in NOTCH_LIG]
    keep_obs = [c for c in keep_obs if c in a.obs.columns]
    light = dict(spatial=np.asarray(a.obsm['spatial'], float),
                 obs=a.obs[keep_obs].copy(), mean_heterotypic=mean_het,
                 coord_type=a.uns.get('coord_type'))
    with open(os.path.join(ps_dir, f'{_safe(sl_key)}_light.pkl'), 'wb') as fh:
        pickle.dump(light, fh)

    # coupling.csv is the _slice_done() checkpoint key -> write it LAST, only after
    # nhood.pkl/ligand.csv/light.pkl all succeeded, so a mid-slice failure leaves the slice
    # correctly marked INCOMPLETE (re-run next time) rather than a silently-partial 'done'.
    cp.to_csv(os.path.join(ps_dir, f'{_safe(sl_key)}_coupling.csv'), index=False)

    # cp can be an empty/column-less DataFrame when a slice has no finite regulon activity
    # (e.g. a Xenium panel where ULM was skipped) -> guard before accessing .kind
    if len(cp) and 'kind' in cp.columns:
        n_lead_sig = int(((cp.kind == 'LEAD') & (cp.perm_p_coord < 0.05)).sum())
        n_lead = int((cp.kind == 'LEAD').sum())
    else:
        n_lead_sig, n_lead = 0, 0
        P(f'  [{sl_key}] no finite regulon activity -> coupling table empty')
    P(f'  [{sl_key}] DONE: {a.n_obs} spots | lead couplings p<0.05 (coord-perm): '
      f'{n_lead_sig}/{n_lead} | mean CellDegree={mean_het:.3f}')
    return dict(coupling=cp, ligand=lg, light=light)


# ======================================================================================
# main
# ======================================================================================
def main():
    ap = argparse.ArgumentParser(description='POWERED spatial immune-stromal coupling (Fig B)')
    ap.add_argument('--slices', type=str, default=None,
                    help='comma-separated slice manifest keys to run (default: all)')
    ap.add_argument('--nperm', type=int, default=None, help='override CFG nperm (smoke test)')
    ap.add_argument('--skip-rctd', action='store_true', help='skip RCTD (isolate coupling core)')
    ap.add_argument('--skip-spotlight', action='store_true', help='skip SPOTlight 2nd deconvolution')
    ap.add_argument('--skip-banksy', action='store_true', help='skip BANKSY')
    ap.add_argument('--force', action='store_true', help='recompute even if per-slice csv exists')
    args = ap.parse_args()
    if args.nperm:
        CFG['nperm'] = args.nperm

    np.random.seed(CFG['seed'])
    os.makedirs(CFG['out_dir'], exist_ok=True)
    P(f'=== POWERED spatial immune coupling | out={CFG["out_dir"]} | nperm={CFG["nperm"]} ===')

    # persist the resolved config for reproducibility
    with open(os.path.join(CFG['out_dir'], 'run_config.json'), 'w') as fh:
        json.dump({k: v for k, v in CFG.items()}, fh, indent=2, default=str)

    net = load_net()
    P(f'CollecTRI net: {net.shape[0]} edges, {net.source.nunique()} TFs')

    manifest = CFG['slices']
    if args.slices:
        want = set(args.slices.split(','))
        manifest = {k: v for k, v in manifest.items() if k in want}
        P(f'restricted to slices: {list(manifest)}')

    all_cp, all_lg, light_all = [], [], {}
    for name, spec in manifest.items():
        P(f'\n[MANIFEST] {name} ({spec["platform"]})')
        try:
            slices = load_slice(name, spec)
        except FileNotFoundError as e:
            P(f'  [{name}] input missing -> SKIP ({e})')
            continue
        except Exception as e:
            P(f'  [{name}] load FAILED -> SKIP: {e}')
            traceback.print_exc()
            continue

        for sl_key, a in slices.items():
            if _slice_done(sl_key) and not args.force:
                P(f'  [{sl_key}] already done -> load per-slice outputs')
                try:
                    cp = pd.read_csv(os.path.join(CFG['out_dir'], 'per_slice',
                                                  f'{_safe(sl_key)}_coupling.csv'))
                    all_cp.append(cp)
                    lgf = os.path.join(CFG['out_dir'], 'per_slice', f'{_safe(sl_key)}_ligand.csv')
                    if os.path.exists(lgf) and os.path.getsize(lgf) > 5:
                        all_lg.append(pd.read_csv(lgf))
                    lpf = os.path.join(CFG['out_dir'], 'per_slice', f'{_safe(sl_key)}_light.pkl')
                    if os.path.exists(lpf) and sl_key not in light_all:
                        with open(lpf, 'rb') as fh:
                            light_all[sl_key] = pickle.load(fh)
                except Exception as e:
                    P(f'    reload failed ({e}); recompute')
                else:
                    continue
            try:
                res = process_slice(sl_key, a, net, args)
            except Exception as e:
                P(f'  [{sl_key}] PROCESS FAILED: {e}')
                traceback.print_exc()
                continue
            if res is not None:
                all_cp.append(res['coupling'])
                if len(res['ligand']):
                    all_lg.append(res['ligand'])
                light_all[sl_key] = res['light']
            del a  # free memory before next slice

    # ---- pool + summarise ----
    if not all_cp:
        P('\nNo slices produced coupling results. Check CFG paths / --slices. Exiting.')
        return
    pooled_cp = pd.concat(all_cp, ignore_index=True)
    pooled_lg = pd.concat(all_lg, ignore_index=True) if all_lg else pd.DataFrame()
    pooled_cp.to_csv(os.path.join(CFG['out_dir'], 'pooled_coupling.csv'), index=False)
    if len(pooled_lg):
        pooled_lg.to_csv(os.path.join(CFG['out_dir'], 'pooled_ligand_sourcing.csv'), index=False)

    # pooled lead summary: mean BV + fraction of slices significant (coord-perm & rand-reg)
    lead = pooled_cp[pooled_cp.regulon.isin(LEADS)]
    from statsmodels.stats.multitest import multipletests
    # BH-FDR across ALL coupling tests in the pooled table (per p-value family).
    for pcol, qcol in [('moranBV_p', 'moranBV_q'), ('perm_p_coord', 'perm_p_coord_q'),
                       ('perm_p_randreg', 'perm_p_randreg_q'), ('spearman_p', 'spearman_q')]:
        pv = pd.to_numeric(pooled_cp[pcol], errors='coerce')
        m = pv.notna()
        if m.sum() >= 1:   # ★2026-07-08 BUG-A: an all-NaN p-family (e.g. perm_p_randreg on Xenium/targeted panels) -> empty array -> multipletests crashes (毁掉整个汇总/FDR/图阶段)
            pooled_cp.loc[m, qcol] = multipletests(pv[m].values, method='fdr_bh')[1]
        else:
            pooled_cp[qcol] = np.nan
    # if a whole p-family was all-NaN the q-column above was never created -> ensure it exists
    for _qc in ['moranBV_q', 'perm_p_coord_q', 'perm_p_randreg_q', 'spearman_q']:
        if _qc not in pooled_cp.columns:
            pooled_cp[_qc] = np.nan
    pooled_cp.to_csv(os.path.join(CFG['out_dir'], 'pooled_coupling.csv'), index=False)
    lead = pooled_cp[pooled_cp.regulon.isin(LEADS)]
    summ = lead.groupby(['substate', 'regulon']).agg(
        n_slices=('slice', 'nunique'),
        mean_moranBV=('moranBV', 'mean'),
        frac_sig_coordperm=('perm_p_coord_q', lambda s: float((pd.to_numeric(s, errors='coerce') < 0.05).mean())),
        frac_sig_randreg=('perm_p_randreg_q', lambda s: float((pd.to_numeric(s, errors='coerce') < 0.05).mean())),
        mean_spearman=('spearman', 'mean')).reset_index()
    summ.to_csv(os.path.join(CFG['out_dir'], 'pooled_lead_summary.csv'), index=False)

    P('\n=== POOLED lead coupling (macrophage sub-state x HES1/FOSB) ===')
    P(summ.round(3).to_string(index=False))
    if len(pooled_lg):
        P('\n=== POOLED Notch-ligand sourcing (mean BV vs HES1; dominant source) ===')
        from statsmodels.stats.multitest import multipletests
        _pv = pd.to_numeric(pooled_lg['p_sim'], errors='coerce')
        _m = _pv.notna()
        if _m.sum() >= 1:   # guard (mirror pooled_cp): an all-NaN p_sim family -> empty array -> multipletests crashes
            pooled_lg.loc[_m, 'p_sim_q'] = multipletests(_pv[_m].values, method='fdr_bh')[1]
        else:
            pooled_lg['p_sim_q'] = np.nan
        srcsumm = pooled_lg.groupby('ligand').agg(
            mean_bv=('moranBV_vs_HES1', 'mean'),
            frac_sig=('p_sim_q', lambda s: float((pd.to_numeric(s, errors='coerce') < 0.05).mean())),
            top_source=('source_compartment', lambda s: s.value_counts().index[0])).reset_index()
        P(srcsumm.round(3).to_string(index=False))

    P('\n[figs] rendering (dot / heatmap / lollipop / spatial-scatter; no bars) ...')
    make_figs(pooled_cp, pooled_lg, light_all)

    # version snapshot (reproducibility iron-law 6)
    import importlib.metadata as ilm
    with open(os.path.join(CFG['out_dir'], 'versions.txt'), 'w') as fh:
        for p in ['scanpy', 'squidpy', 'decoupler', 'anndata', 'esda', 'libpysal',
                  'numpy', 'pandas', 'scipy']:
            try:
                fh.write(f'{p}=={ilm.version(p)}\n')
            except Exception:
                fh.write(f'{p}==NA\n')

    P('\nSPATIAL_IMMUNE_COUPLING_DONE')


if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()          # squidpy nhood_enrichment uses multiprocessing
    main()
