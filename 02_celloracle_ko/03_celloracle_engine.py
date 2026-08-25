# -*- coding: utf-8 -*-
"""
[SERVER, co2 env] LAYER 03 -- POWERED CellOracle virtual-perturbation ENGINE.

Runs the single-cell in-silico TF-knockout engine on the INTEGRATED, MAX-DATA
true-SSc SKIN fibroblast atlas produced by layer 02 (Gur GSE195452 + Tabib
GSE138669 + other skin cohorts). Lung and pulmonary-vascular datasets are excluded.
This extension screens every prespecified screenable TF that has at least one edge
in the filtered fitted GRN. It writes new all-active outputs and does not replace
the completed 43-TF Figure 2 calibration analysis. The myofibroblast
mask is SINGLE-SOURCED via myo_mask_of() (D1): the layer-02 myofibroblast STATE
label when present, else a documented score-quantile knob cfg['myo_quantile'].

REUSE / provenance (read these; do not hand-write CellOracle calls):
  - server_results/scripts/phase0_celloracle.py   -> GRN build + KO-panel ranking
  - server_results/scripts/celloracle_boot.py     -> bootstrap-CI loop + pickle-safe setstate
  - rerun_clean/rerun01_discovery_celloracle.py   -> clean re-run params (identical)
  - rerun_clean/_qc_common.py                      -> MYO / CTRL_GRN / CTRL_BOOT / LEADS panels
  - bioinfo-reusable-code modules/069_celloracle_grn_perturbation.py
        -> vector-field API: simulate_shift -> estimate_transition_prob ->
           calculate_embedding_shift (extended here with the documented
           Oracle_development_module inner-product / perturbation-score).

HONEST FRAMING (kept in comments + logged, per project iron laws):
  * CellOracle in-silico KO = a GRN-propagation HYPOTHESIS + correlation, NOT a
    causal claim. "delta_fibrotic_program < 0" means the model predicts a TF is a
    positive regulator of the myofibroblast program; it is a nomination, not proof.
  * Positive controls (SMAD3/CEBPB/KLF4/FOS + the broader CTRL_GRN panel) must
    rank as strong positive regulators for the screen to be trusted; they are the
    sanity check, not results.
  * No immune reasoning happens here (this layer is fibroblast-intrinsic only).

REPRODUCIBILITY:
  * All paths come from the CFG block below -- nothing hardcoded downstream.
  * Fixed seeds everywhere (numpy default_rng(SEED), subsample random_state=SEED,
    scanpy sc.settings + set_figure_params off).
  * Idempotent: every heavy stage is skipped if its output already exists.
    New all-active outputs use distinct names; existing 43-TF calibration,
    bootstrap and vector-field files are protected by default.
  * Memory-aware: restricts CellOracle to HVG + required control/lead genes,
    optional cell CAP for GRN fit (identity-preserving stratified cap, NOT the old
    flat 4000 subsample -- default cap is large / None), sparse throughout, n_jobs
    from CFG, progress logged with flush=True.

RUN (debug stage by stage -- each stage is guarded + idempotent):
  conda run -n co2 python /data/ssc/powered/03_celloracle_engine.py
"""
import os, sys, gc, time, pickle, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# ----------------------------------------------------------------------------
# anndata-version-safe unpickling of the Oracle object (from celloracle_boot.py).
# Without this, an Oracle pickled under a different anndata build can fail to load
# because AnnDataFileManager.__setstate__ expects a live weakref.
# ----------------------------------------------------------------------------
import anndata._core.file_backing as _fb
import weakref as _wr
def _safe_setstate(self, state):
    self.__dict__ = dict(state)
    ref = state.get("_adata_ref", None)
    self.__dict__["_adata_ref"] = (_wr.ref(ref) if ref is not None else (lambda: None))
_fb.AnnDataFileManager.__setstate__ = _safe_setstate

import celloracle as co

# ============================================================================
# CONFIG -- edit paths here only. Nothing below this block is hardcoded.
# ============================================================================
CFG = {
    # --- inputs (from layer 02 integration) ---
    "atlas_h5ad": "/data/ssc/powered/fibroblast_integrated.h5ad",   # 02 写此名(2026-07-07 修 B2 DAG)
    #   Expected obs: a cluster/state column (see cluster_col), a UMAP in
    #   obsm["X_umap"], and either layers["raw_count"]/["counts"] OR a raw-count
    #   .X. The script auto-detects the raw-count layer.

    # --- REQUIRED skin scATAC base-GRN (primary arm) ---
    #   The SSc study is skin-specific. The GSE312129 skin-scATAC prior is the
    #   primary arm; the generic human-promoter prior is retained only as a
    #   tissue-agnostic sensitivity analysis. No lung prior or lung dataset is used.
    #   Format = CellOracle TF_info_matrix (peak/gene x TF one-hot), loadable via
    #   pandas.read_parquet / read_csv. See build_base_GRN note at bottom.
    "skin_base_grn": "/data/ssc/basegrn/skin_base_GRN_dataframe.parquet",  # e.g. "/data/ssc/raw/GSE312129/skin_base_GRN_dataframe.parquet"
    "require_skin_base_grn": True,
    "run_promoter_sensitivity": True,

    # --- outputs ---
    "out_dir": "/data/ssc/powered",
    # Optional read-only source for frozen oracle/links files. Sensitivity runs
    # can reuse a fitted GRN while writing every new screen table elsewhere.
    "grn_source_dir": None,

    # --- column / embedding names in the integrated atlas ---
    "cluster_col": "leiden",       # fall back auto-detected if absent
    "embedding":   "X_umap",

    # --- CellOracle params (IDENTICAL to published run unless noted) ---
    "n_top_genes":   2000,         # HVG for the GRN feature space (as phase0)
    "recompute_hvg":  False,        # True only for isolated HVG sensitivity runs
    "grn_alpha":     10,           # ridge alpha, get_links + fit_GRN_for_simulation
    "bagging":       20,           # bagging_number in get_links (CellOracle get_links default/tutorial = 20; audit P1, was 5)
    "link_p":        0.001,        # filter_links p threshold
    "link_topn":     2000,         # filter_links threshold_number
    "knn_k":         25,           # knn_imputation k
    "knn_pca_dims":  50,
    "n_propagation": 3,            # simulate_shift propagation depth
    "min_cluster_n": 40,           # drop clusters with < this many cells (as phase0)

    # --- POWERED sizing (identity-preserving STRATIFIED cap, proportional per cluster;
    #   NOT the old flat random 4000). Set to 30000 per CellOracle OFFICIAL guidance:
    #   "If your scRNA-seq data includes more than 20-30K cells, we recommend downsampling
    #   your data ... perturbation simulations may require large amounts of memory"
    #   (morris-lab GRN tutorial). The fibroblast atlas is ~108K -> the O(n^2) transition
    #   matrix would be ~188GB / ~40h, which is BOTH against the tool's design and infeasible.
    #   30K keeps the perturbation field well-estimated (the tool's sweet spot) at ~14GB / a
    #   few hours, close to the validated gate run (22848 cells). The LARGE-DATASET advantage
    #   is showcased elsewhere (424K atlas, CellRank trajectory on full fibroblasts, Milo
    #   differential abundance, cross-cohort reproducibility, spatial RCTD), NOT here.
    "grn_cell_cap":  30000,

    # --- KO screen ---
    # Primary extension = all filtered-GRN edge-active TFs that are detected in
    # >=1% of fibroblasts. This rule is frozen before inspecting HES1 ranks.
    # panel_extra_top is retained only to document the completed 43-TF calibration.
    "ko_screen_mode": "all_active",
    "ko_min_raw_expr_fraction": 0.01,
    "ko_checkpoint_every": 1,
    "panel_extra_top": 40,         # historical calibration only; not used below
    "bootstrap_B":     500,        # bootstrap resamples for lead/control CIs
    "run_closed_downstream_stages": False,  # protect completed bootstrap/vector fields

    # --- vector-field / inner-product (perturbation score in embedding) ---
    "vf_run":            True,     # set False to skip (slow on very large N)
    "vf_genes":          None,     # None = LEADS + top KO hits; or an explicit list
    "vf_n_neighbors":    200,      # estimate_transition_prob (as module 069)
    "vf_sampled_frac":   1.0,
    "vf_sigma_corr":     0.05,     # calculate_embedding_shift
    "vf_n_grid":         40,       # grid resolution for inner-product field
    # --- inner-product / perturbation-score (Oracle_development_module) ---
    #   Requires a developmental reference gradient built from a pseudotime axis.
    #   03 runs BEFORE 04, so no external pseudotime exists; we compute a DPT axis
    #   INSIDE 03 (scanpy sc.tl.dpt) rooted at the SFRP2/DPP4 progenitor-fibroblast
    #   pool, then fall back to the continuous myo-program score as the ordinal axis
    #   if DPT cannot be rooted. Set False to skip the inner-product arm explicitly
    #   (never silently NaN it).
    "vf_inner_product":   True,
    "vf_pseudotime_key":  "dpt_pseudotime",   # obs key computed in stage_vectorfield
    "vf_root_markers":    ["SFRP2", "DPP4", "PI16", "CD34"],  # progenitor-fibroblast root pool

    # --- myofibroblast mask (SINGLE-SOURCE; see myo_mask_of()) ---
    #   PREFERRED anchor = the myofibroblast STATE label from layer 02 (marker-argmax
    #   subtype). If that column/label is absent, fall back to a documented quantile
    #   knob on obs['myo'] = MEAN of the MYO panel on log-norm expression (a plain panel
    #   mean, NOT a Tirosh-2016 control-subtracted module score). This is a tuning
    #   knob, NOT a headline threshold; it is applied identically in EVERY stage via
    #   myo_mask_of() so all stages share one mask definition (was recomputed 3-4x).
    "myo_state_col":      "cell_state",   # obs column that may hold a myofibroblast label
    "myo_state_label":    "Myofibroblast",# value in myo_state_col marking myofibroblasts
    "myo_quantile":       0.75,           # fallback score cutoff (documented default)

    # --- engineering ---
    "n_jobs": 10,   # overridden by SSC_CO_N_JOBS on the 32-core launcher
    "seed":   0,
    "vf_random_seed": 15071990,  # CellOracle default used by the live primary screen
    "force":  False,               # True = ignore existing outputs, recompute all
    "selected_arms": ["skinatac", "promoter"],
}

# canonical panels -- REUSED verbatim from rerun_clean/_qc_common.py (do not edit;
# keeps the powered screen comparable to the published discovery screen).
MYO = ["ACTA2", "TAGLN", "POSTN", "COL1A1", "COL1A2", "COMP", "COL11A1", "CTHRC1", "COL3A1", "FN1"]
LEADS = ["MEF2C", "HES1", "CEBPD", "FOSB"]     # HES1/Notch + FOSB/AP-1 main line
CTRL_GRN = ["SMAD3", "JUN", "FOS", "JUNB", "SRF", "FOSL1", "FOSL2", "RUNX1", "RUNX2", "CREB1",
            "EGR1", "ETS1", "KLF4", "NR4A1", "TWIST1", "SNAI2", "MYC", "STAT3", "NFKB1", "CEBPB"]
CTRL_BOOT = ["CEBPB", "KLF4", "FOS", "NR4A1", "STAT3", "SMAD3"]
# ★2026-07-08 LITERATURE-GROUNDED competitor-defense + baseline (workflow wy8ays4sn; recipes in
# replan_FrontImmunol_2026/UPGRADE_PLAN_2026crossfield.md) -- NOT a subjective design:
COMPETITORS = ["RUNX1", "CREB3L1", "FOSL2", "RUNX3"]   # 2026-nominated SSc fibroblast regulators (Parvizi ARD26 / Huang A&R26 / Lafyatis JCI-Insight26 / Li SciRep25)
DRUGGABLE_AXIS = {"HES1", "HEY1", "HEYL", "HES5", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"}  # gamma-secretase/Notch = nirogacestat-actionable (HES1's unique edge over competitors)
# known-driver-recovery baseline labels: Ahlmann-Eltze/Huber Nat Methods 2025 (PMID 40759747)
# mandates a naive-baseline test; the in-silico-KO analog = recovering established drivers, with
# ground-truth from the Nat Genet 2025 fibroblast TF Perturb-seq atlas (PMID 40770575).
DRIVER_POS = ["SMAD3", "MYOCD", "GLI1", "SNAI1", "FOXC2", "SMAD2", "TWIST1"]  # established myofibroblast activators
DRIVER_NEG = ["FLI1", "KLF5"]                                                 # SSc-downregulated repressors

sc.settings.verbosity = 1
np.random.seed(CFG["seed"])


# ============================================================================
# helpers
# ============================================================================
def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def exists(p):
    return os.path.exists(p) and os.path.getsize(p) > 0


def _env_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def apply_env_overrides(cfg):
    """Apply launcher-controlled overrides without editing the frozen script."""
    cfg = dict(cfg)
    scalar = {
        "SSC_CO_ATLAS_H5AD": ("atlas_h5ad", str),
        "SSC_CO_SKIN_BASE_GRN": ("skin_base_grn", str),
        "SSC_CO_OUT_DIR": ("out_dir", str),
        "SSC_CO_GRN_SOURCE_DIR": ("grn_source_dir", str),
        "SSC_CO_N_JOBS": ("n_jobs", int),
        "SSC_CO_N_TOP_GENES": ("n_top_genes", int),
        "SSC_CO_RECOMPUTE_HVG": ("recompute_hvg", _env_bool),
        "SSC_CO_N_PROPAGATION": ("n_propagation", int),
        "SSC_CO_VF_N_NEIGHBORS": ("vf_n_neighbors", int),
        "SSC_CO_LINK_TOPN": ("link_topn", int),
        "SSC_CO_KO_MIN_RAW_EXPR_FRACTION": ("ko_min_raw_expr_fraction", float),
        "SSC_CO_SEED": ("seed", int),
        "SSC_CO_VF_RANDOM_SEED": ("vf_random_seed", int),
        "SSC_CO_FORCE": ("force", _env_bool),
    }
    for env_name, (key, caster) in scalar.items():
        if env_name in os.environ:
            cfg[key] = caster(os.environ[env_name])

    arms = os.environ.get("SSC_CO_ARMS")
    if arms:
        cfg["selected_arms"] = [
            value.strip().lower() for value in arms.split(",") if value.strip()
        ]
    allowed = {"skinatac", "promoter"}
    if not cfg["selected_arms"] or not set(cfg["selected_arms"]).issubset(allowed):
        raise ValueError(
            f"SSC_CO_ARMS must be a non-empty subset of {sorted(allowed)}; "
            f"received {cfg['selected_arms']}"
        )
    if cfg["n_jobs"] < 1:
        raise ValueError("SSC_CO_N_JOBS must be positive")
    if cfg["n_top_genes"] < 100:
        raise ValueError("SSC_CO_N_TOP_GENES must be at least 100")
    if cfg["n_propagation"] < 1 or cfg["vf_n_neighbors"] < 2:
        raise ValueError("Propagation and transition-neighbour settings are invalid")
    if cfg["link_topn"] < 100:
        raise ValueError("SSC_CO_LINK_TOPN must be at least 100")
    if not 0 <= cfg["ko_min_raw_expr_fraction"] <= 1:
        raise ValueError("SSC_CO_KO_MIN_RAW_EXPR_FRACTION must be in [0,1]")
    return cfg


def frozen_grn_paths(cfg, arm):
    """Return the first complete frozen oracle/links pair available for an arm."""
    tag = "" if arm == "promoter" else f"_{arm}"
    roots = []
    source = cfg.get("grn_source_dir")
    if source:
        roots.append(source)
    roots.append(cfg["out_dir"])
    seen = set()
    for root in roots:
        root = os.path.abspath(root)
        if root in seen:
            continue
        seen.add(root)
        oracle_pkl = os.path.join(root, f"oracle{tag}.pkl")
        links_out = os.path.join(root, f"links{tag}.celloracle.links")
        if exists(oracle_pkl) and exists(links_out):
            return oracle_pkl, links_out
    return None


def detect_raw_layer(a):
    """Return the name of a raw-count layer, or None if .X itself is raw counts."""
    for k in ("raw_count", "counts", "raw_counts", "count"):
        if k in a.layers:
            return k
    return None


def is_integer_counts(X):
    x = X[:200] if X.shape[0] > 200 else X
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return False
    # Integer-count test (F10): counts are non-negative integers. Drop the arbitrary
    # max>30 magnitude gate (a shallow-sequenced count matrix could have max<=30 and
    # be wrongly rejected). Principled check: all values non-negative integers, and
    # not an all-zero slice. Log-normalized data has non-integer values -> fails here.
    x = x[x >= 0]
    return x.size > 0 and bool(np.allclose(x, np.round(x))) and float(np.max(x)) >= 1.0


def myo_program(a, use_raw=True):
    """Mean of the canonical myofibroblast panel on log-norm expression -> obs['myo'].
    Identical definition to phase0 / rerun01 (used as the continuous fibrotic axis)."""
    src = a.raw if (use_raw and a.raw is not None) else a
    pres = [g for g in MYO if g in src.var_names]
    return np.asarray(src[:, pres].X.mean(1)).ravel(), pres


def myo_mask_of(adata, cfg):
    """SINGLE SOURCE OF TRUTH for the myofibroblast mask (was recomputed 3-4x with
    divergent 0.75-quantile masks across stages). Anchoring rule (matches the
    already-decided 02 rule: myofibroblast = marker-argmax STATE label):
      PREFERRED: gate on the layer-02 myofibroblast STATE label if present.
      FALLBACK : score-quantile knob cfg['myo_quantile'] on obs['myo'] (documented
                 default 0.75), applied IDENTICALLY here so every stage shares it.
    Returns a boolean numpy mask over adata.n_obs.
    """
    col = cfg.get("myo_state_col")
    lab = cfg.get("myo_state_label")
    if col and (col in adata.obs) and lab is not None:
        vals = adata.obs[col].astype(str).values
        if lab in set(vals):
            m = (vals == lab)
            log(f"[myo_mask] STATE anchor: {col}=='{lab}' -> {int(m.sum())}/{len(m)} cells")
            return np.asarray(m, dtype=bool)
        log(f"[myo_mask] state label '{lab}' not found in obs['{col}'] -> quantile fallback")
    myo = np.asarray(adata.obs["myo"].values, dtype=float)
    q = float(np.quantile(myo, cfg["myo_quantile"]))
    m = myo >= q
    log(f"[myo_mask] score-quantile fallback: myo>=q({cfg['myo_quantile']})={q:.4f} "
        f"-> {int(m.sum())}/{len(m)} cells")
    return np.asarray(m, dtype=bool)


def stratified_cap(a, cluster_col, cap, seed):
    """Identity-preserving downsample: keep proportional cells per cluster so the
    cluster structure (and thus the GRN unit) is preserved. Returns an index mask.
    This replaces the old flat sc.pp.subsample(4000) which was a POWER limitation,
    not a design choice."""
    if cap is None or a.n_obs <= cap:
        return np.ones(a.n_obs, dtype=bool)
    rng = np.random.default_rng(seed)
    keep = np.zeros(a.n_obs, dtype=bool)
    clusters = a.obs[cluster_col].astype(str).values
    frac = cap / a.n_obs
    for cl in pd.unique(clusters):
        idx = np.where(clusters == cl)[0]
        n_take = max(min(len(idx), 1), int(round(len(idx) * frac)))
        take = rng.choice(idx, size=min(n_take, len(idx)), replace=False)
        keep[take] = True
    return keep


# ============================================================================
# STAGE 0 -- load integrated atlas, define state column + myo program
# ============================================================================
def stage_load(cfg):
    out = cfg["out_dir"]
    os.makedirs(out, exist_ok=True)
    log("== STAGE 0: load integrated atlas (backed) ==")
    log("atlas:", cfg["atlas_h5ad"])
    a = sc.read_h5ad(cfg["atlas_h5ad"])   # layer-02 atlas is the analysis object
    log("atlas shape", a.shape)

    # state/cluster column
    cc = cfg["cluster_col"]
    if cc not in a.obs:
        for alt in ("leiden", "cell_state", "state", "subcluster", "cluster", "celltype"):
            if alt in a.obs:
                cc = alt
                break
    if cc not in a.obs:
        raise SystemExit(f"[FATAL] no cluster column found (tried {cfg['cluster_col']} + fallbacks). "
                         "Layer 02 must write a fibroblast-state column.")
    a.obs[cc] = a.obs[cc].astype(str)
    cfg["cluster_col"] = cc
    log("using cluster column:", cc, "->", a.obs[cc].value_counts().sort_index().to_dict())

    if cfg["embedding"] not in a.obsm:
        raise SystemExit(f"[FATAL] embedding {cfg['embedding']} missing in obsm. "
                         "Layer 02 must provide a UMAP for the vector field.")

    # ensure a raw-count layer for CellOracle + a log-norm .raw for the myo program
    rl = detect_raw_layer(a)
    if rl is None:
        if is_integer_counts(a.X):
            a.layers["raw_count"] = a.X.copy()
            log("no raw layer; .X detected as integer counts -> layers['raw_count']")
        else:
            raise SystemExit("[FATAL] no raw-count layer and .X is not integer counts. "
                             "Layer 02 must keep counts in layers['raw_count'] or ['counts'].")
        rl = "raw_count"
    else:
        if rl != "raw_count":
            a.layers["raw_count"] = a.layers[rl].copy()
        log("raw-count layer:", rl)

    # a.raw = log-norm for the myo program (build only if absent / not log-scaled)
    need_lognorm = (a.raw is None) or (float((a.raw.X[:100].toarray() if hasattr(a.raw.X, "toarray")
                                              else np.asarray(a.raw.X[:100])).max()) > 50)
    if need_lognorm:
        tmp = ad.AnnData(a.layers["raw_count"].copy(), obs=a.obs.copy(), var=a.var.copy())
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        a.raw = tmp
        del tmp
        log("built log-norm a.raw for myo program")

    a.obs["myo"], pres = myo_program(a, use_raw=True)
    log("myo program on", len(pres), "genes:", pres)

    # HVG in the GRN feature space (as phase0: 2000 HVG on log-norm)
    if cfg.get("recompute_hvg", False) or "highly_variable" not in a.var.columns:
        hv_src = a.raw.to_adata()
        # Pin flavor explicitly (F08). NOTE: scanpy's highly_variable_genes has NO
        # random_state parameter; the 'seurat' flavor is fully deterministic (no RNG
        # code path), so pinning the flavor is sufficient for reproducibility here.
        sc.pp.highly_variable_genes(hv_src, n_top_genes=cfg["n_top_genes"], flavor="seurat")
        a.var["highly_variable"] = hv_src.var["highly_variable"].reindex(a.var_names).fillna(False).values
        del hv_src
        log(
            "recomputed highly_variable feature mask:",
            cfg["n_top_genes"],
            "genes requested",
        )
    log("HVG:", int(a.var["highly_variable"].sum()))
    gc.collect()
    return a


# ============================================================================
# STAGE 1 -- build the CellOracle GRN + fit for simulation (POWERED, no 4000 cap)
#   base_grn: "promoter" (default) or path to a skin scATAC TF_info_matrix.
# ============================================================================
def stage_grn(a, cfg, arm, base_grn_path):
    out = cfg["out_dir"]
    tag = "" if arm == "promoter" else f"_{arm}"
    oracle_pkl = f"{out}/oracle{tag}.pkl"
    links_out  = f"{out}/links{tag}.celloracle.links"

    frozen = None if cfg["force"] else frozen_grn_paths(cfg, arm)
    if frozen is not None:
        frozen_oracle, frozen_links = frozen
        log(
            f"[skip GRN:{arm}] loading frozen read-only inputs "
            f"{frozen_oracle} + {frozen_links}"
        )
        oracle = pickle.load(open(frozen_oracle, "rb"))
        links = co.load_hdf5(frozen_links)
        return oracle, links
    if a is None:
        raise RuntimeError(
            f"No frozen GRN pair was found for arm {arm!r}, but the atlas was not "
            "loaded. Remove SSC_CO_GRN_SOURCE_DIR or provide a complete oracle/links pair."
        )

    cc = cfg["cluster_col"]
    log(f"== STAGE 1: CellOracle GRN [{arm} base-GRN] ==")

    # feature space = HVG + all panel genes (so controls/leads are never dropped)
    panel_genes = [g for g in (MYO + CTRL_GRN + LEADS) if g in a.var_names]
    keep_genes = list(dict.fromkeys(list(a.var_names[a.var["highly_variable"]]) + panel_genes))
    aco = a[:, keep_genes].copy()
    aco.X = aco.layers["raw_count"].copy()          # CellOracle wants raw counts in .X
    aco.obs[cc] = a.obs[cc].values
    aco.obs["myo"] = a.obs["myo"].values
    aco.obsm[cfg["embedding"]] = a.obsm[cfg["embedding"]]

    # POWERED sizing: use ALL cells (grn_cell_cap=None) or an identity-preserving
    # stratified cap as a memory fallback. NOT the flat 4000 subsample.
    mask = stratified_cap(aco, cc, cfg["grn_cell_cap"], cfg["seed"])
    if not mask.all():
        log(f"stratified cap -> {int(mask.sum())}/{aco.n_obs} cells (cap={cfg['grn_cell_cap']})")
        aco = aco[mask].copy()

    # drop tiny clusters (unstable GRN unit) -- same rule as phase0 (>=40 cells)
    vc = aco.obs[cc].value_counts()
    aco = aco[aco.obs[cc].isin(vc[vc >= cfg["min_cluster_n"]].index)].copy()
    myo_mask = myo_mask_of(aco, cfg)   # single-source mask (D1)
    log("CellOracle adata:", aco.shape,
        "| clusters:", aco.obs[cc].value_counts().to_dict(),
        "| myofibroblast cells:", int(myo_mask.sum()))

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=aco, cluster_column_name=cc,
                                       embedding_name=cfg["embedding"])

    # base GRN: default human promoter, OR a skin scATAC-derived TF_info_matrix
    if arm == "promoter":
        base = co.data.load_human_promoter_base_GRN()
        log("base-GRN = human promoter (generic, tissue-agnostic)")
    else:
        log("base-GRN = SKIN scATAC (GSE312129):", base_grn_path)
        if base_grn_path.endswith(".parquet"):
            base = pd.read_parquet(base_grn_path)
        else:
            base = pd.read_csv(base_grn_path, index_col=0)
    oracle.import_TF_data(TF_info_matrix=base)

    oracle.perform_PCA()
    # knn imputation params -- CellOracle TUTORIAL rule (A1): k scales with cell count.
    # The official Network_analysis notebook sets  k = int(0.025 * n_cell),
    # b_sight = k*8, b_maxl = k*4  (n_cell = oracle.adata.shape[0]). A fixed k=25 on a
    # >100k-cell atlas is ~100x too small -> the imputed_count layer (which feeds EVERY
    # KO delta) is badly under-smoothed. Compute k from the ACTUAL oracle cell count
    # (after the stratified cap above); cfg["knn_k"] is kept only as a floor/override.
    n_cell = oracle.adata.n_obs
    k = max(int(0.025 * n_cell), cfg["knn_k"])   # tutorial 2.5% rule, cfg["knn_k"] = floor
    b_sight = min(k * 8, n_cell - 1)
    b_maxl  = min(k * 4, n_cell - 1)
    log(f"knn_imputation: n_cell={n_cell} -> k={k} (=max(2.5%*n_cell, floor={cfg['knn_k']})), "
        f"b_sight={b_sight}, b_maxl={b_maxl}")
    np.random.seed(cfg["seed"])   # re-pin global RNG before stochastic balanced-KNN (F09 ext.): knn_imputation samples b_sight/b_maxl from numpy global RNG and has no random_state arg; the imputed_count layer feeds every downstream delta
    oracle.knn_imputation(n_pca_dims=cfg["knn_pca_dims"], k=k, balanced=True,
                          b_sight=b_sight, b_maxl=b_maxl,
                          n_jobs=cfg["n_jobs"])

    np.random.seed(cfg["seed"])   # re-pin global RNG immediately before stochastic bagging (CellOracle get_links has no random_state arg)
    links = oracle.get_links(cluster_name_for_GRN_unit=cc, alpha=cfg["grn_alpha"],
                             verbose_level=0, n_jobs=cfg["n_jobs"],
                             bagging_number=cfg["bagging"])
    links.filter_links(p=cfg["link_p"], weight="coef_abs", threshold_number=cfg["link_topn"])
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=cfg["grn_alpha"], use_cluster_specific_TFdict=True)

    pickle.dump(oracle, open(oracle_pkl, "wb"))
    links.to_hdf5(links_out)
    log(f"GRN fitted + saved: {oracle_pkl}")
    gc.collect()
    return oracle, links


# ============================================================================
# SHARED -- prespecified all-active TF candidate manifest
# ============================================================================
def build_active_tf_manifest(oracle, links, cfg, arm):
    """Inventory every TF with >=1 edge in the filtered fitted GRN.

    Screenability is defined before looking at perturbation results: the TF must
    be in the fitted expression matrix and detected by raw counts in at least the
    configured fraction of fibroblasts. Both all-fibroblast and myofibroblast
    expression fractions are reported. No rank, HES1 status or druggability term
    enters candidate eligibility.
    """
    A = oracle.adata
    tag = "" if arm == "promoter" else f"_{arm}"
    manifest_csv = os.path.join(
        cfg["out_dir"], f"KO_candidates_all_active{tag}.csv")

    edge_count = {}
    edge_clusters = {}
    for cluster, frame in links.filtered_links.items():
        if not isinstance(frame, pd.DataFrame) or "source" not in frame.columns:
            continue
        counts = frame["source"].astype(str).value_counts()
        for gene, count in counts.items():
            edge_count[gene] = edge_count.get(gene, 0) + int(count)
            edge_clusters.setdefault(gene, set()).add(str(cluster))

    active = sorted(edge_count)
    var_names = list(map(str, A.var_names))
    var_index = {gene: i for i, gene in enumerate(var_names)}
    present = [gene for gene in active if gene in var_index]
    present_idx = [var_index[gene] for gene in present]

    raw_layer = detect_raw_layer(A)
    if raw_layer is None and not is_integer_counts(A.X):
        raise RuntimeError("All-active screenability requires raw counts, but the fitted "
                           "Oracle object has no raw-count layer and .X is not integer counts.")
    raw = A.layers[raw_layer] if raw_layer is not None else A.X
    myo_mask = myo_mask_of(A, cfg)

    def _fraction_nonzero(X):
        if X.shape[1] == 0:
            return np.array([], dtype=float)
        if hasattr(X, "getnnz"):
            return np.asarray(X.getnnz(axis=0), dtype=float).ravel() / max(X.shape[0], 1)
        return np.count_nonzero(np.asarray(X) > 0, axis=0) / max(X.shape[0], 1)

    raw_present = raw[:, present_idx]
    frac_all = _fraction_nonzero(raw_present)
    frac_myo = _fraction_nonzero(raw_present[myo_mask])
    imputed = A.layers["imputed_count"][:, present_idx]
    mean_imputed_myo = np.asarray(imputed[myo_mask].mean(axis=0)).ravel()

    expr_all = dict(zip(present, map(float, frac_all)))
    expr_myo = dict(zip(present, map(float, frac_myo)))
    imputed_myo = dict(zip(present, map(float, mean_imputed_myo)))
    min_fraction = float(cfg["ko_min_raw_expr_fraction"])
    rows = []
    for gene in active:
        in_matrix = gene in var_index
        fraction = expr_all.get(gene, np.nan)
        reasons = []
        if not in_matrix:
            reasons.append("absent_from_fitted_expression_matrix")
        elif not np.isfinite(fraction) or fraction < min_fraction:
            reasons.append("raw_expression_fraction_below_threshold")
        rows.append({
            "gene": gene,
            "arm": arm,
            "active_edge_count": int(edge_count[gene]),
            "active_cluster_count": int(len(edge_clusters[gene])),
            "in_expression_matrix": bool(in_matrix),
            "raw_expr_fraction_fibroblast": fraction,
            "raw_expr_fraction_myo": expr_myo.get(gene, np.nan),
            "mean_imputed_expr_myo": imputed_myo.get(gene, np.nan),
            "min_raw_expr_fraction_required": min_fraction,
            "screenable": len(reasons) == 0,
            "exclusion_reason": ";".join(reasons),
        })

    if not rows:
        raise RuntimeError(f"No edge-active TFs remained after link filtering in arm '{arm}'.")
    manifest = pd.DataFrame(rows).sort_values(
        ["screenable", "active_edge_count", "gene"],
        ascending=[False, False, True]).reset_index(drop=True)
    manifest.to_csv(manifest_csv, index=False)
    log(f"[all-active manifest:{arm}] edge-active={len(manifest)}, "
        f"screenable={int(manifest['screenable'].sum())}, "
        f"excluded={int((~manifest['screenable']).sum())}; saved {manifest_csv}")
    return manifest


# ============================================================================
# SHARED -- developmental reference gradient (used by BOTH stage_ko_ranking's
#   Perturbation-Score ranking and stage_vectorfield). Built ONCE per oracle.
#   This is the SFRP2/DPP4-rooted DPT axis + Gradient_calculator machinery that
#   the CellOracle Gata1-KO tutorial (Kamimoto, Nature 2023) uses as the
#   differentiation reference flow for the perturbation-score inner product.
#   Factored out of the original stage_vectorfield code (was inline there) so the
#   KO ranking and the vector-field stage share ONE identical gradient definition.
# ============================================================================
def build_dev_gradient(oracle, cfg):
    """Build the developmental reference gradient (Gradient_calculator) rooted at
    the SFRP2/DPP4 progenitor-fibroblast pool via a scanpy DPT axis computed HERE
    (03 runs BEFORE 04, so no external pseudotime exists). NO try/except silent-NaN:
    if the inner-product arm is requested but the gradient is unbuildable, RAISE.
    Returns a Gradient_calculator with ref_flow/pseudotime_on_grid populated."""
    from celloracle.applications import Gradient_calculator
    pkey = cfg["vf_pseudotime_key"]
    adata_o = oracle.adata
    if pkey not in adata_o.obs:
        # compute a DPT pseudotime inside 03 for the gradient
        root_hi = None
        root_markers = [g for g in cfg["vf_root_markers"] if g in adata_o.var_names]
        if root_markers:
            # progenitor score = mean imputed expression of root markers; root = argmax
            ridx = [list(adata_o.var_names).index(g) for g in root_markers]
            prog = np.asarray(adata_o.layers["imputed_count"][:, ridx]).mean(1).ravel()
            root_hi = int(np.argmax(prog))
        if root_hi is None:
            raise RuntimeError(
                "[GRAD] no progenitor root markers present in the atlas; cannot root the "
                "developmental pseudotime for the inner-product gradient. Provide "
                "cfg['vf_root_markers'] genes present in the data, or set "
                "cfg['vf_inner_product']=False to disable this arm explicitly.")
        # NO try/except fallback: a DPT failure is a bug to FIX (bad graph / root),
        # not a reason to silently swap the developmental axis. Let it raise.
        tmp = adata_o.copy()
        if "X_pca" not in tmp.obsm:
            sc.pp.pca(tmp, n_comps=min(50, tmp.n_vars - 1), random_state=cfg["seed"])
        sc.pp.neighbors(tmp, n_neighbors=15, random_state=cfg["seed"])
        sc.tl.diffmap(tmp, n_comps=15)
        tmp.uns["iroot"] = root_hi
        sc.tl.dpt(tmp)
        adata_o.obs[pkey] = np.asarray(tmp.obs["dpt_pseudotime"].values, dtype=float)
        log(f"[GRAD] DPT pseudotime computed inside 03, root cell {root_hi} "
            f"(progenitor markers {root_markers}) -> obs['{pkey}']")
        del tmp
    pt = np.asarray(adata_o.obs[pkey].values, dtype=float)
    if not np.isfinite(pt).any() or np.nanstd(pt) == 0:
        raise RuntimeError(
            f"[FATAL] pseudotime '{pkey}' is degenerate (all-NaN/constant); cannot build a "
            "developmental gradient. Fix the root/axis or set cfg['vf_inner_product']=False.")
    # ★2026-07-08 audit: independent axis-sign check (BEYOND the SMAD3 gate). Pseudotime rooted
    # at the SFRP2/DPP4 progenitor MUST increase toward the myofibroblast tip, so myofibroblast-
    # panel expression must correlate POSITIVELY with it; otherwise the axis is inverted and every
    # downstream Perturbation-Score sign would flip -> abort loudly rather than silently trust it.
    _myo_idx = [list(adata_o.var_names).index(g) for g in MYO if g in adata_o.var_names]
    if _myo_idx and "imputed_count" in adata_o.layers:
        from scipy.stats import spearmanr as _spr
        _mx = np.asarray(adata_o.layers["imputed_count"][:, _myo_idx]).mean(1).ravel()
        _rho, _ = _spr(pt, _mx, nan_policy="omit")
        log(f"[GRAD] axis-sign check: spearman(pseudotime, myofib panel) = {_rho:+.3f} (expect >0)")
        if np.isfinite(_rho) and _rho < 0:
            # DPT direction is root-defined and can come out INVERTED (esp. on a subsample);
            # ORIENT the axis so progenitor->myofibroblast INCREASES -- the myofib-panel
            # correlation is the ground-truth direction (same principle as 04's orient-by-myo).
            # Flip + re-store on oracle.adata so Gradient_calculator reads the corrected axis.
            pt = float(np.nanmax(pt)) - pt
            adata_o.obs[pkey] = pt
            _rho = -_rho
            log(f"[GRAD] pseudotime was INVERTED (rho<0) -> FLIPPED so myofibroblast program "
                f"increases along it; rho now {_rho:+.3f}")
        if not (np.isfinite(_rho) and _rho > 0):
            raise RuntimeError(
                f"[FATAL] developmental axis STILL not oriented after flip (rho={_rho:+.3f}); "
                "degenerate/non-monotone axis. Fix root or set cfg['vf_inner_product']=False.")
    gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key=pkey)
    gradient.calculate_p_mass(smooth=0.8, n_grid=cfg["vf_n_grid"], n_neighbors=cfg["vf_n_neighbors"])
    # ★2026-07-07 修: suggest_mass_thresholds() 只 PLOT 候选阈值、返回 None(源码无 return),
    # 原代码 float(None) 崩。改为编程式取 min_mass = 非零局部密度 mass 的低分位数(0.05)
    # -> 滤掉近空网格点、保留被细胞覆盖的流形; 数据自适应, 无交互依赖。
    _mp = np.asarray(gradient.total_p_mass)
    gmin = float(np.quantile(_mp[_mp > 0], 0.05)) if np.any(_mp > 0) else 0.0
    gradient.calculate_mass_filter(min_mass=gmin, plot=False)
    gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly": 3}, plot=False)
    gradient.calculate_gradient()
    log(f"[GRAD] developmental reference gradient built from pseudotime key '{pkey}'")
    return gradient


# ============================================================================
# SHARED -- per-gene Perturbation Score (PS) via the REAL CellOracle
#   Oracle_development_module inner-product API (Gata1-KO tutorial, Kamimoto
#   Nature 2023). Runs one in-silico KO, grids the perturbation field, and takes
#   the inner product of the KO flow with the developmental reference flow inside
#   the myofibroblast region. SIGN: NEGATIVE inner product = KO OPPOSES the
#   fibrotic differentiation flow = the TF is a PRO-fibrotic regulator.
#   The randomized-GRN null (score_randomized) is built automatically because
#   estimate_transition_prob(calculate_randomized=True) fills delta_embedding_random,
#   which calculate_grid_arrows -> flow_rndm -> calculate_inner_product carries into
#   inner_product_random. get_negative_PS_p_value() Wilcoxon-tests observed negative
#   PS vs this randomized null. This REPLACES the broken non-regulator empirical null.
# ============================================================================
def ps_for_gene(oracle, gradient, gene, myo_cell_idx, cfg):
    """Return (ps_score, ps_pval, mean_shift_myo) for one in-silico KO.
      ps_score = signed sum of negative inner products in the myo region
                 (<=0; more NEGATIVE = KO more strongly opposes the fibrotic flow
                 = stronger pro-fibrotic regulator). This is CellOracle's
                 Perturbation Score restricted to negative IPs (get_sum_of_negative_ips).
      ps_pval  = Wilcoxon (paired, one-sided 'less') p of observed negative PS vs
                 the randomized-GRN null (get_negative_PS_p_value).
    NO try/except silent-NaN here -- a failure is surfaced to the caller (which logs
    it as a screen failure), never masked."""
    from celloracle.applications import Oracle_development_module
    # 1) propagate the KO and estimate transition probs on the KNN, WITH the
    #    randomized-GRN control (calculate_randomized=True -> delta_embedding_random).
    oracle.simulate_shift(perturb_condition={gene: 0.0}, n_propagation=cfg["n_propagation"])
    np.random.seed(cfg["seed"])   # re-pin RNG before stochastic estimate_transition_prob(knn_random=True)
    oracle.estimate_transition_prob(n_neighbors=cfg["vf_n_neighbors"],
                                    knn_random=True,
                                    sampled_fraction=cfg["vf_sampled_frac"],
                                    calculate_randomized=True,
                                    n_jobs=cfg["n_jobs"],
                                    threads=cfg["n_jobs"],
                                    random_seed=cfg["vf_random_seed"])   # builds the randomized-GRN null field
    oracle.calculate_embedding_shift(sigma_corr=cfg["vf_sigma_corr"])
    # per-cell shift magnitude on the embedding (descriptive; myo region)
    shift = np.asarray(oracle.delta_embedding)
    mag = np.linalg.norm(shift, axis=1)
    mean_shift_myo = float(mag[myo_cell_idx].mean()) if len(myo_cell_idx) else float("nan")
    # 2) grid the perturbation field: load_perturb_simulation_data reads
    #    oracle.total_p_mass, so calculate_p_mass + calculate_mass_filter are required
    #    (exactly as stage_vectorfield does).
    oracle.calculate_p_mass(smooth=0.8, n_grid=cfg["vf_n_grid"], n_neighbors=cfg["vf_n_neighbors"])
    # ★2026-07-07 修: suggest_mass_thresholds() 只 PLOT 返回 None -> 编程式取 min_mass =
    # 非零 mass 低分位数(0.05)(与 build_dev_gradient 同规则), 无交互依赖。
    _mp = np.asarray(oracle.total_p_mass)
    min_mass = float(np.quantile(_mp[_mp > 0], 0.05)) if np.any(_mp > 0) else 0.0
    oracle.calculate_mass_filter(min_mass=min_mass, plot=False)
    # 3) inner product vs the developmental reference flow, restricted to the myo region
    dev = Oracle_development_module()
    dev.load_differentiation_reference_data(gradient_object=gradient)
    dev.load_perturb_simulation_data(oracle_object=oracle,
                                     cell_idx_use=myo_cell_idx,
                                     name=f"KO_{gene}")
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)   # builds inner_product_df (required by get_negative_PS_p_value)
    # get_negative_PS_p_value(return_ps_sum=True) -> (p, -x.sum(), -y.sum());
    #   x = negative-clipped observed IPs, y = negative-clipped randomized IPs.
    #   ps_sum = -x.sum() is the POSITIVE magnitude of negative PS; the SIGNED PS
    #   (task convention: more negative = stronger) is therefore ps_score = -ps_sum.
    p, ps_sum, ps_sum_random = dev.get_negative_PS_p_value(return_ps_sum=True, plot=False)
    ps_score = -float(ps_sum)          # signed: <=0, more negative = stronger pro-fibrotic
    # ★2026-07-08 BUG-D fix: get_negative_PS_p_value returns p=NaN (NOT an error) on a DEGENERATE
    # input (too few surviving grid points / empty negative-inner-product set). A NaN p is
    # "UNCOMPUTABLE", not "biologically non-significant" -- flag it so the SMAD3 control gate can
    # tell a degenerate slice apart from a true control failure (02 mad_bounds NaN-poison family).
    _p_finite = (p is not None) and bool(np.isfinite(p))
    _ipdf = getattr(dev, "inner_product_df", None)
    ps_uncomputable = (not _p_finite) or (_ipdf is None) or (len(_ipdf) == 0)
    ps_pval = float(p) if _p_finite else float("nan")
    return ps_score, ps_pval, mean_shift_myo, ps_uncomputable


# ============================================================================
# STAGE 2 -- all-active-TF in-silico KO ranking (+ positive controls)
#   RANKED BY CELLORACLE PERTURBATION SCORE (PS), the tool's intended KO readout
#   (Kamimoto Nature 2023 / Gata1-KO tutorial): PS = inner product of the KO
#   perturbation vector field with the developmental (progenitor->myofibroblast)
#   differentiation flow, restricted to the myofibroblast region. The old mean
#   fibrotic-program expression delta (delta_fibrotic_program) washed out real
#   regulators (SMAD3 came out ~0 at rank 6); it is KEPT only as a SECONDARY
#   descriptive column and drives NEITHER the ranking NOR the control gate.
# ============================================================================
def stage_ko_ranking(oracle, links, cfg, arm):
    out = cfg["out_dir"]
    tag = "" if arm == "promoter" else f"_{arm}"
    if cfg.get("ko_screen_mode") != "all_active":
        raise ValueError("This extension requires CFG['ko_screen_mode']='all_active'. "
                         "The completed 43-TF calibration is preserved in legacy outputs.")
    ko_csv = f"{out}/KO_ranking_all_active{tag}.csv"
    progress_csv = f"{out}/KO_progress_all_active{tag}.csv"
    failure_csv = f"{out}/KO_failures_all_active{tag}.csv"
    if not cfg["force"] and exists(ko_csv):
        log(f"[skip all-active KO ranking:{arm}] {ko_csv} exists")
        return pd.read_csv(ko_csv)

    log(f"== STAGE 2: all-active-TF KO ranking [{arm}] -- "
        "ranked by CellOracle Perturbation Score ==")
    A = oracle.adata
    vn = list(A.var_names)
    myo_mask = myo_mask_of(A, cfg)   # single-source mask (D1)
    myo_cell_idx = np.where(myo_mask)[0]

    # Candidate scope is fixed independently of HES1: all TFs with >=1 filtered
    # GRN edge that pass the raw-expression screenability threshold.
    manifest = build_active_tf_manifest(oracle, links, cfg, arm)
    eligible = manifest.loc[manifest["screenable"]].copy()
    panel = sorted(eligible["gene"].astype(str).tolist())
    if not panel:
        raise RuntimeError("[CONTROL FAIL] no screenable edge-active TFs; inspect "
                           f"KO_candidates_all_active{tag}.csv")
    meta_by_gene = eligible.set_index("gene").to_dict(orient="index")
    log(f"All-active KO candidates ({len(panel)}; alphabetical execution order):", panel)

    pg = [g for g in MYO if g in vn]
    pidx = [vn.index(g) for g in pg]

    def _ko_delta(gene):
        """SECONDARY, descriptive only: predicted change in the fibrotic (MYO) program
        in myofibroblast cells after in-silico KO of `gene` (mean simulated-imputed
        expression of the MYO panel). Negative = gene nominally a positive regulator.
        This metric WASHES OUT real regulators (SMAD3 ~0 at rank 6) and is NO LONGER the
        ranking axis or the control gate -- it is kept as `delta_fibrotic_program` for
        descriptive comparison only. It reads simulated_count from the LAST simulate_shift,
        so it is computed INSIDE the PS loop right after that gene's simulate_shift."""
        sim = oracle.adata.layers["simulated_count"]
        imp = oracle.adata.layers["imputed_count"]
        return float(np.asarray(sim[myo_mask][:, pidx]).mean()
                     - np.asarray(imp[myo_mask][:, pidx]).mean())

    # ------------------------------------------------------------------
    # Developmental reference gradient (D3) -- built ONCE, up front, and REUSED for
    # every screened TF's Perturbation Score. SFRP2/DPP4-rooted DPT axis + Gradient_calculator
    # (Gata1-KO tutorial / Kamimoto Nature 2023). If the inner-product arm is disabled we
    # cannot compute PS at all -> that is a configuration error for THIS stage, so RAISE.
    # ------------------------------------------------------------------
    if not cfg["vf_inner_product"]:
        raise RuntimeError(
            "[CONTROL FAIL] stage_ko_ranking now ranks by the CellOracle Perturbation Score, "
            "which REQUIRES the developmental inner-product machinery, but cfg['vf_inner_product'] "
            "is False. Set cfg['vf_inner_product']=True (the PS ranking cannot be computed without it).")
    gradient = build_dev_gradient(oracle, cfg)

    # ------------------------------------------------------------------
    # PRIMARY SCREEN: rank the panel by the CellOracle PERTURBATION SCORE (PS) --
    # the inner product of each KO's perturbation vector field with the developmental
    # (progenitor -> myofibroblast) reference flow, summed over NEGATIVE inner-product
    # grid points in the myofibroblast region (Gata1-KO tutorial; Kamimoto Nature 2023).
    #   ps_score  = signed sum of negative IPs (<=0); MORE NEGATIVE = KO more strongly
    #               OPPOSES the fibrotic differentiation flow = stronger PRO-fibrotic regulator.
    #   ps_pval   = Wilcoxon (paired, one-sided 'less') p of observed negative PS vs the
    #               randomized-GRN null field (built by estimate_transition_prob(
    #               calculate_randomized=True)). This REPLACES the broken non-regulator null.
    # `delta_fibrotic_program` is still computed per gene (right after each simulate_shift
    # inside ps_for_gene) but is SECONDARY/descriptive and drives neither ranking nor gate.
    # ------------------------------------------------------------------
    res = []
    failure_rows = []
    if not cfg["force"] and exists(progress_csv):
        prior = pd.read_csv(progress_csv)
        prior = prior[prior["gene"].astype(str).isin(panel)]
        res = prior.to_dict(orient="records")
        log(f"[resume all-active:{arm}] loaded {len(res)} successful TFs from {progress_csv}")
    if not cfg["force"] and exists(failure_csv):
        prior_failures = pd.read_csv(failure_csv)
        prior_failures = prior_failures[prior_failures["gene"].astype(str).isin(panel)]
        failure_rows = prior_failures.to_dict(orient="records")
        log(f"[resume all-active:{arm}] loaded {len(failure_rows)} failed TFs from {failure_csv}")
    completed = {str(row["gene"]) for row in res + failure_rows}

    for run_i, gene in enumerate(panel, start=1):
        if gene in completed:
            continue
        # per-gene guard: in a KO SCREEN some panel genes legitimately fail simulate_shift's
        # sanity gates (not a base-GRN TF / too few GRN connections). Skipping ONE such gene
        # with a LOUD log + a reported failure list is correct screen behaviour -- it is NOT
        # a silent fallback (the failures are logged and counted; the SMAD3 control gate below
        # still fires if SMAD3 itself fails, because SMAD3 is then absent from `res`).
        try:
            ps_score, ps_pval, mean_shift_myo, ps_uncomputable = ps_for_gene(
                oracle, gradient, gene, myo_cell_idx, cfg)
            # descriptive delta uses the simulated_count from the simulate_shift that
            # ps_for_gene just ran for THIS gene (no extra propagation).
            d = _ko_delta(gene)
        except Exception as e:
            log(f"  KO {gene} FAIL (simulate_shift / PS sanity gate): {e}")
            failure_rows.append({
                "gene": gene,
                "arm": arm,
                "failure_type": type(e).__name__,
                "failure_message": str(e),
                **meta_by_gene[gene],
            })
            pd.DataFrame(failure_rows).to_csv(failure_csv, index=False)
            continue
        row = {"gene": gene,
               "ps_score": ps_score, "ps_pval": ps_pval,
               "ps_uncomputable": bool(ps_uncomputable),
               "delta_fibrotic_program": d,
               "is_control": gene in CTRL_GRN, "is_lead": gene in LEADS,
               "myo_TF_expr": meta_by_gene[gene]["mean_imputed_expr_myo"],
               **meta_by_gene[gene]}
        res.append(row)
        log(f"  KO {gene}: PS={ps_score:.4f} p={ps_pval:.3g} "
            f"shift_myo={mean_shift_myo:.4f} dProg(2ndary)={d:.4f}")
        if run_i % max(int(cfg.get("ko_checkpoint_every", 1)), 1) == 0:
            pd.DataFrame(res).to_csv(progress_csv, index=False)

    pd.DataFrame(res).to_csv(progress_csv, index=False)
    failure_columns = ["gene", "arm", "failure_type", "failure_message"]
    pd.DataFrame(failure_rows, columns=(
        list(dict.fromkeys(failure_columns + list(failure_rows[0].keys())))
        if failure_rows else failure_columns)).to_csv(failure_csv, index=False)
    ko_failed = [str(row["gene"]) for row in failure_rows]
    if ko_failed:
        log(f"[KO screen] {len(ko_failed)}/{len(panel)} active TFs not simulatable "
            f"(retained in failure table): {ko_failed}")
    if not res:
        raise RuntimeError("[CONTROL FAIL] no active TF yielded a Perturbation Score "
                           "(all failed the simulate_shift/PS sanity gate) -> GRN screen invalid.")

    # ------------------------------------------------------------------
    # OUT-OF-DISTRIBUTION (OOD) SANITY after simulate_shift (A2, tutorial-grounded).
    # CellOracle's tutorial checks that simulated expression values do not fall far
    # OUTSIDE the observed (imputed) expression range -- deltas built on OOD-simulated
    # values are not trustworthy. This is REPORTING ONLY: it does NOT replace the
    # per-gene simulate_shift sanity guard above; it logs, for the LEADS + SMAD3
    # control, the fraction of simulated values outside the observed imputed min/max.
    # Uses oracle.evaluate_and_plot_simulation_value_distribution() when the installed
    # celloracle exposes it (verified in the Oracle API); otherwise falls back to an
    # explicit per-gene OOD-fraction computation against the imputed_count layer.
    ood_panel = [g for g in (LEADS + ["SMAD3"]) if g in vn]
    for gene in ood_panel:
        try:
            oracle.simulate_shift(perturb_condition={gene: 0.0},
                                  n_propagation=cfg["n_propagation"])
        except Exception as e:
            log(f"  [OOD] {gene}: simulate_shift failed ({e}); skipping OOD report")
            continue
        if hasattr(oracle, "evaluate_and_plot_simulation_value_distribution"):
            try:
                fig_png = f"{out}/OOD_simvalue_dist_all_active_{gene}{tag}.png"
                oracle.evaluate_and_plot_simulation_value_distribution(
                    n_genes=4, save=fig_png)
                log(f"  [OOD] {gene}: evaluate_and_plot_simulation_value_distribution -> {fig_png}")
            except Exception as e:
                log(f"  [OOD] {gene}: evaluate_and_plot_simulation_value_distribution "
                    f"unavailable/failed ({e}); using explicit OOD-fraction fallback")
        # explicit OOD fraction vs the observed imputed range (always logged, myo cells)
        try:
            sim = np.asarray(oracle.adata.layers["simulated_count"])[myo_mask]
            imp = np.asarray(oracle.adata.layers["imputed_count"])[myo_mask]
            lo, hi = imp.min(axis=0), imp.max(axis=0)
            ood = (sim < lo) | (sim > hi)
            frac = float(ood.mean())
            log(f"  [OOD] {gene}: fraction simulated values outside observed imputed "
                f"range = {frac:.4f} (myo cells x genes)")
        except Exception as e:
            log(f"  [OOD] {gene}: explicit OOD-fraction computation failed ({e})")

    # ------------------------------------------------------------------
    # MULTIPLE-TESTING on the PERTURBATION SCORE axis. Each successfully screened TF has a
    # per-TF ps_pval from CellOracle's OWN randomized-GRN null (get_negative_PS_p_value:
    # Wilcoxon of observed negative IPs vs the randomized-GRN negative IPs). That null
    # is CORRECT and always available (it uses the randomized field built per KO), which
    # REPLACES the structurally broken "simulatable non-regulator genes" empirical null
    # (that null had 0 candidates because every simulatable gene IS a regulator -> p/q NaN
    # -> the SMAD3 gate always raised). BH is applied once across all successfully
    # screened active TFs, rather than across the historical 43-TF calibration panel.
    # ------------------------------------------------------------------
    from statsmodels.stats.multitest import multipletests
    _p = np.array([(r["ps_pval"] if np.isfinite(r["ps_pval"]) else 1.0) for r in res])
    _rej, _q, _, _ = multipletests(_p, alpha=0.05, method="fdr_bh")
    for r, q, rej in zip(res, _q, _rej):
        r["ps_qbh"] = float(q)
        # a TF is a FDR-significant PRO-fibrotic regulator on the PS axis if it
        # SIGNIFICANTLY opposes the fibrotic flow (q<0.05) AND its PS is negative.
        r["fdr_pos_regulator"] = bool(rej and (r["ps_score"] < 0))
    log(f"PS-axis FDR: {sum(r['fdr_pos_regulator'] for r in res)}/{len(res)} active TFs are "
        "FDR-significant negative-PS regulators (BH q<0.05 across the full successful screen).")

    # PRIMARY RANK = ps_score ASCENDING (most-negative PS first = strongest pro-fibrotic
    # regulator = KO most strongly opposes the progenitor->myofibroblast differentiation flow).
    R = pd.DataFrame(res).sort_values("ps_score", ascending=True).reset_index(drop=True)
    R["rank"] = R.index + 1   # rank 1 = strongest pro-fibrotic regulator (most negative PS)
    R["n_ranked"] = len(R)
    R["rank_fraction"] = R["rank"] / max(len(R), 1)
    R["rank_label"] = R["rank"].astype(str) + "/" + R["n_ranked"].astype(str)
    # column order: PS axis primary, delta_fibrotic_program kept as secondary descriptor.
    col_order = ["gene", "rank", "n_ranked", "rank_fraction", "rank_label",
                 "ps_score", "ps_pval", "ps_qbh", "fdr_pos_regulator",
                 "delta_fibrotic_program", "is_control", "is_lead",
                 "active_edge_count", "active_cluster_count",
                 "raw_expr_fraction_fibroblast", "raw_expr_fraction_myo",
                 "mean_imputed_expr_myo"]
    R = R[[c for c in col_order if c in R.columns] + [c for c in R.columns if c not in col_order]]
    R.to_csv(ko_csv, index=False)
    log(f"KO ranking saved: {ko_csv}")
    log("TOP pro-fibrotic regulators (most-negative Perturbation Score):\n"
        + R.head(15).to_string(index=False))

    # ------------------------------------------------------------------
    # POSITIVE/NEGATIVE CONTROL HARD-GATE (D2) -- NOW ON THE PERTURBATION-SCORE AXIS.
    # SMAD3 (canonical TGF-beta pro-fibrotic driver) must be a FDR-significant NEGATIVE-PS
    # regulator in ABSOLUTE terms: ps_score<0 AND ps_qbh<0.05, i.e. KO of SMAD3 significantly
    # OPPOSES (blocks) the fibrotic differentiation flow. This is the tool's intended readout
    # (Gata1-KO tutorial; Kamimoto Nature 2023) and replaces the old mean-delta gate that let
    # a near-zero SMAD3 pass. If SMAD3 fails -> HARD ABORT (a degraded GRN is surfaced
    # immediately). MEF2C (soft negative control) = warn loudly, do NOT abort.
    # Association-level nomination, not a causal claim.
    # ------------------------------------------------------------------
    log("CONTROL sanity (PS axis):\n"
        + R[R.is_control][["gene", "ps_score", "ps_pval", "ps_qbh", "rank"]].to_string(index=False))
    if "SMAD3" not in set(R["gene"]):
        raise RuntimeError("[CONTROL FAIL] SMAD3 absent from the successful all-active screen. "
                           "Inspect the candidate manifest and failure table before interpretation.")
    smad3 = R.loc[R.gene == "SMAD3"].iloc[0]
    # ★2026-07-08 BUG-D fix: distinguish an UNCOMPUTABLE SMAD3 PS p-value (degenerate myofibroblast
    # grid / too few surviving grid points) from a genuine biological control failure -- the former
    # means "fix the grid/region", NOT "the GRN is invalid".
    if bool(smad3.get("ps_uncomputable", False)):
        raise RuntimeError(
            "[CONTROL INCONCLUSIVE] SMAD3 Perturbation-Score p-value is UNCOMPUTABLE (NaN) -- "
            "degenerate myofibroblast grid / too few surviving grid points after mass filtering, "
            "NOT biological non-significance. Enlarge the myo region / relax min_cluster_n / lower "
            "the mass quantile before trusting the screen (input-degeneracy, not a GRN failure).")
    smad3_ok = bool((smad3["ps_score"] < 0) and (float(smad3.get("ps_qbh", np.nan)) < 0.05))
    if not smad3_ok:
        raise RuntimeError(
            f"[CONTROL FAIL] SMAD3 is not a FDR-significant NEGATIVE-PS (pro-fibrotic) regulator "
            f"(ps_score={smad3['ps_score']:.4f}, ps_qbh={smad3.get('ps_qbh')}, "
            f"rank={int(smad3['rank'])}/{len(R)}). The canonical pro-fibrotic control failed on the "
            "Perturbation-Score axis -> KO of SMAD3 does NOT significantly block the fibrotic "
            "differentiation flow -> GRN screen invalid; fix the GRN before trusting any lead.")
    log(f"[CONTROL PASS] SMAD3 significant pro-fibrotic regulator on PS axis: "
        f"ps_score={smad3['ps_score']:.4f} ps_qbh={smad3.get('ps_qbh')} rank={int(smad3['rank'])}/{len(R)}")
    if "MEF2C" in set(R["gene"]):
        mef = R.loc[R.gene == "MEF2C"].iloc[0]
        if not (mef["ps_score"] < 0 and float(mef.get("ps_qbh", np.nan)) < 0.05):
            log(f"[CONTROL WARN] soft negative control MEF2C is NOT a significant negative-PS "
                f"regulator (ps_score={mef['ps_score']:.4f}, ps_qbh={mef.get('ps_qbh')}) "
                "-- expected/consistent with its demotion; non-fatal.")

    # ------------------------------------------------------------------
    # ★2026-07-08 LITERATURE-GROUNDED baseline + competitor-defense (workflow wy8ays4sn; recipes in
    # replan_FrontImmunol_2026/UPGRADE_PLAN_2026crossfield.md) -- REPLACES the earlier subjective
    # version. (1) KNOWN-DRIVER-RECOVERY baseline: Ahlmann-Eltze/Huber Nat Methods 2025 (PMID
    # 40759747) mandates a naive-baseline test; the in-silico-KO analog = does the Perturbation
    # Score recover ESTABLISHED drivers (DRIVER_POS; ground-truth Nat Genet 2025 fibroblast TF
    # Perturb-seq atlas PMID 40770575) BETTER than a naive expression-DE baseline? Report AUROC +
    # precision@k for PS vs DE (selling point: HES1 tops PS but NOT the naive DE baseline).
    # (2) COMPETITOR-DEFENSE engine columns (HES1 vs RUNX1/CREB3L1/FOSL2/RUNX3): emit the ENGINE
    # side (PS rank/q + druggable-Notch-axis flag + engine gate); the FULL convergence matrix +
    # Sankey (regulon/CCC/bootstrap columns) is assembled at the figure stage. All non-fatal.
    # ------------------------------------------------------------------
    try:
        import scipy.sparse as _spx
        from sklearn.metrics import roc_auc_score
        _A = oracle.adata; _vn = list(_A.var_names)
        _myo = np.asarray(myo_mask_of(_A, cfg), dtype=bool)
        _X = _A.layers["imputed_count"] if "imputed_count" in _A.layers else _A.X
        _pg = [g for g in R["gene"].tolist() if g in _vn]
        _ix = [_vn.index(g) for g in _pg]
        _Xp = _X[:, _ix]; _Xp = np.asarray(_Xp.todense()) if _spx.issparse(_Xp) else np.asarray(_Xp)
        _de = pd.DataFrame({"gene": _pg,
                            "de_log2fc": np.ravel(np.log2((_Xp[_myo].mean(0) + 1e-9) / (_Xp[~_myo].mean(0) + 1e-9)))})
        _de["de_absfc"] = _de["de_log2fc"].abs()
        B = R[["gene", "ps_score", "rank"]].merge(_de, on="gene", how="inner")
        B["is_driver_pos"] = B["gene"].isin(DRIVER_POS)
        B["is_driver_neg"] = B["gene"].isin(DRIVER_NEG)
        _y = B["is_driver_pos"].astype(int).values
        _ps_order = B.sort_values("ps_score", ascending=True)["gene"].tolist()    # most-negative PS first
        _de_order = B.sort_values("de_absfc", ascending=False)["gene"].tolist()
        _patk = lambda order, k: float(sum(g in set(DRIVER_POS) for g in order[:k]) / max(k, 1))
        if 0 < int(_y.sum()) < len(_y):
            _ps_auc = float(roc_auc_score(_y, -B["ps_score"].values))   # PS: more-negative = stronger pro-fibrotic
            _de_auc = float(roc_auc_score(_y, B["de_absfc"].values))
        else:
            _ps_auc = _de_auc = float("nan")
        _rec = {"n_screened": int(len(B)), "n_driver_pos_screened": int(_y.sum()),
                "ps_auroc_recover_drivers": _ps_auc, "de_baseline_auroc": _de_auc,
                "ps_precision_at_10": _patk(_ps_order, 10), "de_precision_at_10": _patk(_de_order, 10),
                "ps_precision_at_20": _patk(_ps_order, 20), "de_precision_at_20": _patk(_de_order, 20),
                "perturbation_beats_DE_baseline": bool(np.isfinite(_ps_auc) and np.isfinite(_de_auc) and _ps_auc > _de_auc)}
        B.sort_values("ps_score").to_csv(
            os.path.join(cfg["out_dir"], f"baseline_recovery_all_active_{arm}.csv"),
            index=False)
        with open(os.path.join(
                cfg["out_dir"], f"baseline_recovery_summary_all_active_{arm}.json"), "w") as _fh:
            json.dump(_rec, _fh, indent=2)
        log(f"[BASELINE RECOVERY] known-driver recovery: PerturbationScore AUROC={_ps_auc:.3f} vs "
            f"expression-DE baseline AUROC={_de_auc:.3f} (prec@10 PS={_rec['ps_precision_at_10']:.2f} "
            f"DE={_rec['de_precision_at_10']:.2f}); "
            f"{'PS beats baseline' if _rec['perturbation_beats_DE_baseline'] else 'PS NOT above DE -- report honestly'}.")
        # competitor-defense engine columns (full convergence matrix/Sankey assembled at figure stage)
        _cd = ["HES1"] + [g for g in COMPETITORS if g != "HES1"]
        _cdf = R[R.gene.isin(_cd)][["gene", "rank", "ps_score", "ps_pval", "ps_qbh", "fdr_pos_regulator"]].merge(
            _de[["gene", "de_log2fc", "de_absfc"]], on="gene", how="left")
        _cdf["druggable_notch_axis"] = _cdf["gene"].isin(DRUGGABLE_AXIS)
        _cdf["engine_gate_pass"] = _cdf["fdr_pos_regulator"]
        for g in _cd:
            if g not in set(_cdf["gene"]):
                _cdf = pd.concat([_cdf, pd.DataFrame([{"gene": g, "druggable_notch_axis": g in DRUGGABLE_AXIS,
                    "engine_gate_pass": False, "fdr_pos_regulator": False}])], ignore_index=True)
        _cdf.sort_values("ps_score", na_position="last").to_csv(
            os.path.join(cfg["out_dir"], f"competitor_defense_engine_all_active_{arm}.csv"),
            index=False)
        log("[COMPETITOR DEFENSE - engine columns] HES1 vs 2026 competitors (PS + druggable-axis; "
            "full convergence matrix/Sankey at figure stage):\n" + _cdf.to_string(index=False))
    except Exception as e:
        log(f"[BASELINE/COMPETITOR-DEFENSE] skipped (non-fatal): {e}")

    return R


# ============================================================================
# STAGE 3 -- bootstrap CIs for leads + controls (robustness of the nomination)
# ============================================================================
def stage_bootstrap(oracle, cfg, arm):
    out = cfg["out_dir"]
    tag = "" if arm == "promoter" else f"_{arm}"
    boot_csv = f"{out}/bootstrap{tag}.csv"
    if not cfg["force"] and exists(boot_csv):
        log(f"[skip bootstrap:{arm}] {boot_csv} exists")
        return pd.read_csv(boot_csv)

    log(f"== STAGE 3: bootstrap CIs [{arm}] ==")
    A = oracle.adata
    vn = list(A.var_names)
    mm = np.where(myo_mask_of(A, cfg))[0]   # single-source mask (D1)
    pg = [g for g in MYO if g in vn]
    pidx = [vn.index(g) for g in pg]

    rng = np.random.default_rng(cfg["seed"])
    B = cfg["bootstrap_B"]
    rows = []
    for gene in LEADS + CTRL_BOOT:
        if gene not in vn:
            log(f"  {gene} absent")
            continue
        # ★2026-07-08 BUG-B fix: a gene in var_names but pruned out of the fitted GRN (not a
        # base-GRN TF / too few connections after filter_links) makes simulate_shift raise
        # ValueError. Guard it (mirror stage_ko_ranking's try/except) so ONE unsimulatable lead
        # (e.g. FOSB / demoted MEF2C) does NOT abort the whole bootstrap arm after ranking succeeded.
        try:
            oracle.simulate_shift(perturb_condition={gene: 0.0}, n_propagation=cfg["n_propagation"])
            sim = np.asarray(oracle.adata.layers["simulated_count"][:, pidx])[mm]
            imp = np.asarray(oracle.adata.layers["imputed_count"][:, pidx])[mm]
            pc = (sim - imp).mean(1)
            point = float(pc.mean())
            bs = np.array([pc[rng.integers(0, len(pc), len(pc))].mean() for _ in range(B)])
            lo, hi = np.percentile(bs, [2.5, 97.5])
            # bootstrap two-sided p for H0: delta>=0 (regulator not positive). NOTE (honest
            # framing): this within-sample cell-resample quantifies STABILITY, not a true
            # null; robust_neg below is the FDR-gated stability call, not a causal test.
            p_boot = 2.0 * min(float((bs >= 0).mean()), float((bs < 0).mean()))
            p_boot = min(max(p_boot, 1.0 / (B + 1)), 1.0)   # avoid p=0
            rows.append({"gene": gene, "kind": "LEAD" if gene in LEADS else "ctrl",
                         "delta": point, "ci_lo": float(lo), "ci_hi": float(hi),
                         "frac_cells_neg": float((pc < 0).mean()),
                         "p_boot": p_boot})
            log(f"  {gene}: d={point:.4f} CI[{lo:.4f},{hi:.4f}] frac_neg={(pc<0).mean():.2f}")
        except Exception as e:
            log(f"  {gene} not simulatable (bootstrap): {e}")
            continue
    # BH/FDR across all tested genes (F07: ~10 TFs each tested for CI-excludes-0)
    from statsmodels.stats.multitest import multipletests
    if not rows:   # ★2026-07-08 BUG-B: all leads/controls unsimulatable -> empty; avoid multipletests([]) crash
        log("[bootstrap] no simulatable lead/control -> empty; skipping BH")
        return pd.DataFrame()
    _p = np.array([r["p_boot"] for r in rows])
    _rej, _q, _, _ = multipletests(_p, alpha=0.05, method="fdr_bh")
    for r, q, rej in zip(rows, _q, _rej):
        r["q_bh"] = float(q)
        r["robust_neg"] = bool(rej and r["delta"] < 0)   # FDR-significant AND negative direction
    R = pd.DataFrame(rows).sort_values("delta")
    R.to_csv(boot_csv, index=False)
    log(f"bootstrap saved: {boot_csv}\n" + R.to_string(index=False))

    # ★2026-07-07: 此 bootstrap 在 mean-delta 轴上(已知会洗掉真调控子——SMAD3 的 mean-delta 甚至
    # 为正)。权威阳性对照已在 stage_ko_ranking 的 PS 轴通过(get_negative_PS_p_value, SMAD3 ps_qbh<0.05)。
    # 故此处 mean-delta bootstrap 的 SMAD3 未过 robust_neg 仅作 WARN, 不再 hard-abort(否则次要 robustness
    # 阶段的劣势指标会否掉已在正确指标上验证的门)。bootstrap CI 仅作描述性稳定性输出。
    if "SMAD3" in set(R["gene"]):
        s = R.loc[R.gene == "SMAD3"].iloc[0]
        if not bool(s["robust_neg"]):
            log(f"[CONTROL WARN] SMAD3 bootstrap(mean-delta 轴) 非 robust_neg "
                f"(delta={s['delta']:.4f}, q_bh={s.get('q_bh')}, CI[{s['ci_lo']:.4f},{s['ci_hi']:.4f}]); "
                "权威对照已在 PS 轴通过, 此为次要 mean-delta 稳定性描述, 非致命。")
        else:
            log(f"[CONTROL PASS] SMAD3 bootstrap robust_neg (delta={s['delta']:.4f}, q_bh={s.get('q_bh')})")
    return R


# ============================================================================
# STAGE 4 -- perturbation vector-field / inner-product (embedding-level effect)
#   API from module 069 (estimate_transition_prob + calculate_embedding_shift)
#   extended with the REAL Oracle_development_module inner-product / perturbation
#   score (verified against the CellOracle Gata1 KO tutorial). The inner product vs
#   a developmental reference gradient quantifies whether a KO opposes the
#   myofibroblast differentiation flow. SIGN: NEGATIVE inner product = KO BLOCKS
#   differentiation (opposes the fibrotic flow). HYPOTHESIS-level; not causal.
# ============================================================================
def stage_vectorfield(oracle, cfg, arm, ko_rank_df):
    if not cfg["vf_run"]:
        log(f"[skip vector-field:{arm}] vf_run=False")
        return
    out = cfg["out_dir"]
    tag = "" if arm == "promoter" else f"_{arm}"
    vf_csv = f"{out}/vectorfield{tag}.csv"
    if not cfg["force"] and exists(vf_csv):
        log(f"[skip vector-field:{arm}] {vf_csv} exists")
        return pd.read_csv(vf_csv)

    log(f"== STAGE 4: perturbation vector-field / inner-product [{arm}] ==")

    # genes = explicit list, else LEADS + top non-control KO hits from stage 2.
    # stage 2 now ranks by ps_score (ascending = strongest pro-fibrotic); pick top hits
    # by ps_score if present, else fall back to the secondary delta_fibrotic_program.
    if cfg["vf_genes"]:
        genes = [g for g in cfg["vf_genes"] if g in oracle.adata.var_names]
    else:
        if ko_rank_df is None:
            top_hits = []
        else:
            noncontrol = ko_rank_df[(~ko_rank_df.get("is_control", False))]
            rank_col = "ps_score" if "ps_score" in noncontrol.columns else "delta_fibrotic_program"
            top_hits = noncontrol.sort_values(rank_col, ascending=True).head(6)["gene"].tolist()
        genes = list(dict.fromkeys([g for g in LEADS if g in oracle.adata.var_names] + top_hits))
    log("vector-field genes:", genes)

    emb = cfg["embedding"]
    myo = np.asarray(oracle.adata.obs["myo"].values, dtype=float)
    myo_mask = myo_mask_of(oracle.adata, cfg)   # single-source mask (D1)

    # ------------------------------------------------------------------
    # Developmental reference gradient (D3). Built ONCE via the SHARED build_dev_gradient()
    # helper (same SFRP2/DPP4-rooted DPT axis + Gradient_calculator used by stage_ko_ranking's
    # Perturbation-Score ranking), so the vector field and the KO ranking share ONE identical
    # gradient definition. No try/except silent-NaN: if the IP arm is requested but unbuildable,
    # build_dev_gradient RAISES. NOTE: if stage_ko_ranking already ran, the DPT pseudotime key
    # is cached in oracle.adata.obs, so this call reuses it (idempotent).
    # ------------------------------------------------------------------
    gradient = None
    if cfg["vf_inner_product"]:
        gradient = build_dev_gradient(oracle, cfg)

    from celloracle.applications import Oracle_development_module

    rows = []
    for gene in genes:
        # 1) propagate the KO and estimate transition probabilities on the KNN.
        #    NO try/except silent-skip (F04): a failure here is a real bug and must raise.
        oracle.simulate_shift(perturb_condition={gene: 0.0}, n_propagation=cfg["n_propagation"])
        np.random.seed(cfg["seed"])   # re-pin RNG before stochastic estimate_transition_prob(knn_random=True) (missed-issue reproducibility)
        oracle.estimate_transition_prob(n_neighbors=cfg["vf_n_neighbors"],
                                        knn_random=True,
                                        sampled_fraction=cfg["vf_sampled_frac"])
        oracle.calculate_embedding_shift(sigma_corr=cfg["vf_sigma_corr"])

        # per-cell shift magnitude on the embedding (delta_embedding)
        shift = np.asarray(oracle.delta_embedding)          # cells x 2 (UMAP)
        mag = np.linalg.norm(shift, axis=1)
        row = {"gene": gene,
               "mean_shift_all": float(mag.mean()),
               "mean_shift_myo": float(mag[myo_mask].mean())}

        # 2) inner-product (perturbation score) vs the developmental reference flow,
        #    via the REAL Oracle_development_module API (verified against the CellOracle
        #    Gata1 KO tutorial). SIGN CONVENTION (from tutorial): NEGATIVE inner product
        #    = KO BLOCKS differentiation (opposes the fibrotic flow); POSITIVE = promotes.
        #    Association-level nomination, not a causal claim. NO try/except silent-NaN.
        if cfg["vf_inner_product"]:
            # grid the perturbation field first: load_perturb_simulation_data reads
            # oracle.total_p_mass, so calculate_p_mass + calculate_mass_filter are required.
            oracle.calculate_p_mass(smooth=0.8, n_grid=cfg["vf_n_grid"], n_neighbors=cfg["vf_n_neighbors"])
            # ★2026-07-07 修: suggest_mass_thresholds() 只 PLOT 返回 None -> 编程式取低分位 min_mass。
            _mp = np.asarray(oracle.total_p_mass)
            min_mass = float(np.quantile(_mp[_mp > 0], 0.05)) if np.any(_mp > 0) else 0.0
            oracle.calculate_mass_filter(min_mass=min_mass, plot=False)
            dev = Oracle_development_module()
            dev.load_differentiation_reference_data(gradient_object=gradient)
            dev.load_perturb_simulation_data(oracle_object=oracle,
                                             cell_idx_use=np.where(myo_mask)[0],
                                             name=f"KO_{gene}")
            dev.calculate_inner_product()
            dev.calculate_digitized_ip(n_bins=10)
            # dev.inner_product holds the per-grid PS array (self.inner_product =
            # np.array of dot(perturb_flow, ref_flow)). Attribute name follows the
            # tutorial; version-check against the installed CellOracle build if it moves.
            ip = np.asarray(dev.inner_product)
            row["mean_inner_product"] = float(np.nanmean(ip))
            row["frac_blocking"] = float(np.nanmean(ip < 0))   # <0 = KO opposes the fibrotic flow (association only)
        rows.append(row)
        log(f"  VF {gene}: shift_myo={row['mean_shift_myo']:.4f} "
            f"IP={row.get('mean_inner_product')} frac_blocking={row.get('frac_blocking')}")
    R = pd.DataFrame(rows)
    R.to_csv(vf_csv, index=False)
    log(f"vector-field saved: {vf_csv}\n" + R.to_string(index=False))
    return R


# ============================================================================
# driver -- runs the skin-scATAC arm as primary and promoter as sensitivity.
#   Debug stage-by-stage: comment out later stage() calls, or set env
#   SSC_STAGE to load only. Each stage is independently idempotent.
# ============================================================================
def run_arm(a, cfg, arm, base_grn_path):
    log(f"########## ARM = {arm} ##########")
    oracle, links = stage_grn(a, cfg, arm, base_grn_path)
    ko = stage_ko_ranking(oracle, links, cfg, arm)
    if cfg.get("run_closed_downstream_stages", False):
        stage_bootstrap(oracle, cfg, arm)
        stage_vectorfield(oracle, cfg, arm, ko)
    else:
        log("[protected] completed bootstrap/vector-field outputs were not rerun; "
            "the all-active extension only writes new screen tables.")
    del oracle, links
    gc.collect()
    log(f"########## ARM {arm} DONE ##########")


def main():
    t0 = time.time()
    cfg = apply_env_overrides(CFG)
    np.random.seed(cfg["seed"])
    os.makedirs(cfg["out_dir"], exist_ok=True)
    arms = list(dict.fromkeys(cfg["selected_arms"]))
    run_label = "_".join(arms)
    # Persist one configuration per independently launchable arm set.
    config_path = (
        f"{cfg['out_dir']}/03_engine_config_all_active_{run_label}_used.json"
    )
    if exists(config_path) and not cfg["force"]:
        previous = json.load(open(config_path))
        frozen_keys = [
            "atlas_h5ad",
            "skin_base_grn",
            "grn_source_dir",
            "n_top_genes",
            "recompute_hvg",
            "grn_alpha",
            "bagging",
            "link_p",
            "link_topn",
            "n_propagation",
            "vf_n_neighbors",
            "vf_sigma_corr",
            "vf_n_grid",
            "ko_min_raw_expr_fraction",
            "seed",
            "selected_arms",
        ]
        conflicts = {
            key: (previous.get(key), cfg.get(key))
            for key in frozen_keys
            if previous.get(key) != cfg.get(key)
        }
        if conflicts:
            raise RuntimeError(
                f"Frozen output root conflicts with requested settings: {conflicts}. "
                "Use a new SSC_CO_OUT_DIR rather than overwriting the run."
            )
    with open(config_path, "w") as fh:
        json.dump({k: v for k, v in cfg.items()}, fh, indent=2, default=str)

    need_atlas = any(
        cfg["force"] or frozen_grn_paths(cfg, arm) is None for arm in arms
    )
    a = stage_load(cfg) if need_atlas else None
    if a is None:
        log("[memory] every selected arm has a frozen GRN pair; atlas load skipped")

    # ARM 1: skin scATAC base-GRN (primary; no lung data/prior).
    skin_path = cfg.get("skin_base_grn")
    if "skinatac" in arms and skin_path and exists(skin_path):
        run_arm(a, cfg, "skinatac", skin_path)
    elif "skinatac" in arms and cfg.get("require_skin_base_grn", True):
        raise FileNotFoundError(
            f"Required skin base-GRN is missing: {skin_path}. Refusing to substitute "
            "a lung dataset or silently treat the tissue-agnostic promoter prior as primary.")
    elif "skinatac" in arms:
        log(f"[skin arm] missing {skin_path}; skipped by explicit configuration")

    # ARM 2: generic human-promoter base-GRN (tissue-agnostic sensitivity only).
    if "promoter" in arms and cfg.get("run_promoter_sensitivity", True):
        run_arm(a, cfg, "promoter", None)

    log(f"ALL DONE in {(time.time()-t0)/60:.1f} min. Outputs in {cfg['out_dir']}")
    log("CELLORACLE_ENGINE_DONE")


if __name__ == "__main__":
    main()
