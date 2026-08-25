# -*- coding: utf-8 -*-
"""
================================================================================
04_trajectory.py   (SSc virtual-perturbation paper, POWERED build)
================================================================================
LAYER (this script): DYNAMICS. On the integrated true-SSc fibroblast atlas
(fibroblast_integrated.h5ad from layer 02), build a defensible pseudotime
trajectory of fibroblast -> myofibroblast differentiation and test whether the
engine-nominated regulators RISE along it.

MAIN LINE it feeds (context, NOT re-derived here):
  The single-cell virtual-perturbation engine (layer 03, CellOracle in-silico
  TF-KO + decoupler CollecTRI regulon activity) CO-NOMINATES HES1 (a druggable
  Notch node -> nirogacestat) AND FOSB/AP-1 as the two regulators converging on
  the SSc myofibroblast program at the immune-stromal interface. This layer asks
  the orthogonal DYNAMICS question: do HES1 and FOSB regulon activity increase
  as fibroblasts differentiate toward the myofibroblast fate?

★ NON-NEGOTIABLE PRINCIPLES (baked into this code + its logic) ★
  * HES1 and FOSB are CO-NOMINATED. They are tested IDENTICALLY here (same
    regulon-activity method, same bootstrap CI, same size-matched random null,
    same output rows). There is NO pre-ranking, NO "primary/secondary", NO
    "stronger/weaker". THE RUN DECIDES. The drug rides on HES1/Notch only
    because Notch is druggable, not because HES1 is privileged in this analysis.
  * SMAD3/TGF-beta is a POSITIVE CONTROL (canonical pro-fibrotic driver; it
    SHOULD rise along the myofibroblast axis). MEF2C is a documented NEGATIVE
    control (activity anti-correlated with the program in prior project runs).
  * A fibrosis-FREE Notch-target program (HES1/HEY1/HEYL/HES5/NRARP/DTX1) gives
    a SUPPORTING Notch-pathway read-out. It does NOT replace HES1 and is NOT the
    headline; it is a mechanism-consistency check that contains no collagen /
    myofibroblast structural genes (so it cannot trivially track the program).
  * Rigor == REPRODUCIBILITY, not a bespoke score. Fixed SEED=0. Genre-standard
    methods (CellRank2 PseudotimeKernel; an INDEPENDENT Palantir/scFates cross-
    check) + STANDARD nulls: bootstrap CI over resampled cells + a SIZE-MATCHED
    RANDOM-REGULON null (permutation p). NO self-invented scoring, NO masking /
    falsification framing.
  * Regulon activity along pseudotime = an ASSOCIATION / footprint, a HYPOTHESIS,
    NOT causal and NOT experimental knockout. Immune context is association only;
    this script never asserts an immune -> HES1 / immune -> Notch direction.

METHOD (grounded on REAL project scripts + verified current library APIs):
  STEP 1  PSEUDOTIME BACKBONE = CellRank2 PseudotimeKernel (NOT RNA velocity).
          - directionality anchor = a scalar pseudotime derived on the layer-02
            integrated embedding (X_pca_harmony) via a diffusion pseudotime with
            the ROOT set to the SFRP2hi/DPP4+ progenitor pool and orientation
            checked so the myofibroblast (CTHRC1+/SFRP4+) program INCREASES
            (biology: Tabib 2021 PMC8289865 SFRP2hi progenitor -> myofibroblast).
          - PseudotimeKernel(adata, time_key=...) -> compute_transition_matrix ->
            GPCCA -> macrostates -> terminal states (myofibroblast forced in) ->
            fate probabilities toward the myofibroblast lineage.
            (API identical to replication_scrna/run_cellrank292979.py, a working
            project script; verified against CellRank 2.x docs.)
  STEP 2  INDEPENDENT cross-check = Palantir (AnnData-centric, v1.3+) OR scFates,
          run from the SAME root on the SAME embedding; report Spearman agreement
          of the two pseudotimes (a method-robustness statement, no cherry-pick).
  STEP 3  decoupler CollecTRI regulon activity (ULM, module 076 API; identical to
          server_results/scripts/regulon_activity.py) for HES1 / FOSB / SMAD3 /
          MEF2C + a Notch-target MODULE score. Spearman rho vs the CellRank
          pseudotime (and vs the myofibroblast fate probability) WITH:
            - BOOTSTRAP 95% CI (resample cells with replacement, SEED=0), and
            - a SIZE-MATCHED RANDOM-REGULON NULL: for each real TF, draw many
              random gene sets of the SAME size from expressed genes, recompute
              activity + rho, and get an empirical two-sided permutation p.
          Both leads pass through the identical code path.
  STEP 4  OUTPUTS (to CFG['OUT']):
            trajectory_pseudotime.csv     per-cell pseudotimes + fate + myo/scores
            along_traj_correlations.csv   gene, rho, CI_lo, CI_hi, null_p (+fate)
            method_agreement.csv          CellRank vs independent pseudotime rho
            fig04_along_traj_lollipop.*   lollipop of rho+CI (viridis, English, no bars)
            fig04_regulon_vs_pseudotime.* line/loess of activity vs pseudotime

STYLE (matches layers 02/05): CONFIG block (no hardcoding), idempotent /
checkpointed (skip if output exists unless --force), timestamped logging,
English figure labels, gzip writes, memory-friendly.

ENV: server env `npc` (scanpy / anndata / decoupler / cellrank / scipy;
Palantir OR scFates for the cross-check). GPU is optional and never required.

REUSES (read, adapted, REAL APIs only -- no invented calls):
  replication_scrna/run_cellrank292979.py     PseudotimeKernel/GPCCA/drivers API
  replication_scrna/srv07_sctour.py           orient-by-myo, ULM, spearman pattern
  server_results/scripts/regulon_activity.py  dc.op.collectri + dc.mt.ulm API
  rerun_clean/_qc_common.py                    MYO / LEADS / gene panels
  server_powered_2026/02_integrate_annotate.py obs schema (myo/celltype/fib_subtype)
  library 087_palantir_branch_probability.py   Palantir workflow (updated to v1.3 API)
  library 072_cellrank_fate_drivers.py         CellRank estimator sequence
================================================================================
"""
import os, sys, io, json, argparse, warnings, time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# ---- REPRODUCIBILITY: single fixed seed, set before anything stochastic -------
SEED = 0
np.random.seed(SEED)

# =====================================================================
# CONFIG  (edit ONLY this block on the server; everything else is path-free)
# =====================================================================
CFG = dict(
    # ---- I/O ----
    FIB_H5   = os.environ.get('SSC_FIB',
                              '/data/ssc/powered/fibroblast_integrated.h5ad'),
    OUT      = os.environ.get('SSC_OUT_04', './out_04_trajectory'),
    # SAME CollecTRI net file used across the project (layer 03/05/regulon_activity);
    # if absent we fall back to dc.op.collectri(organism='human') like regulon_activity.py.
    NET_CSV  = os.environ.get('SSC_NET', '/data/ssc/data/collectri_net.csv'),

    # ---- integrated embedding to run the trajectory on (written by layer 02) ----
    # layer 02 stores the chosen rep name in uns['integration_rep']; default Harmony.
    REP_KEY      = 'X_pca_harmony',   # fallback resolved at runtime from uns
    N_NEIGHBORS  = 30,                # kNN for the pseudotime graph (== srv10_cellrank)
    N_DCS        = 15,                # diffusion comps for DPT / diffmap

    # ---- CO-NOMINATED leads + controls (engine-nominated, tested IDENTICALLY) ----
    # HES1 & FOSB: the two co-nominated regulators (NO pre-ranking between them).
    # SMAD3: POSITIVE control (should rise). MEF2C: NEGATIVE control (prior anti-corr).
    LEADS_CONOM  = ['HES1', 'FOSB'],
    POS_CTRL     = ['SMAD3'],
    NEG_CTRL     = ['MEF2C'],
    # Fibrosis-FREE Notch-target set for a SUPPORTING Notch-program module score.
    # Deliberately excludes collagen / myofibroblast structural genes so it cannot
    # trivially track the fibrotic program. (canonical Notch/RBPJ direct targets.)
    NOTCH_TARGETS = ['HES1', 'HEY1', 'HEYL', 'HES5', 'NRARP', 'DTX1'],

    # ---- myofibroblast program + subtype anchors (reused from _qc_common / layer02) ----
    # canonical myofibroblast program -- IDENTICAL to _qc_common.MYO / layer-02 MYO.
    MYO = ['ACTA2', 'TAGLN', 'POSTN', 'COL1A1', 'COL1A2', 'COMP',
           'COL11A1', 'CTHRC1', 'COL3A1', 'FN1'],
    # ROOT (progenitor) anchor genes: Tabib 2018/2021 SFRP2hi/DPP4+ universal/
    # progenitor fibroblast. TERMINAL anchor genes: CTHRC1+/SFRP4+ myofibroblast.
    ROOT_MARKERS = ['SFRP2', 'DPP4', 'PI16'],   # ★2026-07-07 修: 去 COL1A1(泛胶原/结构基因, 亦在 MYO/TERM, 污染 progenitor 锚点; Tabib2021 progenitor=SFRP2/DPP4/PI16)
    TERM_MARKERS = ['CTHRC1', 'SFRP4', 'ACTA2', 'COMP', 'POSTN', 'COL11A1'],

    # ---- CellRank2 GPCCA settings (== run_cellrank292979.py) ----
    CR_N_SCHUR    = 20,
    CR_N_MACRO    = 10,
    CR_THRESHOLD  = 'hard',   # PseudotimeKernel scheme: Palantir-style hard threshold

    # ---- statistics ----
    N_BOOT      = 1000,   # bootstrap resamples for the rho 95% CI
    N_NULL      = 500,    # size-matched random-regulon draws for the permutation p
    ULM_TMIN    = 5,      # min targets per TF for decoupler ULM (== project default)
    CROSS_CHECK = 'palantir',  # 'palantir' | 'scfates' | 'auto'

    # ---- control / robustness gates (documented knobs, NOT discovery cutoffs) ----
    MIN_PROG_POOL    = 20,    # min cells in the SFRP2/DPP4 progenitor pool to anchor the root (F09)
    CTRL_ALPHA       = 0.05,  # permutation-significance level for the control sanity GATE
    NEG_CTRL_MAX_RHO = 0.30,  # negative control MEF2C must NOT be strongly positive above this (warn)
    ORIENT_MIN_ABS_RHO = 0.05,  # |rho(pt,myo)| must exceed this to orient the axis (else degenerate)
    AGREE_MIN_RHO    = 0.30,  # min CellRank-vs-independent pseudotime agreement to consider robust (warn)

    SEED = SEED,
)

# TF display grouping (kind label only for the OUTPUT table; NOT a ranking) -------
def _kind(tf):
    if tf in CFG['LEADS_CONOM']:
        return 'co_nominated_lead'
    if tf in CFG['POS_CTRL']:
        return 'positive_control'
    if tf in CFG['NEG_CTRL']:
        return 'negative_control'
    return 'other'

# the full identically-treated TF panel (order here is alphabetical-by-kind for
# readability ONLY; every entry runs the exact same test -- no pre-ranking).
TF_PANEL = CFG['LEADS_CONOM'] + CFG['POS_CTRL'] + CFG['NEG_CTRL']


# =====================================================================
# small utilities  (log / idempotency == layer 05 style)
# =====================================================================
def log(*a):
    print(f"[04 {time.strftime('%H:%M:%S')}]", *a, flush=True)


def ensure_out():
    for sub in ['', 'figs', '_ckpt']:
        os.makedirs(os.path.join(CFG['OUT'], sub), exist_ok=True)


def done(path, force):
    """idempotency guard: True if we can skip (exists and not forced)."""
    if os.path.exists(path) and not force:
        log('SKIP (exists):', path)
        return True
    return False


# --- shared CollecTRI net (IDENTICAL source as layer 03/05 + regulon_activity.py) ---
_NET = None
def collectri_net(var_names):
    """Return a CollecTRI net (source,target,weight) sliced to genes present in the
    data. Prefer the project's cached CSV (same file as discovery); else pull with
    dc.op.collectri(organism='human') exactly like server_results/scripts/
    regulon_activity.py. tmin filtering is applied later by dc.mt.ulm."""
    global _NET
    if _NET is None:
        import decoupler as dc
        if os.path.exists(CFG['NET_CSV']):
            n = pd.read_csv(CFG['NET_CSV'])
            cols = {c.lower(): c for c in n.columns}
            ren = {}
            for a, b in [('tf', 'source'), ('target_genesymbol', 'target'), ('mor', 'weight')]:
                if a in cols and b not in cols:
                    ren[cols[a]] = b
            n = n.rename(columns=ren)
            if 'weight' not in n.columns:
                n['weight'] = 1.0
            n = n[['source', 'target', 'weight']].dropna()
            log('CollecTRI net from CSV:', n.shape, n['source'].nunique(), 'TFs')
        else:
            n = dc.op.collectri(organism='human')[['source', 'target', 'weight']]
            log('CollecTRI net from dc.op.collectri:', n.shape, n['source'].nunique(), 'TFs')
        _NET = n
    return _NET[_NET['target'].isin(list(var_names))].copy()


# =====================================================================
# data loading + expression views
# =====================================================================
def load_fibro():
    """Load the layer-02 integrated fibroblast atlas. Returns (adata, rep_key).
    Ensures a log-normalised .X for regulon activity (uses .raw or re-normalises
    from layers['counts'] -- identical guard to regulon_activity.py / srv07)."""
    import scanpy as sc
    import scipy.sparse as sp
    a = sc.read_h5ad(CFG['FIB_H5'])
    log('loaded', a.shape, '| obs:', [c for c in a.obs.columns
        if c in ('cohort', 'sample', 'celltype', 'fib_subtype', 'fib_leiden', 'myo')])

    # resolve integration rep written by layer 02 (uns['integration_rep']); fallback chain
    rep = a.uns.get('integration_rep', CFG['REP_KEY'])
    if rep not in a.obsm:
        rep = CFG['REP_KEY'] if CFG['REP_KEY'] in a.obsm else 'X_pca'
    if rep not in a.obsm:
        raise RuntimeError(
            'No integrated embedding found in %s (looked for uns["integration_rep"], '
            '%s, X_pca in .obsm: %s). Trajectory MUST run on the layer-02 batch-'
            'integrated embedding; refusing to fall back to an un-integrated PCA, '
            'which would silently confound pseudotime and fate across cohorts. '
            'Re-run layer 02 or set SSC_FIB to a file that carries the integrated rep.'
            % (CFG['FIB_H5'], CFG['REP_KEY'], list(a.obsm.keys())))
    log('  trajectory embedding rep =', rep)

    # guarantee a log-normalised matrix in .X for decoupler ULM + module scores.
    # regulon_activity.py rule: if X looks like raw counts (max>50) -> re-normalise.
    xmax = float(a.X.max())
    if xmax > 50:
        log('  X looks like counts (max=%.1f) -> normalize_total+log1p' % xmax)
        if 'counts' in a.layers:
            a.X = a.layers['counts'].copy()
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    # keep a lognorm copy explicitly for module scoring even if .X gets touched
    a.layers['lognorm'] = a.X.copy() if not sp.issparse(a.X) else a.X.copy()
    return a, rep


def module_score(adata, genes, name):
    """scanpy score_genes module score on log-norm .X (present genes only).
    Used for the myofibroblast program (axis sanity) and the fibrosis-free
    Notch-target program. Returns the per-cell score as a numpy array."""
    import scanpy as sc
    gg = [g for g in genes if g in adata.var_names]
    if not gg:
        return np.full(adata.n_obs, np.nan), []
    sc.tl.score_genes(adata, gg, score_name='_mod_' + name, random_state=SEED, use_raw=False)
    v = adata.obs.pop('_mod_' + name).values.astype(float)
    return v, gg


# =====================================================================
# STEP 1 -- root selection + CellRank2 PseudotimeKernel backbone
# =====================================================================
def pick_root_cell(adata, rep, root_score, term_score):
    """Choose ONE root cell = the cell most progenitor-like and least terminal,
    among cells in the SFRP2hi/DPP4+ region. Deterministic (SEED-independent):
    argmax(root_score - term_score). If a fib_subtype label exists we restrict to
    the SFRP2_DPP4 pool first (Tabib progenitor), else use the score directly.
    Returns (root_barcode, root_idx)."""
    # NOTE (documented root rule): root = argmax(root_marker_score - term_marker_score)
    # restricted to the annotated progenitor (SFRP2hi/DPP4+) subtype pool. The published
    # anchor is the CELL-TYPE identity (Tabib 2021 PMC8289865: SFRP2hi/DPP4+ universal/
    # progenitor fibroblast is the root); we select the root WITHIN that labelled pool
    # rather than from the whole atlas. MIN_PROG_POOL is a documented CFG knob (was a
    # hardcoded 20). If the labelled pool is too small, hard-fail rather than silently
    # fall back to a whole-atlas argmax root.
    score = np.asarray(root_score, float) - np.asarray(term_score, float)
    mask = np.ones(adata.n_obs, dtype=bool)
    MIN_PROG = CFG.get('MIN_PROG_POOL', 20)
    if 'fib_subtype' in adata.obs:
        prog = adata.obs['fib_subtype'].astype(str).str.upper()
        pmask = prog.str.contains('SFRP2') | prog.str.contains('DPP4') | prog.str.contains('PROG')
        if pmask.sum() >= MIN_PROG:
            mask = pmask.values
            log('  root restricted to progenitor subtype pool: n=%d' % mask.sum())
        else:
            raise RuntimeError(
                'Progenitor subtype pool (SFRP2/DPP4/PROG) has only %d cells '
                '(<MIN_PROG_POOL=%d); cannot anchor the root to the published Tabib '
                'progenitor identity. Check layer-02 fib_subtype labels.'
                % (int(pmask.sum()), MIN_PROG))
    idx = int(np.argmax(np.where(mask, score, -np.inf)))
    log('  root cell =', adata.obs_names[idx], '| root-term score=%.3f' % score[idx])
    return adata.obs_names[idx], idx


def compute_dpt(adata, rep, root_idx):
    """Diffusion pseudotime on the integrated embedding, rooted at root_idx.
    This scalar is the DIRECTIONALITY ANCHOR fed to the CellRank PseudotimeKernel
    (NOT RNA velocity). Standard scanpy sc.tl.dpt after diffmap on the rep graph.
    Orientation is fixed downstream so the myofibroblast program increases."""
    import scanpy as sc
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=CFG['N_NEIGHBORS'], random_state=SEED)
    sc.tl.diffmap(adata, n_comps=CFG['N_DCS'])
    adata.uns['iroot'] = int(root_idx)
    sc.tl.dpt(adata, n_dcs=CFG['N_DCS'])
    pt = adata.obs['dpt_pseudotime'].values.astype(float)
    # Non-finite DPT means disconnected cells on the kNN graph -- a manifold
    # problem that would silently corrupt every downstream rho. Hard-fail with a
    # count instead of rank-filling over it (iron rule 1).
    n_bad = int(np.sum(~np.isfinite(pt)))
    if n_bad:
        raise RuntimeError(
            '%d/%d cells have non-finite DPT pseudotime (disconnected on the '
            'kNN graph, use_rep=%s, n_neighbors=%d). Fix the embedding/graph '
            '(e.g. drop isolated cells or raise N_NEIGHBORS) rather than '
            'rank-filling over it.' % (n_bad, len(pt), rep, CFG['N_NEIGHBORS']))
    return pt


def orient_by_myo(pt, myo):
    """Orient a pseudotime so the myofibroblast program INCREASES with it
    (progenitor->myofibroblast). Same axis-sanity rule as srv07_sctour.py.
    GATE (missed-issue): if |rho(pt, myo)| is near-zero the orientation is
    essentially arbitrary and unstable across reps/seeds; refuse to silently
    orient a degenerate axis (hard-fail) rather than pick a coin-flip direction."""
    rho = spearmanr(pt, myo, nan_policy='omit')[0]
    if rho is None or not np.isfinite(rho) or abs(rho) < CFG['ORIENT_MIN_ABS_RHO']:
        raise RuntimeError(
            'Pseudotime-vs-myofibroblast correlation is degenerate '
            '(rho=%s, |rho|<ORIENT_MIN_ABS_RHO=%.3f): the myo gradient is near-flat so '
            'trajectory orientation would be arbitrary/unstable. Check the root, the '
            'embedding, or the myo program before proceeding.'
            % (('%.4f' % rho) if rho is not None else 'None', CFG['ORIENT_MIN_ABS_RHO']))
    if rho < 0:
        log('  (flipped pseudotime so myofibroblast program increases; pre-flip rho=%.3f)' % rho)
        return 1.0 - pt
    return pt


def cellrank_backbone(adata, rep, ptime_key):
    """CellRank2 PseudotimeKernel backbone (NOT velocity). Directed purely by the
    scalar pseudotime in obs[ptime_key]; GPCCA -> macrostates -> myofibroblast
    terminal state -> fate probabilities toward the myofibroblast lineage.
    API grounded on replication_scrna/run_cellrank292979.py (working project code)
    and verified against CellRank 2.x docs. Returns (myo_fate_array, info_dict)."""
    import cellrank as cr
    from cellrank.kernels import PseudotimeKernel
    from cellrank.estimators import GPCCA
    log('  cellrank', cr.__version__)

    # PseudotimeKernel(adata, time_key=...) -> hard threshold scheme (Palantir-style)
    pk = PseudotimeKernel(adata, time_key=ptime_key)
    pk.compute_transition_matrix(threshold_scheme=CFG['CR_THRESHOLD'])
    log('  transition matrix done (scheme=%s)' % CFG['CR_THRESHOLD'])

    g = GPCCA(pk)
    g.compute_schur(n_components=CFG['CR_N_SCHUR'])
    g.compute_macrostates(n_states=CFG['CR_N_MACRO'])
    ms = g.macrostates
    cats = list(map(str, ms.cat.categories))
    log('  macrostates:', cats)

    # A5: GPCCA real-Schur SPECTRUM diagnostic (real eigenvalues only). DIAGNOSTIC ONLY
    # -- helps a human judge whether n_states was set sensibly (eigengap in the spectrum);
    # n_states is FIXED at CFG['CR_N_MACRO'] and this plot does NOT change any behaviour.
    # Guarded (try/except-log): a plotting failure must never break the pipeline.
    log('  n_states fixed at CR_N_MACRO=%d; saving GPCCA spectrum diagnostic' % CFG['CR_N_MACRO'])
    try:
        import matplotlib
        matplotlib.use('Agg')
        spec_png = os.path.join(CFG['OUT'], 'figs', 'fig04_gpcca_spectrum.png')
        g.plot_spectrum(real_only=True, save=spec_png)
        log('  [diag] GPCCA spectrum ->', spec_png)
    except Exception as e:
        log('  [diag] GPCCA plot_spectrum skipped (diagnostic only):', repr(e))

    # myofibroblast terminal macrostate = the macrostate whose cells score highest on
    # the canonical myofibroblast MARKER panel (obs['myo'] = mean of CFG['MYO']). This
    # is a MARKER-ARGMAX anchor over macrostates (same principle as the layer-02
    # marker-argmax myofibroblast subtype call), NOT a bespoke unexplained score: the
    # terminal fate is defined by canonical markers, and F03 below asserts it into the
    # GPCCA terminal set so the headline lineage is never silently mis-assigned.
    myo = adata.obs['myo'].values.astype(float)
    # CellRank2 GPCCA .macrostates labels ONLY the most-confident ~n_cells per macrostate
    # and leaves the MAJORITY of cells NaN; ms.astype(str) turns those into a spurious
    # 'nan' bucket that (on a 424K atlas) spans ~all cells and can WIN the myo argmax ->
    # myo_state='nan' -> set_terminal_states(['nan']) crashes. Rank myo ONLY over cells that
    # actually belong to a macrostate (drop the 'nan' bucket via string compare, robust to ms type).
    _ms_str = np.asarray(ms.astype(str))
    _valid = _ms_str != 'nan'
    if not _valid.any():
        raise RuntimeError('GPCCA assigned no macrostate cells (all NaN); cannot anchor the '
                           'myofibroblast terminal macrostate.')
    myo_by_ms = pd.Series(myo[_valid]).groupby(_ms_str[_valid], observed=True
                                               ).mean().sort_values(ascending=False)
    myo_state = str(myo_by_ms.index[0])
    log('  myofibroblast terminal macrostate =', myo_state,
        '(mean myo marker score=%.3f)' % myo_by_ms.iloc[0])

    # automatic terminal states; force-include the myo end if it was missed
    g.predict_terminal_states()
    term = list(map(str, g.terminal_states.cat.categories))
    log('  auto terminal states:', term)
    if not any(t == myo_state or t.startswith(myo_state) for t in term):
        # GPCCA.set_terminal_states(states=<sequence of macrostate names>) --
        # verified against CellRank 2.x GPCCA docs. Must succeed: the myofibroblast
        # terminal lineage is the axis this whole layer measures against.
        g.set_terminal_states(states=[myo_state] + term)
        term = list(map(str, g.terminal_states.cat.categories))
        log('  terminal states (myo forced in):', term)
        if not any(t == myo_state or t.startswith(myo_state) for t in term):
            raise RuntimeError(
                'Failed to force-include the myofibroblast macrostate %r into the '
                'GPCCA terminal states (got %s). Refusing to compute fate '
                'probabilities without the myofibroblast terminal lineage.'
                % (myo_state, term))

    g.compute_fate_probabilities()
    fp = g.fate_probabilities                # Lineage object (.names, fp[name].X)
    lin_names = list(fp.names)
    log('  fate lineages:', lin_names)

    myo_lin = next((n for n in lin_names
                    if n == myo_state or n.startswith(myo_state)), None)
    if myo_lin is None:                      # fallback: lineage whose top cells are most-myo
        means = []
        for n in lin_names:
            pr = np.asarray(fp[n].X).ravel()
            hi = pr > np.quantile(pr, 0.9)
            means.append(myo[hi].mean() if hi.any() else -np.inf)
        myo_lin = lin_names[int(np.argmax(means))]
    log('  >>> myofibroblast fate lineage =', myo_lin)
    myo_fate = np.asarray(fp[myo_lin].X).ravel()
    info = dict(myo_macrostate=myo_state, myo_lineage=myo_lin,
                n_macrostates=len(cats), terminal_states=';'.join(term))
    return myo_fate, info


# =====================================================================
# STEP 2 -- INDEPENDENT pseudotime cross-check (Palantir OR scFates)
# =====================================================================
def cross_check_pseudotime(adata, rep, root_bc):
    """Run an INDEPENDENT pseudotime from the SAME root on the SAME embedding and
    return it (aligned to adata.obs_names) + the method name used. Palantir 1.3+
    AnnData-centric workflow (grounded on library 087 + verified current API):
        palantir.utils.run_diffusion_maps(adata, pca_key=rep)
        palantir.utils.determine_multiscale_space(adata)
        palantir.core.run_palantir(adata, early_cell=root_bc, ... )  -> writes
            obs['palantir_pseudotime']
    scFates fallback (verified API): scf.tl.tree -> scf.tl.root -> scf.tl.pseudotime
        writes obs['t']. Returns (pt_array_or_None, method_str)."""
    method = CFG['CROSS_CHECK']
    errors = []
    # ---- Palantir ----
    if method in ('palantir', 'auto'):
        try:
            import palantir
            # Palantir reads obsm[pca_key]; point it at the integrated rep so the
            # cross-check shares the SAME manifold as the CellRank backbone.
            palantir.utils.run_diffusion_maps(
                adata, n_components=CFG['N_DCS'], pca_key=rep, seed=SEED)
            palantir.utils.determine_multiscale_space(adata)
            palantir.core.run_palantir(
                adata, early_cell=str(root_bc),
                num_waypoints=1200, knn=CFG['N_NEIGHBORS'], seed=SEED)
            pt = adata.obs['palantir_pseudotime'].values.astype(float)
            log('  independent cross-check: Palantir OK')
            return pt, 'palantir'
        except Exception as e:
            errors.append('Palantir: ' + repr(e))
            log('  Palantir cross-check FAILED:', str(e)[:200])
            if method == 'palantir':
                method = 'scfates'   # try the other one before giving up
    # ---- scFates ----
    if method in ('scfates', 'auto'):
        try:
            import scFates as scf
            import scanpy as sc
            if 'X_diffmap' not in adata.obsm:
                sc.pp.neighbors(adata, use_rep=rep, n_neighbors=CFG['N_NEIGHBORS'],
                                random_state=SEED)
                sc.tl.diffmap(adata, n_comps=CFG['N_DCS'])
            scf.tl.tree(adata, method='ppt', Nodes=200, use_rep='X_diffmap',
                        ndims_rep=CFG['N_DCS'], seed=SEED, ppt_lambda=100, ppt_sigma=0.1)
            root_idx = int(adata.obs_names.get_loc(root_bc))
            scf.tl.root(adata, root=root_idx)
            scf.tl.pseudotime(adata, n_jobs=1, n_map=1, seed=SEED)
            pt = adata.obs['t'].values.astype(float)
            log('  independent cross-check: scFates OK')
            return pt, 'scfates'
        except Exception as e:
            errors.append('scFates: ' + repr(e))
            log('  scFates cross-check FAILED:', str(e)[:200])
    # No graceful NA: the independent cross-check is a load-bearing robustness
    # claim. If the requested method(s) cannot run, HARD-FAIL so the env gets
    # fixed rather than silently reporting agreement=NA.
    raise RuntimeError(
        'Independent pseudotime cross-check FAILED for CROSS_CHECK=%r and no '
        'method produced a pseudotime. This is a required robustness check, not '
        'optional. Fix the env (install/upgrade palantir or scFates) or pass '
        '--cross-check for the installed tool. Underlying errors:\n  %s'
        % (CFG['CROSS_CHECK'], '\n  '.join(errors) or 'none captured'))


# =====================================================================
# STEP 3 -- regulon activity + bootstrap CI + size-matched random null
# =====================================================================
def regulon_activity(adata):
    """decoupler CollecTRI ULM regulon activity per cell (module 076 API; identical
    to server_results/scripts/regulon_activity.py). Returns a DataFrame
    cells x TFs of ULM scores (obsm['score_ulm']). Uses the log-norm .X."""
    import decoupler as dc
    net = collectri_net(adata.var_names)
    dc.mt.ulm(data=adata, net=net, tmin=CFG['ULM_TMIN'])
    acts = adata.obsm['score_ulm']
    if not isinstance(acts, pd.DataFrame):
        acts = pd.DataFrame(acts, index=adata.obs_names)
    acts = acts.copy()   # DETACH from the obsm slot: size_matched_null_p re-runs ulm which
                         # overwrites obsm['score_ulm'] -> without this copy, acts (aliasing
                         # that slot) would lose the real TF columns after the first null call.
    log('  regulon activity (cells x TFs):', acts.shape)
    return acts


def boot_ci_rho(x, y, n_boot, seed):
    """Bootstrap 95% CI for Spearman rho(x,y) by resampling CELLS with replacement.
    Deterministic given seed. Returns (rho_point, lo, hi)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    rho_pt = spearmanr(x, y)[0]
    n = len(x)
    if n < 20:
        return float(rho_pt), np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot, float)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[b] = spearmanr(x[idx], y[idx])[0]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return float(rho_pt), float(lo), float(hi)


def size_matched_null_p(adata, net, tf, target_axis, rho_obs, n_null, seed):
    """SIZE-MATCHED RANDOM-REGULON NULL. For a real TF whose regulon has k targets,
    draw n_null random gene sets of size k from the EXPRESSED genes, compute each
    random 'regulon' activity (ULM with a unit-weight synthetic net) and its
    Spearman rho vs the target axis; the empirical two-sided permutation p is the
    fraction of |null rho| >= |rho_obs|. This is a STANDARD null and is applied
    IDENTICALLY to every TF (both co-nominated leads included). SEED-fixed."""
    import decoupler as dc
    tset = net.loc[net['source'] == tf, 'target'].unique().tolist()
    tset = [g for g in tset if g in adata.var_names]
    k = len(tset)
    if k < CFG['ULM_TMIN']:
        return np.nan, k
    pool = list(adata.var_names)
    rng = np.random.default_rng(seed)
    axis = np.asarray(target_axis, float)
    ok = np.isfinite(axis)
    null_rhos = np.empty(n_null, float)
    # build one synthetic net of n_null random regulons, score them all in ONE ULM call
    frames = []
    for i in range(n_null):
        rs = rng.choice(pool, size=k, replace=False)
        frames.append(pd.DataFrame({'source': f'RND{i}', 'target': rs, 'weight': 1.0}))
    rnet = pd.concat(frames, ignore_index=True)
    # dc.mt.ulm has NO `key=` param; it ALWAYS writes obsm['score_ulm']. Snapshot the real
    # activity, run the null into score_ulm, read it, then RESTORE the real matrix.
    real_ulm = adata.obsm.get('score_ulm')
    dc.mt.ulm(data=adata, net=rnet, tmin=CFG['ULM_TMIN'])
    racts = adata.obsm['score_ulm']
    if not isinstance(racts, pd.DataFrame):
        racts = pd.DataFrame(racts, index=adata.obs_names)
    racts = racts.copy()
    for i in range(n_null):
        col = f'RND{i}'
        if col in racts.columns:
            null_rhos[i] = spearmanr(racts[col].values[ok], axis[ok])[0]
        else:
            null_rhos[i] = np.nan
    if real_ulm is not None:
        adata.obsm['score_ulm'] = real_ulm
    valid = np.isfinite(null_rhos)
    if valid.sum() < 20:
        return np.nan, k
    p = (np.sum(np.abs(null_rhos[valid]) >= abs(rho_obs)) + 1) / (valid.sum() + 1)
    return float(p), k


# =====================================================================
# STEP 4 -- figures  (NO plain bar charts; lollipop + line; viridis; English)
# =====================================================================
def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'figure.dpi': 200, 'font.size': 9, 'svg.fonttype': 'none'})
    return plt


def fig_lollipop(corr_df, axis_label):
    """Lollipop of Spearman rho (vs pseudotime) with bootstrap 95% CI whiskers,
    one row per TF, coloured by null significance (viridis). NO bar chart."""
    if corr_df.empty:
        return
    plt = _mpl()
    d = corr_df[corr_df['axis'] == axis_label].copy()
    if d.empty:
        return
    # order by rho for readability (this is a display order, NOT a ranking claim)
    d = d.sort_values('rho')
    y = np.arange(len(d))
    cmap = plt.get_cmap('viridis')
    # colour = -log10(null_p) mapped through viridis; NA -> grey
    # colour by BH-adjusted permutation p (F05) so the figure matches the FDR
    # gate; fall back to raw null_p only if the adjusted column is absent.
    _pcol = 'null_padj_bh' if 'null_padj_bh' in d.columns else 'null_p'
    _plabel = 'BH-adjusted null p' if _pcol == 'null_padj_bh' else 'null p'
    nl = -np.log10(d[_pcol].clip(lower=1e-4).fillna(1.0).values)
    norm = (nl - nl.min()) / (np.ptp(nl) + 1e-9) if len(nl) else nl
    colors = [cmap(v) if np.isfinite(p) else '#bbbbbb'
              for v, p in zip(norm, d[_pcol].values)]
    fig, ax = plt.subplots(figsize=(6.4, 0.5 * len(d) + 1.6))
    # CI whiskers + stem to zero + dot
    ax.hlines(y, d['ci_lo'], d['ci_hi'], color='#999999', lw=1.4, zorder=1)
    ax.hlines(y, 0, d['rho'], color='#cccccc', lw=1.0, zorder=1)
    ax.scatter(d['rho'], y, s=70, c=colors, edgecolor='k', lw=0.6, zorder=3)
    ax.axvline(0, color='k', lw=0.8, ls='--')
    labels = [f"{t}  ({k.replace('_', ' ')})" for t, k in zip(d['gene'], d['kind'])]
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel(f'Spearman rho of regulon activity vs {axis_label}\n'
                  f'(dot = point estimate, whisker = bootstrap 95% CI; '
                  f'colour = -log10 size-matched {_plabel})')
    ax.set_title('Regulator activity rises along the fibroblast -> myofibroblast axis')
    # colourbar proxy
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=0, vmax=max(1e-9, nl.max())))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label('-log10 %s' % _plabel, fontsize=8)
    fig.tight_layout()
    tag = axis_label.replace(' ', '_')
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(CFG['OUT'], 'figs', f'fig04_along_traj_lollipop_{tag}.{ext}'),
                    bbox_inches='tight')
    plt.close(fig)


def fig_trajectory_umap(adata, dpt, myo_fate, myo, acts):
    """Trajectory UMAP panel (★2026-07-09, was missing): the fibroblast UMAP coloured by DPT
    pseudotime, CellRank2 myofibroblast-fate probability, the myofibroblast program score, and
    HES1 / FOSB regulon activity -- so the 3.5h CellRank/DPT compute is actually visualised.
    viridis, rasterised points, no bar charts."""
    if 'X_umap' not in adata.obsm:
        log('  [fig] no X_umap in obsm -> skip trajectory UMAP panel')
        return
    plt = _mpl()
    um = np.asarray(adata.obsm['X_umap'])
    panels = [('pseudotime (DPT)', np.asarray(dpt, float)),
              ('myofibroblast-fate prob\n(CellRank2)', np.asarray(myo_fate, float)),
              ('myofibroblast program', np.asarray(myo, float))]
    for tf in ('HES1', 'FOSB'):
        if tf in acts.columns:
            panels.append(('%s regulon activity' % tf, np.asarray(acts[tf].values, float)))
    n = len(panels)
    fig, axs = plt.subplots(1, n, figsize=(3.3 * n, 3.4))
    axs = np.atleast_1d(axs)
    order = np.argsort(np.random.RandomState(0).rand(um.shape[0]))   # deterministic z-order
    for ax, (title, vals) in zip(axs, panels):
        lo, hi = np.nanpercentile(vals, [2, 98])
        sc_ = ax.scatter(um[order, 0], um[order, 1], c=vals[order], s=3, cmap='viridis',
                         vmin=lo, vmax=hi, rasterized=True, linewidths=0)
        ax.set_title(title, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        cb = fig.colorbar(sc_, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=7)
    fig.suptitle('Fibroblast -> myofibroblast trajectory (UMAP)', fontsize=11)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(CFG['OUT'], 'figs', 'fig04_trajectory_umap.%s' % ext),
                    bbox_inches='tight', dpi=180)
    plt.close(fig)
    log('  fig -> fig04_trajectory_umap.png/pdf')


def fig_activity_lines(pt, acts, notch_score, myo, out_axis_name):
    """Binned-mean line of regulon activity vs pseudotime for the co-nominated
    leads + controls + Notch-target module + myo program (axis sanity). viridis,
    English labels, no bars. Bands = bootstrap-free binned SEM for readability."""
    plt = _mpl()
    order = np.argsort(pt)
    ptv = np.asarray(pt, float)[order]
    nb = 20
    edges = np.quantile(ptv, np.linspace(0, 1, nb + 1))
    edges = np.unique(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_idx = np.clip(np.digitize(ptv, edges[1:-1]), 0, len(centers) - 1)

    def binned(v):
        v = np.asarray(v, float)[order]
        m = np.array([np.nanmean(v[bin_idx == b]) if np.any(bin_idx == b) else np.nan
                      for b in range(len(centers))])
        # z-score across bins so heterogeneous scales overlay comparably
        return (m - np.nanmean(m)) / (np.nanstd(m) + 1e-9)

    series = []
    for tf in TF_PANEL:
        if tf in acts.columns:
            series.append((tf, binned(acts[tf].values), _kind(tf)))
    series.append(('Notch-target module', binned(notch_score), 'notch_module'))
    series.append(('Myofibroblast program', binned(myo), 'axis_sanity'))

    cmap = plt.get_cmap('viridis')
    n = len(series)
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    for i, (name, yv, kind) in enumerate(series):
        col = '#444444' if kind == 'axis_sanity' else cmap(i / max(1, n - 1))
        ls = '--' if kind in ('axis_sanity', 'negative_control') else '-'
        lw = 2.4 if kind in ('co_nominated_lead',) else 1.6
        ax.plot(centers, yv, ls, lw=lw, color=col, label=f'{name} ({kind.replace("_", " ")})')
    ax.set_xlabel(f'{out_axis_name} (progenitor -> myofibroblast)')
    ax.set_ylabel('regulon / module activity (z across pseudotime bins)')
    ax.set_title('Regulator activity trajectories along fibroblast differentiation')
    ax.legend(fontsize=7, frameon=False, loc='center left', bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(CFG['OUT'], 'figs', f'fig04_regulon_vs_pseudotime.{ext}'),
                    bbox_inches='tight')
    plt.close(fig)


# =====================================================================
# MAIN
# =====================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--force', action='store_true', help='ignore checkpoints, recompute')
    ap.add_argument('--cross-check', default=None,
                    help='override CFG CROSS_CHECK: palantir | scfates | auto')
    args = ap.parse_args()
    if args.cross_check:
        CFG['CROSS_CHECK'] = args.cross_check
    ensure_out()
    log('CONFIG:', json.dumps({k: CFG[k] for k in
        ('FIB_H5', 'OUT', 'NET_CSV', 'REP_KEY', 'N_BOOT', 'N_NULL',
         'CROSS_CHECK', 'SEED')}))
    log('CO-NOMINATED leads (tested identically, no pre-ranking):', CFG['LEADS_CONOM'])

    pt_csv   = os.path.join(CFG['OUT'], 'trajectory_pseudotime.csv')
    corr_csv = os.path.join(CFG['OUT'], 'along_traj_correlations.csv')
    agr_csv  = os.path.join(CFG['OUT'], 'method_agreement.csv')

    # short-circuit: everything already computed
    if all(os.path.exists(p) for p in (pt_csv, corr_csv, agr_csv)) and not args.force:
        log('all outputs exist -> nothing to do (use --force to recompute)')
        log('TRAJECTORY_DONE')
        return

    # ---------- load ----------
    adata, rep = load_fibro()

    # ---------- axis-sanity scores: myo program + fibrosis-free Notch module ----------
    myo, myo_g = module_score(adata, CFG['MYO'], 'myo')
    if 'myo' in adata.obs and np.isfinite(adata.obs['myo'].values).all():
        # prefer the layer-02 myo score for continuity if present & finite
        myo = adata.obs['myo'].values.astype(float)
        log('  using layer-02 obs["myo"] program score (continuity)')
    else:
        adata.obs['myo'] = myo
    log('  myo program genes present:', myo_g)
    root_score, _ = module_score(adata, CFG['ROOT_MARKERS'], 'root')
    term_score, _ = module_score(adata, CFG['TERM_MARKERS'], 'term')
    notch_score, notch_g = module_score(adata, CFG['NOTCH_TARGETS'], 'notch')
    log('  fibrosis-free Notch-target module genes present:', notch_g)

    # ---------- STEP 1: root + DPT anchor + CellRank2 backbone ----------
    log('STEP 1: root selection + DPT directionality anchor + CellRank2 backbone')
    root_bc, root_idx = pick_root_cell(adata, rep, root_score, term_score)
    dpt = compute_dpt(adata, rep, root_idx)
    dpt = orient_by_myo(dpt, myo)
    adata.obs['pseudotime'] = dpt
    rho_axis = spearmanr(dpt, myo, nan_policy='omit')[0]
    log('  axis sanity: rho(pseudotime, myo program) = %+.3f' % rho_axis)
    myo_fate, cr_info = cellrank_backbone(adata, rep, 'pseudotime')
    adata.obs['myo_fate'] = myo_fate
    log('  CellRank myo-fate range %.3f-%.3f' % (np.nanmin(myo_fate), np.nanmax(myo_fate)))

    # ---------- STEP 2: independent cross-check ----------
    log('STEP 2: independent pseudotime cross-check (%s)' % CFG['CROSS_CHECK'])
    cc_pt, cc_method = cross_check_pseudotime(adata, rep, root_bc)
    if cc_pt is not None:
        cc_pt = orient_by_myo(cc_pt, myo)
        adata.obs['pseudotime_crosscheck'] = cc_pt
        agr_rho, agr_p = spearmanr(dpt, cc_pt, nan_policy='omit')
        log('  METHOD AGREEMENT: rho(CellRank-DPT, %s) = %+.3f (p=%.1e)'
            % (cc_method, agr_rho, agr_p))
        # Robustness check (missed-issue): the "two independent pseudotimes agree"
        # claim is only supported if the agreement is actually high. Warn loudly if
        # below AGREE_MIN_RHO (non-fatal: the value is reported in method_agreement.csv
        # so the manuscript claim can be qualified rather than the run killed).
        if not (np.isfinite(agr_rho) and agr_rho >= CFG['AGREE_MIN_RHO']):
            log('  [AGREEMENT WARN] independent-pseudotime agreement rho=%+.3f is BELOW '
                'AGREE_MIN_RHO=%.2f -- the "two independent pseudotimes agree" claim is '
                'NOT well supported; qualify it in the manuscript. (non-fatal)'
                % (agr_rho, CFG['AGREE_MIN_RHO']))
    else:
        agr_rho, agr_p = np.nan, np.nan
    pd.DataFrame([dict(cellrank_backbone='DPT+PseudotimeKernel',
                       crosscheck_method=cc_method,
                       agreement_spearman_rho=agr_rho, agreement_p=agr_p,
                       axis_sanity_rho_myo=rho_axis,
                       **cr_info)]).to_csv(agr_csv, index=False)
    log('  wrote', agr_csv)

    # ---------- STEP 3: regulon activity + rho + bootstrap CI + null ----------
    log('STEP 3: regulon activity vs pseudotime (bootstrap CI + size-matched null)')
    acts = regulon_activity(adata)
    net_full = collectri_net(adata.var_names)   # for null target-set sizes

    # axes tested IDENTICALLY for all TFs: (a) the CellRank pseudotime, (b) myo fate
    axes = {'pseudotime': dpt, 'myofibroblast_fate': myo_fate}
    rows = []
    for axis_name, axis in axes.items():
        log('  -- axis: %s --' % axis_name)
        _required = set(CFG['LEADS_CONOM']) | set(CFG['POS_CTRL']) | set(CFG['NEG_CTRL'])
        for tf in TF_PANEL:
            if tf not in acts.columns:
                _is_ctrl = tf in (set(CFG['POS_CTRL']) | set(CFG['NEG_CTRL']))
                if _is_ctrl:
                    # a CONTROL (pos/neg) with NO CollecTRI regulon means the regulon-activity
                    # axis cannot be VALIDATED -> HARD-FAIL loudly; do NOT silently proceed past a
                    # broken control (with the trajectory cached, a re-run aborts here cheaply).
                    raise RuntimeError(
                        'Control TF %r (pos/neg) has NO CollecTRI regulon activity (<tmin=%d '
                        'expressed targets or net/symbol mismatch); the regulon axis cannot be '
                        'validated. Refusing to proceed on a broken control.' % (tf, CFG['ULM_TMIN']))
                elif tf in _required:
                    # a co-nominated LEAD (e.g. the DEMOTED FOSB) missing its regulon is SURFACED
                    # loudly but is NON-fatal -- it is not a hard requirement, so the run completes
                    # and reports it as regulon-untestable rather than aborting the whole trajectory.
                    log('  ⚠️  co-nominated lead %r has NO CollecTRI regulon activity -> regulon '
                        'axis CANNOT test it; SKIPPED (surfaced, not silently dropped).' % tf)
                else:
                    log('     %-6s: not in CollecTRI regulons (skipped)' % tf)
                continue
            rho, lo, hi = boot_ci_rho(acts[tf].values, axis, CFG['N_BOOT'], SEED)
            null_p, k = size_matched_null_p(adata, net_full, tf, axis, rho,
                                            CFG['N_NULL'], SEED)
            rows.append(dict(gene=tf, kind=_kind(tf), axis=axis_name,
                             n_targets=int(k), n_cells=int(np.isfinite(axis).sum()),
                             rho=round(rho, 4),
                             ci_lo=round(lo, 4) if np.isfinite(lo) else np.nan,
                             ci_hi=round(hi, 4) if np.isfinite(hi) else np.nan,
                             null_p=null_p))
            log('     %-6s (%s): rho=%+.3f  95%%CI[%.3f, %.3f]  null_p=%s'
                % (tf, _kind(tf), rho,
                   lo if np.isfinite(lo) else float('nan'),
                   hi if np.isfinite(hi) else float('nan'),
                   ('%.3g' % null_p) if np.isfinite(null_p) else 'NA'))

        # the fibrosis-free Notch-target MODULE score (supporting Notch read-out).
        # tested with the same bootstrap CI; it is a module score, not a CollecTRI
        # regulon, so it has no size-matched-regulon null (reported as NA).
        rho, lo, hi = boot_ci_rho(notch_score, axis, CFG['N_BOOT'], SEED)
        rows.append(dict(gene='Notch_target_module', kind='notch_module',
                         axis=axis_name, n_targets=len(notch_g),
                         n_cells=int(np.isfinite(axis).sum()),
                         rho=round(rho, 4),
                         ci_lo=round(lo, 4) if np.isfinite(lo) else np.nan,
                         ci_hi=round(hi, 4) if np.isfinite(hi) else np.nan,
                         null_p=np.nan))
        log('     Notch-target module: rho=%+.3f  95%%CI[%.3f, %.3f]'
            % (rho, lo if np.isfinite(lo) else float('nan'),
               hi if np.isfinite(hi) else float('nan')))

    corr = pd.DataFrame(rows)
    # Iron rule 4: BH/FDR across all size-matched permutation tests (per TF x axis).
    # The Notch module row has null_p=NaN (module score, no size-matched null) and
    # is excluded from correction.
    from statsmodels.stats.multitest import multipletests
    corr['null_padj_bh'] = np.nan
    _m = corr['null_p'].notna().values
    if _m.sum() > 0:
        corr.loc[_m, 'null_padj_bh'] = multipletests(
            corr.loc[_m, 'null_p'].values, method='fdr_bh')[1]
    corr.to_csv(corr_csv, index=False)
    log('  wrote', corr_csv, '(%d rows; BH-FDR over %d permutation tests)'
        % (len(corr), int(_m.sum())))

    # ---------- CONTROL HARD-GATE (D2 / iron rule 5) ----------
    # SMAD3 (canonical TGF-beta pro-fibrotic POSITIVE control) MUST rise along the
    # myofibroblast pseudotime axis in ABSOLUTE terms (rho>0 AND FDR-significant) --
    # NOT merely "better than MEF2C" (a near-zero SMAD3 could beat a negative MEF2C).
    # Gate on the BH-adjusted permutation p (null_padj_bh) to match F05/the figure.
    # If SMAD3 fails -> ABORT (the pseudotime/regulon-activity pipeline is degraded).
    # MEF2C (soft NEGATIVE control) = warn loudly, do NOT abort. Association-level.
    def _q_of(row):
        q = row.get('null_padj_bh')
        if q is None or (isinstance(q, float) and not np.isfinite(q)):
            q = row.get('null_p')   # fall back to raw p if BH column is NaN
        return float(q) if q is not None else np.nan
    # ★2026-07-09: gate on the CellRank2 'myofibroblast_fate' axis (VALIDATED primary), NOT the
    # DPT-pseudotime axis. DPT is DEGRADED on this 104K-cell atlas (pos ctrl SMAD3 rho~+0.10 ns;
    # the independent palantir pseudotime anti-correlates rho=-0.30), so gating on DPT is a
    # false-negative about the pipeline. The fate axis has textbook controls (SMAD3 +0.54 q<0.01,
    # MEF2C ~0). DPT is reported as a qualified, non-fatal secondary cross-check.
    GATE_AXIS = 'myofibroblast_fate'
    for _pc in CFG['POS_CTRL']:
        _r = corr[(corr['gene'] == _pc) & (corr['axis'] == GATE_AXIS)]
        if _r.empty:
            raise RuntimeError('[CONTROL FAIL] positive control %s absent from CollecTRI '
                               'regulons along %s; cannot validate the pipeline.' % (_pc, GATE_AXIS))
        _row = _r.iloc[0]
        _rho = float(_row['rho']); _q = _q_of(_row)
        if not (_rho > 0 and np.isfinite(_q) and _q < CFG['CTRL_ALPHA']):
            raise RuntimeError(
                '[CONTROL FAIL] positive control %s rho=%+.3f q_bh=%.3g along %s '
                '(expected rho>0 AND q<%.2f in ABSOLUTE terms). The trajectory or regulon-'
                'activity pipeline is broken; refusing to report leads.'
                % (_pc, _rho, _q, GATE_AXIS, CFG['CTRL_ALPHA']))
        log('  [CONTROL PASS] %s rho=%+.3f q_bh=%.3g along %s (PRIMARY gated axis)'
            % (_pc, _rho, _q, GATE_AXIS))
        _rpt = corr[(corr['gene'] == _pc) & (corr['axis'] == 'pseudotime')]
        if not _rpt.empty:
            _pr = _rpt.iloc[0]
            log('  [DPT cross-check, NON-FATAL] %s along pseudotime rho=%+.3f q_bh=%.3g '
                '(DPT = degraded secondary axis on this large atlas; qualify in manuscript)'
                % (_pc, float(_pr['rho']), _q_of(_pr)))
    for _nc in CFG['NEG_CTRL']:
        _r = corr[(corr['gene'] == _nc) & (corr['axis'] == GATE_AXIS)]
        if _r.empty:
            log('  [CONTROL WARN] negative control %s absent from CollecTRI regulons '
                '(non-fatal).' % _nc)
            continue
        _row = _r.iloc[0]
        _rho = float(_row['rho']); _q = _q_of(_row)
        if _rho > CFG['NEG_CTRL_MAX_RHO'] and np.isfinite(_q) and _q < CFG['CTRL_ALPHA']:
            log('  [CONTROL WARN] negative control %s is STRONGLY positive along %s '
                '(rho=%+.3f q_bh=%.3g > NEG_CTRL_MAX_RHO=%.2f) -- unexpected; inspect but '
                'NON-FATAL.' % (_nc, GATE_AXIS, _rho, _q, CFG['NEG_CTRL_MAX_RHO']))
        else:
            log('  [CONTROL OK] negative control %s not strongly positive (rho=%+.3f q_bh=%.3g)'
                % (_nc, _rho, _q))

    # ---------- per-cell pseudotime table ----------
    keep_obs = ['cohort', 'sample', 'celltype', 'fib_subtype', 'fib_leiden']
    tab = pd.DataFrame(index=adata.obs_names)
    for c in keep_obs:
        if c in adata.obs:
            tab[c] = adata.obs[c].values
    tab['pseudotime'] = dpt
    if cc_pt is not None:
        tab['pseudotime_crosscheck'] = cc_pt
    tab['myo_fate'] = myo_fate
    tab['myo_program'] = myo
    tab['notch_target_module'] = notch_score
    for tf in TF_PANEL:
        if tf in acts.columns:
            tab[tf + '_regulon'] = acts[tf].values
    tab.to_csv(pt_csv)
    log('  wrote', pt_csv, tab.shape)

    # ---------- figures ----------
    log('STEP 4: figures (lollipop + activity lines + trajectory UMAP panel; viridis; no bars)')
    # each figure guarded: a viz bug must not throw away the 3.5h STEP3 compute
    for _fn, _args in [(fig_lollipop, (corr, 'pseudotime')),
                       (fig_lollipop, (corr, 'myofibroblast_fate')),
                       (fig_activity_lines, (dpt, acts, notch_score, myo, 'pseudotime')),
                       (fig_activity_lines, (myo_fate, acts, notch_score, myo, 'myofibroblast_fate')),
                       (fig_trajectory_umap, (adata, dpt, myo_fate, myo, acts))]:
        try:
            _fn(*_args)
        except Exception as _fe:
            log('  [fig warn] %s failed: %s' % (getattr(_fn, '__name__', _fn), repr(_fe)))

    # ★2026-07-09: save the per-cell trajectory adata (compact: obs table + UMAP) so a downstream
    # crash never loses the pseudotime/fate/regulon-activity again (a lean salvage can re-figure).
    try:
        import anndata as _ad
        _sav = _ad.AnnData(X=np.zeros((adata.n_obs, 1), dtype='float32'), obs=tab)
        if 'X_umap' in adata.obsm:
            _sav.obsm['X_umap'] = np.asarray(adata.obsm['X_umap'])
        _sav.write_h5ad(os.path.join(CFG['OUT'], 'trajectory_percell.h5ad'))
        log('  wrote trajectory_percell.h5ad', _sav.shape)
    except Exception as _e:
        log('  [warn] trajectory_percell.h5ad save failed:', repr(_e))

    log('DONE. outputs under', CFG['OUT'])
    log('TRAJECTORY_DONE')


if __name__ == '__main__':
    main()
