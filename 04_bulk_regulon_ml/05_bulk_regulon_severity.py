# -*- coding: utf-8 -*-
"""
============================================================================
05_bulk_regulon_severity.py   (SSc virtual-perturbation paper, POWERED build)
============================================================================
LAYER (this script): BULK regulon-SEVERITY across the full panel + DECONFOUNDED
immune coupling. It does NOT touch the single-cell virtual-KO / spatial layers.

MAIN LINE it feeds (context, not re-derived here):
  Single-cell virtual-perturbation engine (CellOracle in-silico TF-KO +
  decoupler CollecTRI regulon activity) nominates HES1/Notch (druggable ->
  nirogacestat) + FOSB/AP-1 (robust co-regulator + severity biomarker),
  converging on the SSc myofibroblast fibrotic program.

WHAT THIS SCRIPT DOES (3 stages, run/debug stage-by-stage via --stage):
  STAGE A  bulk regulon activity (decoupler CollecTRI ULM) per sample on the
           full bulk panel; HES1 / FOSB / SMAD3(+ctrl) / MEF2C(-ctrl) + a
           Notch-effector activity set vs SSc-status and vs mRSS.
           - per-cohort SSc-vs-HC (Welch t) + per-cohort mRSS (Spearman)
           - platform-as-covariate META (random-effects, DL) across cohorts
           - MIXED model for longitudinal cohorts (GSE76885 MMF trajectory,
             GSE249550 faSScinate TCZ): subject random intercept, timepoint fixed
  STAGE B  DECONFOUNDED immune coupling.
           - IOBR 7-method deconvolution (module 492) -> immune/stromal fractions
             + ESTIMATE purity/immune/stromal scores  (called out to R)
           - PARTIAL Spearman of immune fractions vs HES1/FOSB regulon activity,
             ADJUSTING FOR  ESTIMATE purity + total fibroblast load + total
             immune load + IFN score.  Anchored to mRSS where available.
           - reports HONESTLY whether the raw coupling COLLAPSES after adjustment.
  STAGE C  (optional, --ml) external-cohort ML: reuse the existing 059/Mime
           100+ combo consensus driver (ml_validation/ml_consensus_driver.R).
           This script only re-exports the regulon-activity feature tables the
           driver needs, then shells out to Rscript. Off by default.

HONEST-FRAMING (obey; baked into comments and outputs):
  * regulon activity = TF-target expression footprint. HIGH activity is a
    CORRELATE / HYPOTHESIS, NOT causal. Virtual-KO (upstream engine) is
    in-silico perturbation, not experimental knockout.
  * immune coupling = ASSOCIATION / co-context ONLY. We NEVER claim
    immune -> HES1. Direction is not identifiable from bulk correlation.
  * The raw FOSB/HES1 <-> immune correlations are LARGE for *every* population
    incl. fibroblasts (see ml_validation/mcp_vs_regulon_corr.csv: rho 0.5-0.77).
    That pattern screams a shared cellularity / purity confound -> this whole
    stage exists to test if any coupling SURVIVES purity+load+IFN adjustment.
  * TREM2 macrophages are PROTECTIVE in SSc skin (not pro-fibrotic); do not
    frame any myeloid signal as pro-fibrotic.
  * nirogacestat has NO immune endpoints; anything immune here is hypothesis
    generation for the fibroblast program, not a drug-immune claim.

REPRODUCIBILITY (obey): fixed SEED; every stage checkpoints to disk and is
idempotent (skips if output exists unless --force); ALL paths come from the
CONFIG block below (edit CFG for the server, no hardcoding downstream); logs
progress; memory-friendly (per-cohort streaming, gene x sample frames, no
giant concat).

ENV: run on server env `npc` (scanpy / anndata / decoupler / statsmodels /
pingouin / scipy). STAGE B's deconvolution shells out to `Rscript` with IOBR
installed (module 492). IOBR + ESTIMATE are MANDATORY: if Rscript/IOBR is
unavailable (or ESTIMATE fails), STAGE B HARD-FAILS rather than substituting a
marker-mean surrogate that would fabricate the very purity confounder it adjusts
for (rules #1/#3).

REUSES (read, adapted, real APIs only):
  replication_bulk/compute_bulk_replication.py   (ULM, GPL map, cohort parsing)
  ml_validation/mcp_immune.py                    (MCP-counter + regulon coupling)
  ml_validation/build_gse*_features.py           (decoupler ULM feature export)
  ml_validation/ml_consensus_driver.R            (059/Mime ML, STAGE C)
  library module 492 IOBR deconvolution ; module 076 decoupler
============================================================================
"""
import os, sys, io, gzip, re, json, argparse, subprocess, warnings
from io import StringIO

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# Do NOT blanket-ignore: keep statistical/convergence warnings visible (rules #1/#6).
# Silence only known-benign, specifically-scoped categories if needed.
warnings.simplefilter('default')
# e.g. warnings.filterwarnings('ignore', message='.*use_inf_as_na.*')  # scope narrowly

import numpy as np
import pandas as pd

SEED = 123
np.random.seed(SEED)

# =====================================================================
# CONFIG  (edit ONLY this block on the server; everything else is path-free)
# =====================================================================
CFG = dict(
    # ---- roots ----
    RAW      = os.environ.get('SSC_RAW',  '/data/ssc/raw'),        # downloaded GEO raw lives here
    OUT      = os.environ.get('SSC_OUT',  './out_05_bulk'),        # this script's outputs
    NET_CSV  = os.environ.get('SSC_NET',  '/data/ssc/data/collectri_net.csv'),  # SAME net as discovery
    RSCRIPT  = os.environ.get('RSCRIPT',  'Rscript'),              # for IOBR (module 492) + 059 ML

    # ---- lead TFs + control TFs (engine-nominated; SMAD3=pos ctrl, MEF2C=neg ctrl) ----
    LEADS      = ['HES1', 'FOSB', 'SMAD3', 'MEF2C', 'CEBPD'],
    # Notch-effector activity set: HES1's canonical Notch/HES-HEY module TFs.
    # We report the MEAN regulon activity over those present (a "Notch-effector score").
    NOTCH_SET  = ['HES1', 'HES5', 'HEY1', 'HEY2', 'HEYL', 'RBPJ'],
    # IFN score genes (type-I ISG core) for the deconfounding covariate (STAGE B).
    IFN_GENES  = ['ISG15', 'IFI6', 'IFI44', 'IFI44L', 'IFIT1', 'IFIT3', 'MX1',
                  'MX2', 'OAS1', 'OAS2', 'OAS3', 'RSAD2', 'USP18', 'SIGLEC1'],
    # Fibroblast marker set (for "total fibroblast load" confounder in STAGE B).
    FIBRO_GENES= ['COL1A1', 'COL1A2', 'COL3A1', 'PDGFRA', 'PDGFRB', 'LUM',
                  'DCN', 'FBN1', 'THY1', 'PDPN'],
    # Immune marker set (for "total immune load" confounder; PTPRC=CD45 anchor).
    IMMUNE_GENES=['PTPRC', 'CD3D', 'CD3E', 'CD68', 'CD14', 'CD19', 'MS4A1',
                  'NKG7', 'LYZ', 'ITGAM'],

    SEED = SEED,
)

# ---------------------------------------------------------------------
# COHORT REGISTRY. Each entry declares HOW to load one bulk cohort from
# CFG['RAW']. `loader` picks the parse path proven in compute_bulk_replication.py.
# On the server, point `file` at the actual downloaded name; leave `mrss_col`/
# `long` where the cohort truly has severity / longitudinal structure.
#
# NOTE (honest data provenance): GSE76885 & GSE9285 are 2-colour arrays whose
# SSc-vs-HC positive control (SMAD3) FAILED in prior local replication
# (see replication_bulk/BULK_REPLICATION_SUMMARY.md). They are retained here
# for TRANSPARENCY and screened by the SMAD3 positive-control gate at meta time;
# they are NOT silently dropped and NOT cherry-picked in.
# ---------------------------------------------------------------------
COHORTS = [
    dict(gse='GSE76885',  platform='Agilent-2c', loader='series_matrix',
         file='GSE76885_series_matrix.txt.gz', gpl='GPL6480.txt',  sym='GENE_SYMBOL',
         long=True,  mrss_col=None,  note='MMF trajectory; baseline-restricted; pos-ctrl watched'),
    dict(gse='GSE58095',  platform='Illumina-BA', loader='geoparse',
         file='GSE58095_family.soft.gz',        gpl=None,          sym=None,  grp_mode='control_vs_rest',
         long=False, mrss_col='mrss', note='largest single-platform SSc skin array (merged superseries -> label via control-vs-rest)'),
    dict(gse='GSE130955', platform='RNAseq-FPKM', loader='fpkm_raw_tar',
         file='GSE130955_RAW.tar',               gpl=None,          sym=None,
         long=False, mrss_col='mrss', note='DISCOVERY bulk; pipeline anchor'),
    dict(gse='GSE106358', platform='RNAseq',      loader='geoparse',
         file='GSE106358_family.soft.gz',        gpl=None,          sym=None,
         long=False, mrss_col='mrss', note='SSc skin RNA-seq'),
    dict(gse='GSE9285',   platform='Agilent-2c',  loader='geoparse',
         file='GSE9285_family.soft.gz',          gpl='GPL5981.txt', sym='GENE',
         long=False, mrss_col=None,  note='Milano 2008; tiny HC; pos-ctrl watched'),
    dict(gse='GSE249550', platform='Agilent-2c',  loader='geoparse',
         file='GSE249550_family.soft.gz',        gpl='GPL24205.txt',sym='Symbol',
         long=True,  mrss_col='mrss', note='faSScinate TCZ; subject/timepoint/mrss'),
]

# =====================================================================
# small utilities
# =====================================================================
def log(*a):
    print('[05]', *a, flush=True)

def ensure_out():
    for sub in ['', 'stageA', 'stageB', 'stageC', 'figs', 'features']:
        os.makedirs(os.path.join(CFG['OUT'], sub), exist_ok=True)

def done(path, force):
    """idempotency guard: True if we can skip (exists and not forced)."""
    if os.path.exists(path) and not force:
        log('SKIP (exists):', path)
        return True
    return False

# --- shared CollecTRI net (IDENTICAL file as discovery; see BULK_REPLICATION_SUMMARY) ---
_NET = None
def net():
    global _NET
    if _NET is None:
        n = pd.read_csv(CFG['NET_CSV'])
        # tolerate either raw collectri dump or a pre-sliced 3-col file
        cols = {c.lower(): c for c in n.columns}
        ren = {}
        for a, b in [('tf', 'source'), ('target_genesymbol', 'target'), ('mor', 'weight')]:
            if a in cols and b not in cols:
                ren[cols[a]] = b
        n = n.rename(columns=ren)
        if 'weight' not in n.columns:
            n['weight'] = 1.0
        _NET = n[['source', 'target', 'weight']].dropna()
        log('CollecTRI net loaded:', _NET.shape, 'edges,', _NET['source'].nunique(), 'TFs')
    return _NET

# =====================================================================
# ADAPTED from compute_bulk_replication.py  (GPL map, probe collapse, ULM)
# =====================================================================
def gpl_map(gpl_txt, symcol):
    import GEOparse
    g = GEOparse.get_GEO(filepath=os.path.join(CFG['RAW'], gpl_txt), silent=True)
    t = g.table.copy()
    t['ID'] = t['ID'].astype(str)
    m = t[['ID', symcol]].dropna()
    m = m[~m[symcol].astype(str).isin(['', 'nan', 'NA'])]
    return dict(zip(m['ID'], m[symcol].astype(str)))

def collapse_to_genes(mat, probe2sym):
    """probe x sample -> gene x sample (mean over probes). Same as replication."""
    mat = mat.copy()
    mat.index = mat.index.astype(str).map(lambda p: probe2sym.get(p))
    mat = mat[~mat.index.isna()]
    mat = mat.apply(pd.to_numeric, errors='coerce')
    return mat.groupby(mat.index).mean()

def run_ulm(expr_gene_by_sample, platform=None):
    """gene x sample -> ULM activity sample x TF (decoupler 2.x, module 076 API).
    NaN/inf cleaning identical to replication (2-colour arrays have missing probes).
    A6: for RNA-seq (count/FPKM) cohorts, drop near-zero-expression genes with
    decoupler's edgeR-style filter_by_expr BEFORE ULM. Arrays are already probe-
    filtered upstream, so the filter is applied ONLY to RNA-seq platforms."""
    import anndata as ad
    import decoupler as dc
    e = expr_gene_by_sample.replace([np.inf, -np.inf], np.nan).dropna(how='all')
    e = e.T.fillna(e.mean(axis=1)).T.fillna(0.0)
    A = ad.AnnData(e.T.values.astype('float32'))   # samples x genes (filter_by_expr input)
    A.obs_names = list(e.columns)
    A.var_names = list(e.index)
    # A6: low-expression gene filter for RNA-seq cohorts only (arrays -> skip, probe-filtered).
    # ★2026-07-07 修: 这些 RNA-seq 队列到此已是 log 空间(loader 做过 log1p/log2), 不能用 edgeR
    # 生计数过滤器 dc.pp.filter_by_expr(min_count=10)——log 值 0-8 全被滤→0 基因→dc.mt.ulm 崩,
    # 锚定队列 GSE130955 挂。改用 log 适配的"表达占比"过滤(基因在 ≥10% 样本中 >0), 且只在留够
    # 基因(≥50)时才应用(防清空)。
    if platform and str(platform).lower().startswith('rnaseq'):
        n0 = A.n_vars
        keep = np.asarray((A.X > 0).mean(axis=0)).ravel() >= 0.1
        if int(keep.sum()) >= 50:
            A = A[:, keep].copy()
        print(f'  [A6] RNA-seq expr-fraction filter ({platform}): {n0} -> {A.n_vars} genes')
    n = net()
    n = n[n['target'].isin(A.var_names)]
    np.random.seed(CFG['SEED'])  # explicit at stochastic call site (rule #6)
    dc.mt.ulm(data=A, net=n, tmin=5, verbose=False)
    # pick the score obsm (not pval/padj) — same selection as replication
    k = [x for x in A.obsm if 'ulm' in x.lower()
         and 'pval' not in x.lower() and 'padj' not in x.lower()][0]
    return A.obsm[k]                      # DataFrame: samples x TFs

# =====================================================================
# COHORT LOADERS  (return: expr gene x sample, meta[grp/mrss/subject/time])
# each mirrors a proven branch of compute_bulk_replication.py
# =====================================================================
def _parse_series_matrix(fp):
    with gzip.open(fp, 'rt', encoding='utf-8', errors='replace') as f:
        L = f.readlines()
    beg = next(i for i, l in enumerate(L) if l.startswith('!series_matrix_table_begin'))
    end = next(i for i, l in enumerate(L) if l.startswith('!series_matrix_table_end'))
    tbl = pd.read_csv(StringIO(''.join(L[beg + 1:end])), sep='\t', low_memory=False)
    tbl = tbl.rename(columns={tbl.columns[0]: 'ID'}).set_index('ID')
    meta = {}
    for l in L[:beg]:
        if l.startswith('!Sample_'):
            p = l.rstrip('\n').split('\t')
            meta.setdefault(p[0], []).append([x.strip().strip('"') for x in p[1:]])
    return tbl, meta

def load_GSE76885(co):
    """series-matrix parse; baseline-restricted SSc; subject-level EXPRESSION mean,
    keep per-timepoint too for the mixed model. (adapted from replication)."""
    tbl, meta = _parse_series_matrix(os.path.join(CFG['RAW'], co['file']))
    acc    = meta['!Sample_geo_accession'][0]
    titles = meta['!Sample_title'][0]
    dis = subj = None
    for row in meta.get('!Sample_characteristics_ch1', []):
        if row and 'disease state' in row[0].lower():
            dis = [v.split(':', 1)[1].strip() for v in row]
        if row and row[0].lower().startswith('subject'):
            subj = [v.split(':', 1)[1].strip() for v in row]
    grp, tp = [], []
    for d, t in zip(dis, titles):
        dl, tl = d.lower(), t.lower()
        if 'normal' in dl or 'control' in dl:
            grp.append('HC'); tp.append('base')
        elif 'ssc' in dl:
            grp.append('SSc')
            m = re.search(r'(base|0?6\s*mo|12\s*mo|24\s*mo|36\s*mo)', tl)
            tp.append(m.group(1).replace(' ', '') if m else 'base')
        else:
            grp.append('EXC'); tp.append('na')
    meta_df = pd.DataFrame({'grp': grp, 'subject': subj, 'time': tp}, index=acc)
    expr = collapse_to_genes(tbl, gpl_map(co['gpl'], co['sym']))
    return expr, meta_df

def load_geoparse(co):
    """generic cached-soft loader (GSE58095 / GSE106358 / GSE9285 / GSE249550)."""
    import GEOparse
    g = GEOparse.get_GEO(filepath=os.path.join(CFG['RAW'], co['file']), silent=True)
    mat = g.pivot_samples('VALUE')
    pdt = g.phenotype_data
    # log2 for high-range linear arrays / counts (same guard as build_gse58095)
    if np.nanmax(mat.values.astype(float)) > 50:
        mat = np.log2(mat.clip(lower=1))
    # ---- group (★2026-07-08: scan ALL group/disease-like columns + consolidate; robust to a
    #      MERGED superseries like GSE58095 whose label is spread over batch-specific fields, so the
    #      old single 'disease'-column lookup returned all-EXC. prefer disease/diagnosis cols first) ----
    _gcols = [c for c in pdt.columns
              if any(k in c.lower() for k in ['group', 'disease', 'diagnos', 'condition', 'phenotype'])
              and 'antibody' not in c.lower()]          # 'antibody group' is NOT a disease label
    _gcols.sort(key=lambda c: 0 if any(k in c.lower() for k in ['disease', 'diagnos', 'condition']) else 1)
    if not _gcols and 'title' in pdt.columns:
        _gcols = ['title']
    def _lab_row(i):
        for c in _gcols:
            x = str(pdt.at[i, c]).lower()
            if any(k in x for k in ['healthy', 'control', 'normal']):
                return 'HC'
            if any(k in x for k in ['ssc', 'sclerosis', 'scleroderma']):   # 'ssc' also catches dcSSc/lcSSc
                return 'SSc'
        return 'EXC'
    grp = pd.Series([_lab_row(i) for i in pdt.index], index=pdt.index)
    # control-vs-rest studies (pure SSc-vs-control skin, e.g. GSE58095 merged superseries where most
    # SSc samples carry only batch-specific clinical fields, not an explicit 'SSc' group): any non-
    # control skin biopsy = SSc. Gated by co['grp_mode'] so it is NEVER applied blindly to other cohorts.
    if co.get('grp_mode') == 'control_vs_rest':
        _src = (pdt['source_name_ch1'].astype(str).str.lower()
                if 'source_name_ch1' in pdt.columns else pd.Series(['skin'] * len(pdt), index=pdt.index))
        grp[(grp == 'EXC') & _src.str.contains('skin|biopsy|dermis', na=False)] = 'SSc'
    # ---- mrss / subject / time (best-effort column discovery) ----
    def find(colkeys):
        for c in pdt.columns:
            cl = c.lower()
            if any(k in cl for k in colkeys):
                return pd.to_numeric(pdt[c], errors='coerce') if 'mrss' in colkeys \
                       else pdt[c].astype(str)
        return pd.Series(index=pdt.index, dtype=float if 'mrss' in colkeys else object)
    meta_df = pd.DataFrame(index=pdt.index)
    meta_df['grp']     = grp.values
    meta_df['mrss']    = find(['mrss', 'rodnan', 'skin score']).values
    meta_df['subject'] = find(['subject', 'patient', 'individual']).values
    meta_df['time']    = find(['time point', 'timepoint', 'visit', 'week']).values
    # ---- probe -> symbol ----
    if co['gpl']:
        expr = collapse_to_genes(mat, gpl_map(co['gpl'], co['sym']))
    else:
        gpl = list(g.gpls.values())[0]; ann = gpl.table
        # ★2026-07-08 修: 按偏好顺序选符号列(HGNC Symbol 优先于 Illumina 专有 ILMN_Gene) —— 原来按
        # 列出现顺序取 symcol[0] 会选到 ILMN_Gene, 其名不匹配 CollecTRI 的 HGNC 靶 -> run_ulm
        # "no sources with >5 targets" 崩。
        _avail = {c.lower(): c for c in ann.columns}
        symcol = next((_avail[p] for p in ('symbol', 'gene_symbol', 'genesymbol', 'gene', 'ilmn_gene')
                       if p in _avail), None)
        if symcol is None:
            raise ValueError(f"{co['gse']}: 嵌入 GPL 无可用符号列 (cols={list(ann.columns)[:10]})")
        idcol = 'ID' if 'ID' in ann.columns else ann.columns[0]
        p2s = ann.set_index(idcol)[symcol]
        mat.index = mat.index.map(p2s)
        mat = mat[mat.index.notna() & (mat.index.astype(str) != '')]
        expr = mat.groupby(level=0).max()
    return expr, meta_df

def load_fpkm_raw_tar(co):
    """GSE130955 discovery bulk from per-sample RAW.tar (adapted from
    build_gse130955_features.py). Labels from _P / _N suffix."""
    import tarfile, mygene
    tar = tarfile.open(os.path.join(CFG['RAW'], co['file'])); ser = {}
    for m in tar.getmembers():
        if not m.name.endswith('.txt.gz'):
            continue
        sid = m.name.split('_', 1)[1][:-7]
        df = pd.read_csv(gzip.open(tar.extractfile(m)), sep='\t')
        df['ens'] = df.iloc[:, 0].astype(str).str.split('.').str[0]
        ser[sid] = df.groupby('ens')[df.columns[1]].sum()
    mat = pd.DataFrame(ser)
    mg = mygene.MyGeneInfo()
    q = mg.querymany(mat.index.tolist(), scopes='ensembl.gene', fields='symbol',
                     species='human', as_dataframe=True, verbose=False)
    sym = q['symbol'].dropna(); sym = sym[~sym.index.duplicated()]
    mat = mat.loc[mat.index.isin(sym.index)]; mat.index = sym.loc[mat.index].values
    expr = np.log1p(mat.groupby(level=0).max())      # log1p FPKM (build_gse130955)
    grp = ['SSc' if s.endswith('_P') else ('HC' if s.endswith('_N') else 'EXC')
           for s in expr.columns]
    meta_df = pd.DataFrame({'grp': grp}, index=expr.columns)
    meta_df['mrss'] = np.nan; meta_df['subject'] = expr.columns; meta_df['time'] = 'na'
    return expr, meta_df

LOADERS = {'series_matrix': load_GSE76885, 'geoparse': load_geoparse,
           'fpkm_raw_tar': load_fpkm_raw_tar}

def load_cohort(co):
    fn = LOADERS[co['loader']]
    if co['loader'] == 'series_matrix':          # GSE76885 special-cased
        return load_GSE76885(co)
    return fn(co)

# =====================================================================
# score helpers
# =====================================================================
def zscore_mean(expr_gene_by_sample, genes):
    """per-gene z-score across samples then mean over present genes -> sample score.
    Used for IFN / fibroblast-load / immune-load / Notch confounders and signatures."""
    g = [x for x in genes if x in expr_gene_by_sample.index]
    if not g:
        return pd.Series(np.nan, index=expr_gene_by_sample.columns)
    sub = expr_gene_by_sample.loc[g].apply(pd.to_numeric, errors='coerce')
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)

# =====================================================================
# STAGE A : bulk regulon activity vs SSc-status and vs mRSS (per-cohort + meta + mixed)
# =====================================================================
def stage_A(force):
    from scipy import stats
    out_csv = os.path.join(CFG['OUT'], 'stageA', 'bulk_regulon_severity.csv')
    act_dir = os.path.join(CFG['OUT'], 'features')
    if done(out_csv, force) and not force:
        return pd.read_csv(out_csv)

    rows = []
    failed_cohorts = []
    for co in COHORTS:
        gse = co['gse']
        fpath = os.path.join(CFG['RAW'], co['file'])
        if not os.path.exists(fpath):
            log('MISSING raw, skipping cohort:', gse, '->', fpath); continue
        log(f'== STAGE A :: {gse} ({co["platform"]}) :: {co["note"]} ==')

        # ---- checkpoint per-cohort regulon-activity table (also feeds STAGE C) ----
        act_path = os.path.join(act_dir, f'{gse}_regulon_activity.csv')
        meta_path = os.path.join(act_dir, f'{gse}_meta.csv')
        if os.path.exists(act_path) and os.path.exists(meta_path) and not force:
            act = pd.read_csv(act_path, index_col=0)
            meta = pd.read_csv(meta_path, index_col=0)
            log('  reuse cached regulon activity:', act.shape)
        else:
            try:
                expr, meta = load_cohort(co)
            except Exception as e:
                log(f'  [SKIP cohort] {gse} load_cohort 失败: {e}'); failed_cohorts.append(gse); continue
            log('  gene x sample:', expr.shape,
                '| grp:', meta['grp'].value_counts().to_dict())
            keep = meta.index[meta['grp'].isin(['SSc', 'HC'])]
            # ★2026-07-08: SSc/HC 分类不足(如 grp 全 EXC)或基因-net 重叠不足 -> 跳过该 cohort(响亮),
            # 不 abort 整个 StageA; SMAD3 阳性对照在幸存 cohort 上仍验证(非静默降级: 失败已计入 failed_cohorts)。
            if len(keep) < 6:
                log(f'  [SKIP cohort] {gse}: 仅 {len(keep)} 个 SSc/HC 样本(grp 分类失败?) -> 跳过')
                failed_cohorts.append(gse); continue
            expr_k, meta = expr[keep], meta.loc[keep].copy()
            # per-sample regulon activity (keep raw samples; aggregation handled per-analysis)
            try:
                act = run_ulm(expr_k, platform=co['platform'])   # A6: platform gates the RNA-seq low-count filter
            except Exception as e:
                log(f'  [SKIP cohort] {gse}: run_ulm 失败({e}) -> 基因与 CollecTRI 靶重叠不足, 跳过')
                failed_cohorts.append(gse); continue
            act.to_csv(act_path); meta.to_csv(meta_path)
            # also stash gene scores meta needs for STAGE B / mixed models
            meta['IFN_score']   = zscore_mean(expr_k, CFG['IFN_GENES']).reindex(meta.index)
            meta['fibro_load']  = zscore_mean(expr_k, CFG['FIBRO_GENES']).reindex(meta.index)
            meta['immune_load'] = zscore_mean(expr_k, CFG['IMMUNE_GENES']).reindex(meta.index)
            meta.to_csv(meta_path)
            log('  regulon activity (samples x TFs):', act.shape)

        # ---- Notch-effector score = mean regulon activity over present Notch TFs ----
        notch_present = [t for t in CFG['NOTCH_SET'] if t in act.columns]
        act = act.copy()
        act['NOTCH_EFFECTOR'] = act[notch_present].mean(axis=1) if notch_present else np.nan

        meta = meta.loc[act.index]
        is_ssc = (meta['grp'] == 'SSc').values
        is_hc  = (meta['grp'] == 'HC').values

        # ---- (1) SSc vs HC per TF/score (Welch t) ----
        for tf in CFG['LEADS'] + ['NOTCH_EFFECTOR']:
            if tf not in act.columns:
                continue
            x = act[tf].values[is_ssc]; y = act[tf].values[is_hc]
            if len(x) < 3 or len(y) < 2:
                continue
            t, p = stats.ttest_ind(x, y, equal_var=False, nan_policy='omit')
            rows.append(dict(cohort=gse, platform=co['platform'], analysis='SSc_vs_HC',
                             tf=tf, effect=float(np.nanmean(x) - np.nanmean(y)),
                             t_stat=float(t), p=float(p),
                             n_ssc=int(is_ssc.sum()), n_hc=int(is_hc.sum())))

        # ---- (2) mRSS Spearman (SSc only) where cohort truly has mRSS ----
        if co['mrss_col'] and 'mrss' in meta.columns:
            mr = pd.to_numeric(meta['mrss'], errors='coerce')
            ok = is_ssc & mr.notna().values
            if ok.sum() >= 6:
                for tf in ['FOSB', 'HES1', 'NOTCH_EFFECTOR']:
                    if tf in act.columns:
                        rho, p = stats.spearmanr(act[tf].values[ok], mr.values[ok])
                        rows.append(dict(cohort=gse, platform=co['platform'], analysis='mRSS',
                                         tf=tf, effect=float(rho), p=float(p),
                                         n_ssc=int(ok.sum()), n_hc=0))

        # ---- (3) MIXED model for longitudinal cohorts (subject random intercept) ----
        # HONEST: for trajectory cohorts, timepoints within a subject are NOT
        # independent; a subject-random-intercept model gives the correct SE for
        # the SSc effect and lets us keep on-treatment timepoints as a fixed effect
        # instead of discarding them.
        if co['long']:
            _mixed_ssc_effect(act, meta, gse, co, rows)

    R = pd.DataFrame(rows)
    # ---- BH-FDR within each analysis family (per-TF x per-cohort tests) ----
    if not R.empty and 'p' in R.columns:
        from statsmodels.stats.multitest import multipletests
        R['padj'] = np.nan
        for _an, _idx in R.groupby('analysis').groups.items():
            _p = R.loc[_idx, 'p']
            _ok = _p.notna()
            if _ok.sum() > 0:
                R.loc[_p.index[_ok], 'padj'] = multipletests(
                    _p[_ok].values, method='fdr_bh')[1]
    R.to_csv(out_csv, index=False)
    log('STAGE A rows:', len(R), '-> (BH-FDR per analysis family) ->', out_csv)

    # ---- META across cohorts (random-effects DerSimonian-Laird), platform noted ----
    _meta_analyze(R)
    # ---- figs (dumbbell + dot; NO bars) ----
    _fig_stageA(R)
    return R

def _mixed_ssc_effect(act, meta, gse, co, rows):
    """subject-random-intercept mixed model: activity ~ SSc + time  (+ 1|subject).
    statsmodels MixedLM. Reports the SSc fixed-effect coefficient + p per TF."""
    import statsmodels.formula.api as smf  # required dep; missing -> hard-fail (rule #1)
    m = meta.copy()
    m['SSc'] = (m['grp'] == 'SSc').astype(int)
    m['subject'] = m['subject'].astype(str).replace('nan', np.nan)
    m['time'] = m['time'].astype(str).fillna('na')
    if m['subject'].isna().all() or m['subject'].nunique() < 3:
        log('  [mixed skipped: no usable subject id]'); return
    for tf in CFG['LEADS'] + ['NOTCH_EFFECTOR']:
        if tf not in act.columns:
            continue
        d = m.copy(); d['y'] = act[tf].reindex(d.index).values
        d = d.dropna(subset=['y', 'subject'])
        if d['SSc'].nunique() < 2 or d['subject'].nunique() < 3:
            continue
        formula = 'y ~ SSc + C(time)' if d['time'].nunique() > 1 else 'y ~ SSc'
        fit = smf.mixedlm(formula, d, groups=d['subject']).fit(reml=False, disp=False)
        if not fit.converged:
            raise RuntimeError(f'MixedLM did NOT converge for {gse}/{tf} '
                               f'(formula={formula}); fix the model/data rather '
                               f'than silently dropping this estimate.')
        if 'SSc' not in fit.params.index:
            raise RuntimeError(f'MixedLM for {gse}/{tf} produced no SSc term '
                               f'(design collapsed); investigate, do not skip.')
        rows.append(dict(cohort=gse, platform=co['platform'], analysis='mixed_SSc',
                         tf=tf, effect=float(fit.params['SSc']),
                         p=float(fit.pvalues['SSc']),
                         n_ssc=int((d['SSc'] == 1).sum()),
                         n_hc=int((d['SSc'] == 0).sum())))

def _meta_analyze(R):
    """Random-effects meta (DerSimonian-Laird) of SSc-vs-HC effect per TF across
    cohorts, using per-cohort SE from t and n. Platform is reported alongside as a
    moderator note (not modelled — too few cohorts for meta-regression).
    POSITIVE-CONTROL GATE: cohorts where SMAD3 is not UP are flagged, not dropped."""
    sub = R[R.analysis == 'SSc_vs_HC'].copy()
    if sub.empty:
        return
    # approximate SE of a mean-difference in z-units from t: SE = effect / t; recover
    # t from p is fragile, so approximate SE via effect and n (pooled) -> use a simple
    # inverse-variance with se ~ |effect|/max(z,eps) where z from p (two-sided).
    from scipy.stats import norm
    # Cochrane inverse-variance meta needs the SE of the mean difference from the
    # PRIMARY test, not a back-derivation from p. Welch SE = |effect| / |t|, and
    # t is exact from the two-sample test. Carry t (and per-group SDs/n) out of
    # the Welch step (see F07b) so SE is not reconstructed from p.
    if 't_stat' not in sub.columns:
        raise RuntimeError('meta requires per-cohort Welch t (t_stat). Emit it in '
                           'stage_A SSc_vs_HC rows (see F07b); refusing to '
                           'back-derive SE from p-values (self-created approximation).')
    sub['se'] = (sub['effect'].abs() / sub['t_stat'].abs()).replace([np.inf, 0], np.nan)
    # SMAD3 positive-control status per cohort
    pc = sub[sub.tf == 'SMAD3'].set_index('cohort')['effect']
    sub['pc_pass'] = sub['cohort'].map(lambda c: bool(pc.get(c, -1) > 0))
    out = []
    for tf, g in sub.groupby('tf'):
        g = g.dropna(subset=['se'])
        # meta over ALL cohorts and over positive-control-passing cohorts only
        for tag, gg in [('all', g), ('pc_pass_only', g[g.pc_pass])]:
            if len(gg) == 0:
                continue
            if len(gg) < 2:   # ★2026-07-08 BUG-C: 单队列不是跨队列 meta (k=1 -> tau2=0 -> RE=该单队列), 会假通过 SMAD3 阳性对照门授权 HES1/FOSB 头条
                log(f'  [meta skipped: k={len(gg)} for {tf}/{tag}]')
                continue
            w = 1 / gg['se']**2
            fe = (w * gg['effect']).sum() / w.sum()
            Q = (w * (gg['effect'] - fe)**2).sum()
            df = len(gg) - 1
            tau2 = max(0, (Q - df) / (w.sum() - (w**2).sum() / w.sum())) if df > 0 else 0
            wr = 1 / (gg['se']**2 + tau2)
            re = (wr * gg['effect']).sum() / wr.sum()
            se_re = np.sqrt(1 / wr.sum())
            z = re / se_re
            p = 2 * norm.sf(abs(z))
            I2 = max(0, (Q - df) / Q) * 100 if Q > 0 else 0
            out.append(dict(tf=tf, subset=tag, k=len(gg), RE_effect=re, se=se_re,
                            p=p, I2=I2, cohorts=';'.join(gg['cohort'])))
    M = pd.DataFrame(out)
    mp = os.path.join(CFG['OUT'], 'stageA', 'bulk_regulon_meta.csv')
    M.to_csv(mp, index=False)
    log('META saved ->', mp)
    # ---- POSITIVE-CONTROL GATE (rule #5): SMAD3 must be UP and significant in
    #      the pc_pass_only meta, else the regulon-activity harness is not
    #      validated and downstream HES1/FOSB claims are unsupported. HARD-FAIL. ----
    smad = M[(M.tf == 'SMAD3') & (M.subset == 'pc_pass_only')]
    if smad.empty:
        raise RuntimeError('POSITIVE-CONTROL GATE FAILED: no SMAD3 meta estimate '
                           '(pc_pass_only). Cannot validate the regulon harness.')
    _eff, _p = float(smad['RE_effect'].iloc[0]), float(smad['p'].iloc[0])
    if not (_eff > 0 and _p < 0.05):
        raise RuntimeError(
            f'POSITIVE-CONTROL GATE FAILED: SMAD3 meta RE_effect={_eff:.3f} '
            f'p={_p:.3g} (must be >0 and p<0.05). The bulk regulon-activity '
            f'harness is not validated; HES1/FOSB nominations are unsupported.')
    # negative control: MEF2C must NOT show a strong significant UP signal
    mef = M[(M.tf == 'MEF2C') & (M.subset == 'pc_pass_only')]
    if not mef.empty and float(mef['RE_effect'].iloc[0]) > 0 and float(mef['p'].iloc[0]) < 0.05:
        log('  [WARN neg-ctrl] MEF2C reads significantly UP in meta - inspect harness specificity.')
    log(f'  POSITIVE-CONTROL PASS: SMAD3 meta RE_effect={_eff:.3f} p={_p:.3g}')
    if not M.empty:
        head = M[M.subset == 'pc_pass_only'].sort_values('p')
        log('  META (positive-control-passing cohorts) top:',
            head[['tf', 'k', 'RE_effect', 'p', 'I2']].round(3).to_dict('records')[:4])

# =====================================================================
# STAGE B : DECONFOUNDED immune coupling
# =====================================================================
def stage_B(force):
    from scipy.stats import spearmanr
    out_csv = os.path.join(CFG['OUT'], 'stageB', 'immune_coupling_deconfounded.csv')
    if done(out_csv, force):
        return pd.read_csv(out_csv)

    rows = []
    for co in COHORTS:
        gse = co['gse']
        act_path  = os.path.join(CFG['OUT'], 'features', f'{gse}_regulon_activity.csv')
        meta_path = os.path.join(CFG['OUT'], 'features', f'{gse}_meta.csv')
        if not (os.path.exists(act_path) and os.path.exists(meta_path)):
            # ★2026-07-08: a cohort with NO STAGE-A cache = it failed classification/QC in STAGE A
            # (e.g. GSE9285 unresolved SSc/HC labels) and was LOUDLY skipped there. Do NOT hard-fail
            # STAGE B on it -- LOUDLY skip here too (honest, not silent) and deconfound on the cohorts
            # that DO have a cache. (If NO cohort has a cache, run --stage A first.)
            log(f'  [STAGE B skip: {gse} has no STAGE-A cache (skipped/failed in STAGE A) '
                f'-> deconfounding proceeds WITHOUT it, LOGGED not silent]')
            continue
        act  = pd.read_csv(act_path, index_col=0)
        meta = pd.read_csv(meta_path, index_col=0)
        if not {'HES1', 'FOSB'}.issubset(act.columns):
            continue
        log(f'== STAGE B :: {gse} :: deconfounded immune coupling ==')

        # ---- (i) immune/stromal fractions + ESTIMATE purity ----
        frac = _deconvolve(co, gse, force)   # samples x {cell fractions, ESTIMATE_*}
        if frac is None or frac.empty:
            log('  [no deconvolution for', gse, '- skipping coupling]'); continue

        # align everything on common samples
        idx = act.index.intersection(frac.index).intersection(meta.index)
        act_c, frac_c, meta_c = act.loc[idx], frac.loc[idx], meta.loc[idx]
        # confounders: ESTIMATE purity + total fibroblast load + total immune load + IFN
        conf = pd.DataFrame(index=idx)
        pur_col = [c for c in frac_c.columns if 'purity' in c.lower() or c.lower() == 'tumorpurity']
        conf['purity']      = frac_c[pur_col[0]] if pur_col else np.nan
        conf['fibro_load']  = meta_c.get('fibro_load')
        conf['immune_load'] = meta_c.get('immune_load')
        conf['IFN_score']   = meta_c.get('IFN_score')
        # drop all-NaN confounders (some cohorts lack ESTIMATE purity)
        conf = conf.dropna(axis=1, how='all')
        log('  confounders used:', list(conf.columns))

        # cell-fraction columns only (exclude ESTIMATE scalar scores from the "immune pops")
        estimate_cols = [c for c in frac_c.columns
                         if any(k in c.lower() for k in
                                ['estimate', 'stromalscore', 'immunescore', 'purity'])]
        pop_cols = [c for c in frac_c.columns if c not in estimate_cols]

        for tf in ['HES1', 'FOSB']:
            for pop in pop_cols:
                x = pd.to_numeric(frac_c[pop], errors='coerce')
                y = pd.to_numeric(act_c[tf], errors='coerce')
                # raw Spearman
                m0 = x.notna() & y.notna()
                if m0.sum() < 8:
                    continue
                rho_raw, p_raw = spearmanr(x[m0], y[m0])
                # PARTIAL Spearman adjusting for confounders (residual-on-residual)
                rho_par, p_par, n_used = _partial_spearman(x, y, conf)
                # INTERIM survival flag (D1): the self-created |rho_partial|<0.15
                # magnitude cutoff is DROPPED. Survival = the partial association is
                # significant AND sign-consistent with the raw correlation. This
                # interim per-row flag is OVERWRITTEN below by the BH-FDR-adjusted
                # coupling_survives (F05); no magnitude threshold anywhere.
                collapsed = (p_par >= 0.05) or (np.sign(rho_par) != np.sign(rho_raw))
                rows.append(dict(cohort=gse, tf=tf, population=pop,
                                 rho_raw=round(float(rho_raw), 3), p_raw=float(p_raw),
                                 rho_partial=round(float(rho_par), 3), p_partial=float(p_par),
                                 n=int(n_used),
                                 confounders=';'.join(conf.columns),
                                 coupling_survives=(not collapsed)))

    B = pd.DataFrame(rows)
    # ---- BH-FDR across ALL partial tests (tf x population x cohort) BEFORE
    #      deciding survival; recompute coupling_survives on adjusted p ----
    if not B.empty:
        from statsmodels.stats.multitest import multipletests
        _ok = B['p_partial'].notna()
        B['padj_partial'] = np.nan
        if _ok.sum() > 0:
            B.loc[_ok, 'padj_partial'] = multipletests(
                B.loc[_ok, 'p_partial'].values, method='fdr_bh')[1]
        _okr = B['p_raw'].notna()
        B['padj_raw'] = np.nan
        if _okr.sum() > 0:
            B.loc[_okr, 'padj_raw'] = multipletests(
                B.loc[_okr, 'p_raw'].values, method='fdr_bh')[1]
        # survival = FDR-significant partial assoc with sign preserved (see F03)
        B['coupling_survives'] = (B['padj_partial'] < 0.05) & \
            (np.sign(B['rho_partial']) == np.sign(B['rho_raw']))
    B.to_csv(out_csv, index=False)
    log('STAGE B rows:', len(B), '-> (BH-FDR over all partial tests) ->', out_csv)
    if not B.empty:
        surv = B.groupby(['tf']).agg(n=('population', 'size'),
                                     survived=('coupling_survives', 'sum'))
        log('  coupling survival after adjustment (survived / tested):')
        for tf, r in surv.iterrows():
            log(f'    {tf}: {int(r.survived)} / {int(r.n)} population-links survive '
                f'purity+load+IFN adjustment')
        # HONEST verdict line
        frac_surv = B['coupling_survives'].mean()
        log(f'  VERDICT: {frac_surv:.0%} of raw immune-regulon correlations survive '
            f'deconfounding. If low, the coupling is largely a shared-cellularity/purity '
            f'artifact (report honestly; do NOT claim an immune->HES1 axis).')
        _fig_stageB(B)
    return B

def _partial_spearman(x, y, conf):
    """Partial Spearman rho(x,y | conf) via rank-transform + OLS residualization.
    Prefer pingouin.partial_corr if available (spearman method); else manual."""
    d = pd.concat([x.rename('x'), y.rename('y'), conf], axis=1).dropna()
    covars = [c for c in conf.columns if c in d.columns and d[c].notna().any()]
    if len(d) < max(10, len(covars) + 5):
        return np.nan, np.nan, len(d)
    import pingouin as pg  # required dep; missing -> hard-fail (rule #1)
    r = pg.partial_corr(data=d, x='x', y='y', covar=covars, method='spearman')
    return float(r['r'].iloc[0]), float(r['p-val'].iloc[0]), len(d)

def _reindex_deconv(frac, gse):
    """★2026-07-08: IOBR deconvo_tme returned POSITIONAL integer IDs (1,2,3..) instead of sample
    names for some cohorts (e.g. GSE58095) -> the deconv index did not match the GSM-indexed regulon
    activity -> STAGE B coupling aligned on 0 common samples -> 0 rows for EVERY cohort. Re-attach the
    true sample IDs from the expr matrix that was deconvolved (deconvo_tme preserves eset column order)."""
    ep = os.path.join(CFG['OUT'], 'stageB', f'{gse}_expr_for_deconv.csv')
    if os.path.exists(ep):
        cols = list(pd.read_csv(ep, index_col=0, nrows=0).columns)
        if len(cols) == len(frac):
            frac = frac.copy(); frac.index = cols
    return frac


def _deconvolve(co, gse, force):
    """IOBR 7-method deconvolution + ESTIMATE (module 492) via Rscript.
    Writes a TPM-ish gene x sample matrix, calls R, reads back fractions.
    FALLBACK if Rscript/IOBR unavailable: marker-mean MCP-counter (adapted from
    mcp_immune.py) + an ESTIMATE-style immune/stromal z-score (clearly labelled)."""
    cache = os.path.join(CFG['OUT'], 'stageB', f'{gse}_deconv.csv')
    if os.path.exists(cache) and not force:
        log('  reuse cached deconvolution:', gse)
        return _reindex_deconv(pd.read_csv(cache, index_col=0), gse)

    # rebuild expr for deconvolution (needs raw gene expression, not activity)
    if not os.path.exists(os.path.join(CFG['RAW'], co['file'])):
        return None
    expr, meta = load_cohort(co)
    keep = meta.index[meta['grp'].isin(['SSc', 'HC'])]
    expr = expr[keep]
    # ★2026-07-07 修(P0): 线性模型去卷积(CIBERSORT/EPIC/quanTIseq/ESTIMATE)要求线性 TPM/FPKM,
    # 而 loader 输出 log(RNA-seq=log1p / array=log2)。按数值域检测 log 尺度并反变换回线性
    # (防对已线性数据重复反变换; 否则免疫比例+ESTIMATE 纯度全偏 -> StageB 去混杂静默失效)。
    is_array = not str(co.get('platform', '')).lower().startswith('rnaseq')
    if float(np.nanmax(expr.values)) < 30.0:                 # <30 视为 log 尺度 -> 反变换
        expr = (np.power(2.0, expr) if is_array else np.expm1(expr)).clip(lower=0)
    expr_path = os.path.join(CFG['OUT'], 'stageB', f'{gse}_expr_for_deconv.csv')
    expr.to_csv(expr_path)

    # ---- IOBR via R (module 492 real API: deconvo_tme) is MANDATORY. ----
    # No proxy substitution: a marker-mean surrogate would fabricate the very
    # ESTIMATE-purity confounder STAGE B exists to adjust for. Hard-fail so the
    # env gets fixed (rules #1, #3).
    rscript = _write_iobr_r(gse, expr_path, cache, is_array)
    ok = _try_rscript(rscript)
    if not (ok and os.path.exists(cache)):
        raise RuntimeError(
            f'IOBR deconvolution FAILED for {gse} (Rscript rc!=0 or no output at '
            f'{cache}). STAGE B requires real IOBR deconvo_tme + ESTIMATE purity; '
            f'refusing to substitute a marker-mean surrogate. Install IOBR (module '
            f'492) in the R env pointed to by CFG["RSCRIPT"] and re-run.')
    return _reindex_deconv(pd.read_csv(cache, index_col=0), gse)

def _write_iobr_r(gse, expr_path, out_path, is_array=True):
    """emit an R driver that runs IOBR module-492 methods (mcpcounter/xcell/epic/
    quantiseq/timer/cibersort/estimate) and writes a merged fraction table +
    ESTIMATE scores. Idempotent target = out_path."""
    r_path = os.path.join(CFG['OUT'], 'stageB', f'run_iobr_{gse}.R')
    r = f'''# auto-generated IOBR deconvolution (module 492) for {gse}  [SEED {CFG['SEED']}]
suppressMessages({{library(IOBR); library(tidyverse)}})
set.seed({CFG['SEED']})
eset <- as.matrix(read.csv("{expr_path.replace(os.sep, '/')}", row.names=1, check.names=FALSE))
# IOBR expects genes(rows, symbol) x samples(cols); ★线性 TPM/FPKM(已在 _deconvolve 反 log)。
methods <- c("mcpcounter","xcell","epic","quantiseq","timer","cibersort","estimate")
decon <- lapply(methods, function(m)
  tryCatch(deconvo_tme(eset=eset, method=m, arrays={'TRUE' if is_array else 'FALSE'}, perm=1000),
           error=function(e){{message(m," FAILED: ",conditionMessage(e)); NULL}}))
names(decon) <- methods
failed <- methods[vapply(decon, is.null, logical(1))]
# ESTIMATE is REQUIRED: it supplies the purity/immune/stromal confounder that
# STAGE B adjusts for. If it silently failed, STAGE B would report "coupling
# survives" WITHOUT actually adjusting for purity -> hard-fail (rule #1/#3).
if ("estimate" %in% failed)
  stop("ESTIMATE deconvolution FAILED; STAGE B purity confounder unavailable. Abort.")
decon <- Filter(Negate(is.null), decon)
if (length(decon)==0) quit(status=3)
if (length(failed)>0) message("NOTE non-fatal methods failed: ", paste(failed, collapse=","))
tme <- purrr::reduce(decon, ~ dplyr::inner_join(.x, .y, by="ID"))
rownames(tme) <- tme$ID; tme$ID <- NULL
write.csv(tme, "{out_path.replace(os.sep, '/')}")
cat("IOBR done:", ncol(tme), "features x", nrow(tme), "samples\\n")
'''
    with open(r_path, 'w', encoding='utf-8') as f:
        f.write(r)
    return r_path

def _try_rscript(r_path):
    try:
        p = subprocess.run([CFG['RSCRIPT'], r_path], capture_output=True, text=True,
                           timeout=3600)
        log('  Rscript rc=', p.returncode, '|', (p.stdout or p.stderr)[-200:].strip())
        return p.returncode == 0
    except Exception as e:
        log('  [Rscript call failed]', str(e)[:120]); return False

def _mcp_fallback(expr, cache):
    raise RuntimeError(
        '_mcp_fallback is DISABLED: STAGE B must use real IOBR/ESTIMATE. A '
        'marker-mean surrogate + rank-based pseudo-purity would fabricate the '
        'deconfounding covariate and produce silently-wrong coupling verdicts.')
    # minimal built-in MCP marker sets (subset; full set lives in ml_validation/mcpcounter_genes.txt)
    MCP = {
        'T cells':        ['CD3D', 'CD3E', 'CD3G', 'CD28', 'CD5', 'CD6', 'TRAT1'],
        'CD8 T cells':    ['CD8A', 'CD8B'],
        'Cytotoxic lymph':['GNLY', 'KLRD1', 'NKG7', 'GZMB', 'EOMES'],
        'B lineage':      ['CD19', 'CD79A', 'MS4A1', 'BANK1', 'CD22'],
        'NK cells':       ['KIR2DL1', 'KIR2DL3', 'NCR1'],
        'Monocytic':      ['CD14', 'CD68', 'CSF1R', 'LYZ'],
        'Myeloid DC':     ['CD1A', 'CD1B', 'CLEC10A'],
        'Neutrophils':    ['FCGR3B', 'CSF3R', 'S100A12'],
        'Endothelial':    ['PECAM1', 'VWF', 'CDH5', 'CLDN5'],
        'Fibroblasts':    ['COL1A1', 'COL3A1', 'PDGFRA', 'LUM', 'DCN'],
    }
    frac = {}
    for pop, genes in MCP.items():
        g = [x for x in genes if x in expr.index]
        if len(g) >= 2:
            frac[pop] = expr.loc[g].apply(pd.to_numeric, errors='coerce').mean(axis=0)
    F = pd.DataFrame(frac)
    # ESTIMATE-style surrogates (z-score means) — labelled *_fallback
    F['ImmuneScore_fallback'] = zscore_mean(expr, CFG['IMMUNE_GENES']).reindex(F.index)
    F['StromalScore_fallback'] = zscore_mean(expr, CFG['FIBRO_GENES']).reindex(F.index)
    # crude purity surrogate: -(immune+stromal) rank -> higher = purer parenchyma
    F['TumorPurity_fallback'] = -(F['ImmuneScore_fallback'].rank()
                                  + F['StromalScore_fallback'].rank())
    F.to_csv(cache)
    return F

# =====================================================================
# STAGE C : external-cohort ML (reuse 059/Mime driver). OPTIONAL (--ml).
# =====================================================================
def stage_C(force):
    """Re-export the two regulon-activity feature tables the existing
    ml_consensus_driver.R expects (train=GSE130955, external=GSE58095), then
    shell out to Rscript. We DO NOT re-implement the 100+ combo ML — we reuse it.
    The driver + ml_059_functions.R must be copied next to the feature files."""
    feat = os.path.join(CFG['OUT'], 'features')
    train_gse, ext_gse = 'GSE130955', 'GSE58095'
    need = [f'{train_gse}_regulon_activity.csv', f'{ext_gse}_regulon_activity.csv']
    for n in need:
        if not os.path.exists(os.path.join(feat, n)):
            log('[STAGE C] missing', n, '- run STAGE A first'); return None

    # driver expects *_labels.csv (0/1). derive from cached meta (grp).
    for gse in [train_gse, ext_gse]:
        meta = pd.read_csv(os.path.join(feat, f'{gse}_meta.csv'), index_col=0)
        lab = (meta['grp'] == 'SSc').astype(int)
        lab.rename('SSc').to_csv(os.path.join(feat, f'{gse}_labels.csv'))
    log('[STAGE C] labels exported. Now run the 059/Mime driver:')
    log('  cp ml_validation/ml_consensus_driver.R ml_validation/ml_059_functions.R '
        f'{feat}/ ; cd {feat} && {CFG["RSCRIPT"]} ml_consensus_driver.R')
    # attempt automatically if the driver is already present next to features
    drv = os.path.join(feat, 'ml_consensus_driver.R')
    if os.path.exists(drv) and os.path.exists(os.path.join(feat, 'ml_059_functions.R')):
        if not _try_rscript(drv):
            raise RuntimeError('STAGE C: 059/Mime ML driver returned non-zero. '
                               'Do not report STAGE C as complete on a failed run.')
    else:
        log('  (driver not copied yet -> STAGE C left as an explicit manual step; '
            'NOT counted as completed)')
        return False
    return True

# =====================================================================
# FIGURES  (NO plain bar charts — dumbbell / dot / lollipop / heatmap)
# =====================================================================
def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'figure.dpi': 200, 'font.size': 9, 'svg.fonttype': 'none'})
    return plt

def _fig_stageA(R):
    if R.empty:
        return
    plt = _mpl()
    sub = R[R.analysis == 'SSc_vs_HC']
    if sub.empty:
        return
    # DUMBBELL: per TF, one dot per cohort (effect), colored by significance.
    tfs = ['SMAD3', 'FOSB', 'HES1', 'NOTCH_EFFECTOR', 'CEBPD', 'MEF2C']
    tfs = [t for t in tfs if t in sub.tf.unique()]
    fig, ax = plt.subplots(figsize=(7.2, 0.62 * len(tfs) + 1.4))
    cohorts = sorted(sub.cohort.unique())
    cmap = plt.get_cmap('tab10')
    ccol = {c: cmap(i % 10) for i, c in enumerate(cohorts)}
    for yi, tf in enumerate(tfs):
        g = sub[sub.tf == tf]
        ax.plot([g.effect.min(), g.effect.max()], [yi, yi], color='#bbbbbb', lw=1.2, zorder=1)
        for _, row in g.iterrows():
            ax.scatter(row.effect, yi, s=70 if row.p < 0.05 else 38,
                       color=ccol[row.cohort],
                       edgecolor='k' if row.p < 0.05 else 'none', lw=0.7, zorder=3)
    ax.axvline(0, color='k', lw=0.8, ls='--')
    ax.set_yticks(range(len(tfs))); ax.set_yticklabels(tfs)
    ax.set_xlabel('SSc − HC regulon activity (Welch Δ); filled+outlined = p<0.05')
    ax.set_title('Bulk regulon activity, SSc vs HC (per-cohort dumbbell)')
    handles = [plt.Line2D([0], [0], marker='o', ls='', color=ccol[c], label=c) for c in cohorts]
    ax.legend(handles=handles, fontsize=7, loc='center left', bbox_to_anchor=(1.01, 0.5),
              frameon=False)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(CFG['OUT'], 'figs', f'fig05A_regulon_ssc_dumbbell.{ext}'),
                    bbox_inches='tight')
    plt.close(fig)

    # DOT: mRSS correlations (rho) per TF/cohort
    mr = R[R.analysis == 'mRSS']
    if not mr.empty:
        fig, ax = plt.subplots(figsize=(5.5, 0.5 * mr.tf.nunique() + 1.3))
        for yi, tf in enumerate(sorted(mr.tf.unique())):
            g = mr[mr.tf == tf]
            for _, row in g.iterrows():
                ax.scatter(row.effect, yi, s=90 if row.p < 0.05 else 45,
                           color='#C0392B' if row.effect > 0 else '#2C6FBB',
                           edgecolor='k' if row.p < 0.05 else 'none', lw=0.7)
                ax.annotate(row.cohort, (row.effect, yi), fontsize=6,
                            xytext=(0, 7), textcoords='offset points', ha='center')
        ax.axvline(0, color='k', lw=0.8, ls='--')
        ax.set_yticks(range(mr.tf.nunique())); ax.set_yticklabels(sorted(mr.tf.unique()))
        ax.set_xlabel('Spearman ρ vs mRSS (filled+outlined = p<0.05)')
        ax.set_title('Regulon activity vs skin severity (mRSS)')
        fig.tight_layout()
        for ext in ('png', 'pdf'):
            fig.savefig(os.path.join(CFG['OUT'], 'figs', f'fig05A_mrss_dot.{ext}'),
                        bbox_inches='tight')
        plt.close(fig)

def _fig_stageB(B):
    if B.empty:
        return
    plt = _mpl()
    # DUMBBELL raw->partial rho per (population), faceted by TF, pooled across cohorts (median)
    for tf in ['HES1', 'FOSB']:
        g = B[B.tf == tf]
        if g.empty:
            continue
        piv = g.groupby('population').agg(raw=('rho_raw', 'median'),
                                          par=('rho_partial', 'median'),
                                          surv=('coupling_survives', 'mean')).reset_index()
        piv = piv.sort_values('raw')
        fig, ax = plt.subplots(figsize=(6.6, 0.44 * len(piv) + 1.3))
        y = range(len(piv))
        ax.hlines(list(y), piv.raw, piv.par, color='#cccccc', lw=1.4, zorder=1)
        ax.scatter(piv.raw, y, s=55, color='#888888', label='raw ρ', zorder=3)
        # Colour by whether the FDR-based survival flag holds in the MAJORITY of
        # cohorts; label the 0.5 majority rule explicitly rather than hiding it. The
        # real per-test survival decision lives in the BH-FDR coupling_survives (F05).
        ax.scatter(piv.par, y, s=70,
                   color=['#2ca02c' if s > 0.5 else '#d62728' for s in piv.surv],
                   edgecolor='k', lw=0.6,
                   label='partial ρ (green = survives in majority of cohorts)', zorder=4)
        ax.axvline(0, color='k', lw=0.8, ls='--')
        ax.set_yticks(list(y)); ax.set_yticklabels(piv.population)
        ax.set_xlabel(f'Spearman ρ ( {tf} regulon activity vs immune fraction )')
        ax.set_title(f'{tf}: immune coupling raw → deconfounded\n'
                     f'(adj: ESTIMATE purity + fibro/immune load + IFN)')
        ax.legend(fontsize=7, frameon=False, loc='lower right')
        fig.tight_layout()
        for ext in ('png', 'pdf'):
            fig.savefig(os.path.join(CFG['OUT'], 'figs',
                        f'fig05B_{tf}_coupling_deconfounded.{ext}'), bbox_inches='tight')
        plt.close(fig)

# =====================================================================
# CLI
# =====================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stage', default='AB', help='which stages: any of A,B,C,all (e.g. A, AB, all)')
    ap.add_argument('--ml', action='store_true', help='run STAGE C external ML (059/Mime driver)')
    ap.add_argument('--force', action='store_true', help='ignore checkpoints, recompute')
    args = ap.parse_args()
    ensure_out()
    log('CONFIG:', json.dumps({k: v for k, v in CFG.items()
                               if k in ('RAW', 'OUT', 'NET_CSV', 'RSCRIPT', 'SEED')}))
    stages = 'ABC' if args.stage.lower() == 'all' else args.stage.upper()
    if 'A' in stages:
        stage_A(args.force)
    if 'B' in stages:
        stage_B(args.force)
    if ('C' in stages) or args.ml:
        stage_C(args.force)
    log('DONE. outputs under', CFG['OUT'])

if __name__ == '__main__':
    main()
