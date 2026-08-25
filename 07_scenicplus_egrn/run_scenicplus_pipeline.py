#!/usr/bin/env python
"""
SCENIC+ eGRN pipeline for GSE312129 SSc whole-skin snMULTIOME (14 donors: SSC1-10, HC1-4).
INDEPENDENT second-engine cross-check for the SSc virtual-perturbation paper.
Question: does HES1 / Notch-axis (HEY1/HEYL/HES-family/RBPJ) emerge as an eRegulon
active in the dermal FIBROBLAST/myofibroblast compartment?

Env : source /data/ssc/miniconda3/etc/profile.d/conda.sh && conda activate scenicplus
      scenicplus 1.0a2 / pycisTopic 2.0a0 / pycistarget 1.1 ; macs2 ; Mallet-202108 ; openjdk11
Cores: HARD CAP 6 (base-GRN/CellOracle runs in parallel on the other cores).
Data : /data/ssc/basegrn/frags/GSM*_{sample}_atac_fragments.tsv.gz (+.tbi)  (hg38)
       /data/ssc/basegrn/frags/GSM*_{sample}_filtered_feature_bc_matrix.h5 (multiome GEX+peaks)

STAGE-DRIVEN & CHECKPOINTED (crash-resilient; each stage writes an output, next stage loads it):
  python run_scenicplus_pipeline.py rna        # FAST, network-free  -> can run during audit
  python run_scenicplus_pipeline.py peaks      # LONG: pseudobulk + MACS2 + consensus  (needs audit)
  python run_scenicplus_pipeline.py cistopic   # cisTopic object from consensus peaks
  python run_scenicplus_pipeline.py topics      # LONG: Mallet LDA topic models (n_cpu=6)
  python run_scenicplus_pipeline.py regionsets  # binarize topics + DARs -> region_sets/*.bed
  python run_scenicplus_pipeline.py scenicplus  # LONG: snakemake motif-enrich + eGRN inference
  python run_scenicplus_pipeline.py extract      # HES1/Notch eRegulon per cell type

Parameter choices are TUTORIAL DEFAULTS. Source URLs are cited inline (SRC:).
Tutorials followed:
  pycisTopic ATAC pp: https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html
  SCENIC+ scRNA pp  : https://scenicplus.readthedocs.io/en/latest/human_cerebellum_scRNA_pp.html
  SCENIC+ workflow  : https://scenicplus.readthedocs.io/en/latest/human_cerebellum.html
  pycisTopic API    : https://pycistopic.readthedocs.io/en/latest/api.html
"""
import os, sys, glob, gc, re, pickle, subprocess
import matplotlib; matplotlib.use("Agg")     # headless (nohup) — no display
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------- CONFIG
BASE      = "/data/ssc/scenicplus"
FRAG_DIR  = "/data/ssc/basegrn/frags"          # READ-ONLY (owned by base-GRN); we only READ frags/h5
RNA_DIR   = f"{BASE}/rna"
QC_DIR    = f"{BASE}/qc"
ATAC_DIR  = f"{BASE}/atac"
RS_DIR    = f"{BASE}/region_sets"
CTX_DIR   = f"{BASE}/ctx_db"
REF_DIR   = f"{BASE}/ref"
SCPLUS_DIR= f"{BASE}/scplus"
TMP_DIR   = f"{BASE}/tmp"
LOG_DIR   = f"{BASE}/logs"
for d in (RNA_DIR,QC_DIR,ATAC_DIR,RS_DIR,SCPLUS_DIR,TMP_DIR,LOG_DIR):
    os.makedirs(d, exist_ok=True)

N_CPU        = 6      # cistopic BUILD stays SERIAL (n_cpu=1 per sample); this = pseudobulk/general cap
N_CPU_COMPUTE= 8      # OOM-fix: heavier COMPUTE stages (Mallet, find_diff_features, snakemake) — memory-safe, ~8 to coexist with base-GRN
CISTOPIC_PARTITION = 20  # OOM-fix: chunk fragment reading in create_cistopic_object (was default 5) -> lower peak RAM on big frag files (e.g. SSC7 7.4G)
MIN_DONORS_PER_REGION = 2 # OOM/de-noise: keep consensus peaks reproducible across >=2 donors (methods-documented)
SEED      = 555                                 # SRC: tutorial random_state=555
SEP       = "___"                               # pycisTopic 2.0 default split_pattern (barcode___sample_id)
MACS_PATH = "macs2"
MALLET    = f"{BASE}/Mallet-202108/bin/mallet"
CHROMSIZES= f"{REF_DIR}/hg38.chrom.sizes"
BLACKLIST = f"{REF_DIR}/hg38-blacklist.v2.bed"
RANKINGS_DB = f"{CTX_DIR}/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
SCORES_DB   = f"{CTX_DIR}/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
MOTIF2TF    = f"{CTX_DIR}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

ADATA_ANNOT = f"{RNA_DIR}/adata_annot.h5ad"
CISTOPIC_PKL= f"{ATAC_DIR}/cistopic_obj.pkl"
MODELS_PKL  = f"{ATAC_DIR}/topic_models.pkl"
CISTOPIC_MODEL_PKL = f"{ATAC_DIR}/cistopic_obj_models.pkl"
CONSENSUS_BED = f"{ATAC_DIR}/consensus_regions.bed"

# Coarse skin cell-type marker sets (from task spec). SRC: task grounding.
MARKERS = {
    "Fibroblast"  : ["COL1A1","PDGFRA","LUM"],
    "Myeloid"     : ["CD68","LYZ"],
    "Endothelial" : ["PECAM1","VWF"],
    "Keratinocyte": ["KRT14","KRT5"],
    "Mural"       : ["ACTA2","RGS5"],
    "Lymphocyte"  : ["CD3D","CD3E"],
}
# TFs of interest for the scientific question
NOTCH_TFS = ["HES1","HES2","HES4","HES5","HES6","HES7","HEY1","HEY2","HEYL","RBPJ","RBPJL"]

# --------------------------------------------------------------------- sample discovery
def discover_samples():
    """Map sample_id (SSC1..HC4) -> (h5 path, fragments path). h5 & atac have different GSM ids."""
    frags = {}; h5s = {}
    for f in glob.glob(f"{FRAG_DIR}/GSM*_atac_fragments.tsv.gz"):
        m = re.search(r"GSM\d+_([A-Za-z]+\d+)_atac_fragments", os.path.basename(f)); frags[m.group(1)] = f
    for f in glob.glob(f"{FRAG_DIR}/GSM*_filtered_feature_bc_matrix.h5"):
        m = re.search(r"GSM\d+_([A-Za-z]+\d+)_filtered", os.path.basename(f)); h5s[m.group(1)] = f
    samples = sorted(set(frags) & set(h5s), key=lambda s:(s[:2], int(re.search(r'\d+',s).group())))
    assert len(samples) == 14, f"expected 14 samples, got {len(samples)}: {samples}"
    return {s: (h5s[s], frags[s]) for s in samples}

# ============================================================================ STAGE rna
def stage_rna():
    """Load 14 multiome RNA matrices -> QC -> scrublet -> normalize -> harmony -> leiden ->
       marker-based coarse annotation. Saves adata with obs_names = barcode___sample_id,
       .X = lognorm, .raw = raw counts, obs['cell_type'] & obs['sample_id'] (SCENIC+-ready).
       SRC scRNA pp: https://scenicplus.readthedocs.io/en/latest/human_cerebellum_scRNA_pp.html"""
    import scanpy as sc
    import scrublet as scr
    sc.settings.n_jobs = N_CPU
    samples = discover_samples()
    ads = []
    for s,(h5,_) in samples.items():
        a = sc.read_10x_h5(h5, gex_only=True)      # multiome h5 -> keep Gene Expression only
        a.var_names_make_unique()
        a.obs["sample_id"] = s
        a.obs_names = [f"{bc}{SEP}{s}" for bc in a.obs_names]   # barcode___sample_id (unique, matches ATAC)
        # per-sample scrublet doublet removal (run per sample; never across samples)
        try:
            sob = scr.Scrublet(a.X, random_state=SEED)
            ds, pred = sob.scrub_doublets(min_counts=2, min_cells=3, verbose=False)
            a.obs["doublet_score"] = ds; a.obs["predicted_doublet"] = pred
        except Exception as e:
            a.obs["predicted_doublet"] = False
            print(f"[scrublet warn {s}] {e}")
        ads.append(a)
        print(f"[rna] {s}: {a.n_obs} cells")
    import anndata as ad
    adata = ad.concat(ads, join="outer", index_unique=None); del ads; gc.collect()
    adata.obs_names_make_unique()
    # ---- QC (SRC: scanpy standard + tutorial calculate_qc_metrics with mt) ----
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)       # SRC: scanpy pbmc3k standard
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["pct_counts_mt"] < 15].copy()       # nuclei: lenient mito cutoff
    adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
    print(f"[rna] cells {n0} -> {adata.n_obs} after QC+doublet")
    # ---- preserve raw counts, then normalize (SRC: tutorial adata.raw=adata BEFORE norm) ----
    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=1e4)   # SRC: tutorial
    sc.pp.log1p(adata)                             # SRC: tutorial  -> .X now lognorm (SCENIC+ input)
    # ---- clustering on a SEPARATE scaled copy (keep adata.X lognorm for SCENIC+) ----
    acl = adata.copy()
    sc.pp.highly_variable_genes(acl, min_mean=0.0125, max_mean=3, min_disp=0.5)  # SRC: tutorial
    acl = acl[:, acl.var.highly_variable].copy()
    sc.pp.scale(acl, max_value=10)                 # SRC: tutorial
    sc.tl.pca(acl, n_comps=50, random_state=SEED)
    # harmony batch-integration across 14 donors (harmonypy) before neighbors
    try:
        import scanpy.external as sce
        sce.pp.harmony_integrate(acl, "sample_id", max_iter_harmony=20)
        rep = "X_pca_harmony"
    except Exception as e:
        print(f"[harmony warn] {e}; using X_pca"); rep = "X_pca"
    sc.pp.neighbors(acl, n_neighbors=15, use_rep=rep, random_state=SEED)
    sc.tl.leiden(acl, resolution=1.0, random_state=SEED)   # coarse resolution
    sc.tl.umap(acl, random_state=SEED)
    adata.obs["leiden"] = acl.obs["leiden"].values
    adata.obsm["X_umap"] = acl.obsm["X_umap"]
    # ---- marker-based coarse annotation: score each cell-type gene set, assign cluster by argmax ----
    for ct, genes in MARKERS.items():
        g = [x for x in genes if x in adata.var_names]
        sc.tl.score_genes(adata, g, score_name=f"score_{ct}", use_raw=False)
    score_cols = [f"score_{ct}" for ct in MARKERS]
    clmean = adata.obs.groupby("leiden")[score_cols].mean()
    cl2ct = {cl: score_cols[int(np.argmax(row.values))].replace("score_","")
             for cl,row in clmean.iterrows()}
    adata.obs["cell_type"] = adata.obs["leiden"].map(cl2ct).astype("category")
    # ---- save + QC tables ----
    adata.write(ADATA_ANNOT)
    adata.obs["cell_type"].value_counts().to_csv(f"{QC_DIR}/celltype_counts.csv")
    (adata.obs.groupby(["cell_type","sample_id"]).size().unstack(fill_value=0)
        .to_csv(f"{QC_DIR}/celltype_by_sample.csv"))
    clmean.assign(assigned=[cl2ct[c] for c in clmean.index]).to_csv(f"{QC_DIR}/cluster_marker_scores.csv")
    print("[rna] cell_type counts:\n", adata.obs["cell_type"].value_counts())
    print(f"[rna] saved {ADATA_ANNOT}")

# ========================================================================== STAGE peaks
def get_chromsizes():
    import pyranges as pr
    cs = pd.read_table(CHROMSIZES, header=None, names=["Chromosome","End"]); cs["Start"]=0
    return pr.PyRanges(cs[["Chromosome","Start","End"]])

def stage_peaks():
    """Pseudobulk per cell type (pool all 14 donors) -> MACS2 -> consensus peaks.
       SRC: https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html"""
    import scanpy as sc
    from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
    from pycisTopic.iterative_peak_calling import get_consensus_peaks
    samples = discover_samples()
    fragments_dict = {s: fp for s,(_,fp) in samples.items()}       # sample_id -> fragments path
    adata = sc.read_h5ad(ADATA_ANNOT)
    # cell_data index = barcode___sample_id (matches SEP); export_pseudobulk infers sample via split_pattern
    cell_data = adata.obs[["cell_type","sample_id"]].copy()
    chromsizes = get_chromsizes()
    os.makedirs(f"{ATAC_DIR}/consensus_peak_calling/pseudobulk_bed_files", exist_ok=True)
    os.makedirs(f"{ATAC_DIR}/consensus_peak_calling/pseudobulk_bw_files", exist_ok=True)
    # (A) pseudobulk BED/BW per cell type. SRC: export_pseudobulk tutorial defaults.
    bw_paths, bed_paths = export_pseudobulk(
        input_data=cell_data, variable="cell_type", sample_id_col="sample_id",
        chromsizes=chromsizes,
        bed_path=f"{ATAC_DIR}/consensus_peak_calling/pseudobulk_bed_files",
        bigwig_path=f"{ATAC_DIR}/consensus_peak_calling/pseudobulk_bw_files",
        path_to_fragments=fragments_dict, n_cpu=N_CPU,
        normalize_bigwig=True, split_pattern=SEP, temp_dir=TMP_DIR)
    with open(f"{ATAC_DIR}/bed_paths.pkl","wb") as f: pickle.dump(bed_paths,f)
    # (B) MACS2 per pseudobulk. SRC: tutorial input_format BEDPE, shift 73, ext_size 146, keep_dup all, q 0.05.
    narrow = peak_calling(
        macs_path=MACS_PATH, bed_paths=bed_paths,
        outdir=f"{ATAC_DIR}/consensus_peak_calling/MACS", genome_size="hs",
        n_cpu=N_CPU, input_format="BEDPE", shift=73, ext_size=146,
        keep_dup="all", q_value=0.05)   # AUDIT FIX: pycisTopic 2.0a0 peak_calling has NO _temp_dir (was 1.x remnant, fell into **kwargs) -> removed; MACS2 uses /tmp (87G free)
    # (C) consensus peaks. SRC: tutorial peak_half_width=250 + hg38 blacklist.
    consensus = get_consensus_peaks(
        narrow_peaks_dict=narrow, peak_half_width=250,
        chromsizes=chromsizes, path_to_blacklist=BLACKLIST)
    consensus.to_bed(CONSENSUS_BED)
    print(f"[peaks] consensus regions: {len(consensus)} -> {CONSENSUS_BED}")

# ======================================================================= STAGE cistopic
def stage_cistopic():
    """Build merged cisTopic object from fragments + consensus peaks.
       valid_bc per sample = RNA-annotated multiome barcodes (10x multiome: RNA calls the cells).
       OOM-FIX: build SERIALLY (n_cpu=1), chunk fragment reading (partition=20), gc.collect() after each
       sample to release the per-sample fragment-parsing transient; then filter regions to those
       reproducible across >=MIN_DONORS_PER_REGION donors (de-noise; lowers downstream impute memory;
       preserves cell-type/fibroblast specificity since real enhancers recur across donors — methods-documented).
       SRC: create_cistopic_object_from_fragments / merge  (pycisTopic API)."""
    import scanpy as sc, gc
    import numpy as np
    from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
    samples = discover_samples()
    adata = sc.read_h5ad(ADATA_ANNOT)
    objs = []
    for s,(_,fp) in samples.items():
        # raw barcodes for this sample (strip the ___sample_id tag -> matches fragment file barcodes)
        bcs = [n.split(SEP)[0] for n in adata.obs_names[adata.obs["sample_id"]==s]]
        if len(bcs)==0: continue
        obj = create_cistopic_object_from_fragments(
            path_to_fragments=fp, path_to_regions=CONSENSUS_BED,
            path_to_blacklist=BLACKLIST, valid_bc=bcs, n_cpu=1,
            project=s, split_pattern=SEP, partition=CISTOPIC_PARTITION)   # cell names -> barcode___{s}
        objs.append(obj); print(f"[cistopic] {s}: {len(obj.cell_names)} cells", flush=True); gc.collect()
    cistopic_obj = merge(objs, split_pattern=SEP) if len(objs)>1 else objs[0]
    del objs; gc.collect()
    # attach RNA cell-type metadata (align on barcode___sample_id)
    meta = adata.obs[["cell_type","sample_id"]].copy()
    common = [c for c in cistopic_obj.cell_names if c in meta.index]
    cistopic_obj.add_cell_data(meta.loc[common], split_pattern=SEP)
    # ---- REGION FILTER: keep peaks accessible in cells from >=MIN_DONORS_PER_REGION donors ----
    fm   = cistopic_obj.fragment_matrix                                   # regions x cells (sparse)
    samp = cistopic_obj.cell_data.loc[cistopic_obj.cell_names, "sample_id"].astype(str).values
    binm = (fm > 0)
    n_donors = np.zeros(fm.shape[0], dtype=int)
    for s in np.unique(samp):
        cols = np.where(samp == s)[0]
        n_donors += (np.asarray(binm[:, cols].sum(axis=1)).ravel() > 0).astype(int)
    keep = n_donors >= MIN_DONORS_PER_REGION
    kept = [r for r,k in zip(cistopic_obj.region_names, keep) if k]
    print(f"[cistopic] region filter >= {MIN_DONORS_PER_REGION} donors: {fm.shape[0]} -> {len(kept)} regions", flush=True)
    cistopic_obj = cistopic_obj.subset(cells=None, regions=kept, copy=True)
    with open(CISTOPIC_PKL,"wb") as f: pickle.dump(cistopic_obj,f)
    print(f"[cistopic] final cells={len(cistopic_obj.cell_names)} regions={len(cistopic_obj.region_names)} -> {CISTOPIC_PKL}")

# ========================================================================= STAGE topics
def stage_topics():
    """Mallet LDA topic models (parallel). SPEED TRADE-OFF (methods-documented): the tutorial's full grid
       [2,5,10,15,20,25,30,35,40,45,50] x 500 iters would take ~13-20h on 48113 cells x 407852 regions.
       Trimmed to a 6-point grid spanning 2-50 (ample span for evaluate_models auto-selection) with 250
       iters (pycisTopic uses 150-500; 250 converges well). alpha/eta unchanged; auto model selection unchanged.
       SRC grid/params rationale: https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html"""
    from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models
    os.environ["MALLET_MEMORY"] = "100G"
    with open(CISTOPIC_PKL,"rb") as f: cistopic_obj = pickle.load(f)
    mallet_tmp = f"{TMP_DIR}/mallet"; os.makedirs(mallet_tmp, exist_ok=True)
    models = run_cgs_models_mallet(
        cistopic_obj,
        n_topics=[2,10,20,30,40,50],                          # SPEED: 6-point grid (was 11), full 2-50 span
        n_cpu=10, n_iter=250, random_state=SEED,              # SPEED: n_iter 500->250, Mallet 10 cores (base-GRN=1 core)
        alpha=50, alpha_by_topic=True, eta=0.1, eta_by_topic=False,
        tmp_path=mallet_tmp, save_path=mallet_tmp, mallet_path=MALLET)
    with open(MODELS_PKL,"wb") as f: pickle.dump(models,f)
    # Model selection: NO hard-coded topic count. evaluate_models(select_model=None) auto-selects
    # the best model by Mimno_2011 topic coherence (verified in pycisTopic 2.0 source, line 156-163).
    # SRC: https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html
    model = evaluate_models(models, select_model=None, return_model=True, plot=False)
    assert model is not None, "evaluate_models returned None — model selection failed"  # AUDIT GUARD
    cistopic_obj.add_LDA_model(model)
    with open(CISTOPIC_MODEL_PKL,"wb") as f: pickle.dump(cistopic_obj,f)
    nt = getattr(model, "n_topic", getattr(model, "n_topics", "NA"))
    print(f"[topics] auto-selected model n_topics={nt} -> {CISTOPIC_MODEL_PKL}")

# ===================================================================== STAGE regionsets
def stage_regionsets():
    """Binarize topics (otsu + top3k) + DARs per cell type -> region_sets/*.bed for pycistarget.
       SRC: https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html (region sets)."""
    from pycisTopic.topic_binarization import binarize_topics
    from pycisTopic.diff_features import (impute_accessibility, normalize_scores,
                                          find_highly_variable_features, find_diff_features)
    from pycisTopic.utils import region_names_to_coordinates
    import gc
    with open(CISTOPIC_MODEL_PKL,"rb") as f: cistopic_obj = pickle.load(f)
    # Topic binarization needs NO impute (uses topic-region dist) -> low memory, on ALL cells.
    otsu   = binarize_topics(cistopic_obj, method="otsu")
    top3k  = binarize_topics(cistopic_obj, method="ntop", ntop=3000)
    # ---- OOM-FIX (2nd OOM): impute_accessibility on ALL 48113 cells x 407852 regions is float32 ~78G but
    # peaked to 127.9G -> OOM. Stratified-subsample cells (proportional by cell_type; fibroblast well kept)
    # for impute + DARs. find_diff_features intersects cell_data with imp.cell_names (verified src lines 54-57),
    # so a cell-subset imp is valid. Methods-documented efficiency step; DARs RETAINED on representative cells.
    rng = np.random.default_rng(SEED)
    cd  = cistopic_obj.cell_data
    TARGET_CELLS = 10000                                  # 3rd-OOM was in HVF (95G+33G->128G); 10000 cells -> impute+norm+HVF peak ~79G < RAM; fibroblast still ~1900
    frac = min(1.0, TARGET_CELLS/len(cd))
    sel = []
    for ct, grp in cd.groupby("cell_type", observed=True):
        n = min(len(grp), max(1, int(round(len(grp)*frac))))
        sel += list(rng.choice(grp.index.values, size=n, replace=False))
    print(f"[regionsets] impute on stratified subsample {len(sel)}/{len(cd)} cells; per-type:",
          cd.loc[sel,"cell_type"].value_counts().to_dict(), flush=True)
    imp    = impute_accessibility(cistopic_obj, selected_cells=sel, selected_regions=None, scale_factor=10**6)
    norm   = normalize_scores(imp, scale_factor=10**4)
    varreg = find_highly_variable_features(norm, min_disp=0.05, min_mean=0.0125, max_mean=3,
                                           max_disp=np.inf, n_bins=20, n_top_features=None, plot=False)
    del norm; gc.collect()                                # free the normalized copy before DAR calling
    markers = find_diff_features(cistopic_obj, imp, variable="cell_type", var_features=varreg,
                                 contrasts=None, adjpval_thr=0.05, log2fc_thr=np.log2(1.5),
                                 n_cpu=N_CPU_COMPUTE, split_pattern=SEP)   # cell_data auto-intersected to imp cells
    def dump(d, sub):
        out=f"{RS_DIR}/{sub}"; os.makedirs(out, exist_ok=True)
        for k in d:
            idx=d[k].index if hasattr(d[k],"index") else d[k]
            if len(idx)==0: continue
            region_names_to_coordinates(idx).sort_values(["Chromosome","Start","End"]).to_csv(
                f"{out}/{str(k).replace(' ','_')}.bed", sep="\t", header=False, index=False)
    dump(otsu,"Topics_otsu"); dump(top3k,"Topics_top_3k"); dump(markers,"DARs_cell_type")
    print("[regionsets] wrote Topics_otsu / Topics_top_3k / DARs_cell_type beds")

# ====================================================================== STAGE scenicplus
def stage_scenicplus():
    """SCENIC+ snakemake: init the DEFAULT config (all tutorial-default params kept verbatim),
       SURGICALLY patch only the 6 input paths + temp_dir + n_cpu + biomart mirror (text replace,
       NOT a yaml round-trip -> preserves `1_000`, the lambda string, everything else exactly),
       then run `snakemake --cores 6`.  Outputs land under Snakemake/ (relative), collected by extract.
       SRC: https://scenicplus.readthedocs.io/en/latest/human_cerebellum.html (Running SCENIC+).
       NOTE: default bc_transform_func == identity `lambda x: f'{x}'`; our GEX obs_names already
       equal cisTopic cell_names (barcode___sample), so identity is correct — left untouched."""
    pipe = f"{BASE}/scplus_pipeline"
    cfg_path = f"{pipe}/Snakemake/config/config.yaml"
    if not os.path.exists(cfg_path):
        subprocess.run(["scenicplus","init_snakemake","--out_dir",pipe], check=True)
    txt = open(cfg_path).read()
    repl = {
        'cisTopic_obj_fname: ""'        : f'cisTopic_obj_fname: "{CISTOPIC_MODEL_PKL}"',
        'GEX_anndata_fname: ""'         : f'GEX_anndata_fname: "{ADATA_ANNOT}"',
        'region_set_folder: ""'         : f'region_set_folder: "{RS_DIR}"',
        'ctx_db_fname: ""'              : f'ctx_db_fname: "{RANKINGS_DB}"',
        'dem_db_fname: ""'              : f'dem_db_fname: "{SCORES_DB}"',
        'path_to_motif_annotations: ""' : f'path_to_motif_annotations: "{MOTIF2TF}"',
        'temp_dir: ""'                  : f'temp_dir: "{TMP_DIR}"',
        'n_cpu: 40'                     : f'n_cpu: {N_CPU_COMPUTE}',   # OOM-safe compute -> 8 cores (coexist w/ base-GRN)
        'biomart_host: "http://www.ensembl.org"' : 'biomart_host: "http://asia.ensembl.org"',
    }
    for k,v in repl.items():
        assert txt.count(k) == 1, f"config patch anchor not unique/found: {k!r} (count={txt.count(k)})"
        txt = txt.replace(k, v)
    open(cfg_path,"w").write(txt)
    print(f"[scenicplus] patched {cfg_path}")
    print(f"[scenicplus] launch (LONG, nohup): cd {pipe}/Snakemake && snakemake --cores {N_CPU}")

# ========================================================================= STAGE extract
def stage_extract():
    """Load final scplusmdata.h5mu -> dump eRegulon metadata, flag HES1/Notch-axis eRegulons,
       compute mean AUCell per cell type -> report which cell type(s) each Notch eRegulon is active in.
       Written defensively (SCENIC+ 1.0a2 MuData layout probed at runtime, keys printed).
       NOTE: finalized against the ACTUAL mudata after the run — this is a FAST post-hoc step."""
    import mudata
    mm = glob.glob(f"{BASE}/scplus_pipeline/**/scplusmdata.h5mu", recursive=True) + \
         glob.glob(f"{SCPLUS_DIR}/scplusmdata.h5mu")
    assert mm, "scplusmdata.h5mu not found — run the scenicplus stage first"
    path = mm[0]; mdata = mudata.read(path)
    print(f"[extract] loaded {path}\n  mods={list(mdata.mod)}\n  uns_keys={list(mdata.uns.keys())}")
    def tf_of(name):                       # eRegulon/column name -> TF symbol (prefix before _ / ( / +)
        return re.split(r"[_(+\-]", str(name))[0]
    # ---- eRegulon metadata (TF-region-gene links) ----
    notch_hit = {}
    for key in list(mdata.uns.keys()):
        if "regulon" in key.lower() and "metadata" in key.lower():
            meta = pd.DataFrame(mdata.uns[key])
            meta.to_csv(f"{SCPLUS_DIR}/eRegulon_metadata__{key}.csv", index=False)
            tfcol = "TF" if "TF" in meta.columns else ("Gene_signature_name" if "Gene_signature_name" in meta.columns else None)
            if tfcol:
                nm = meta[meta[tfcol].astype(str).isin(NOTCH_TFS)]
                if len(nm):
                    nm.to_csv(f"{SCPLUS_DIR}/NOTCH_eRegulon_metadata__{key}.csv", index=False)
                    notch_hit[key] = sorted(nm[tfcol].unique().tolist())
    # ---- AUCell per cell -> mean per cell_type ----
    ct = None
    for src in (mdata.obs, ):
        for c in src.columns:
            if c.split(":")[-1] == "cell_type": ct = src[c]; break
    reports = []
    for mod in mdata.mod:
        m = mdata[mod]
        try: df = m.to_df()
        except Exception: continue
        cc = ct
        if cc is None:
            for c in m.obs.columns:
                if c.split(":")[-1] == "cell_type": cc = m.obs[c]; break
        if cc is None: continue
        notch_cols = [c for c in df.columns if tf_of(c) in NOTCH_TFS]
        mean_by_ct = df.groupby(cc.values).mean()
        safe = re.sub(r"[^A-Za-z0-9]+","_",mod)
        mean_by_ct.to_csv(f"{SCPLUS_DIR}/AUCell_mean_by_celltype__{safe}.csv")
        if notch_cols:
            mean_by_ct[notch_cols].to_csv(f"{SCPLUS_DIR}/NOTCH_AUCell_by_celltype__{safe}.csv")
            top = {c: mean_by_ct[c].idxmax() for c in notch_cols}   # cell type with max activity
            reports.append((mod, top))
    print("[extract] Notch-axis eRegulons in metadata:", notch_hit if notch_hit else "NONE")
    print("[extract] Notch eRegulon AUCell -> most-active cell type per eRegulon:")
    for mod, top in reports: print(f"  [{mod}] {top}")
    print(f"[extract] CSVs (eRegulon_metadata / NOTCH_* / AUCell_mean_by_celltype) under {SCPLUS_DIR}")

# =============================================================================== main
STAGES = {"rna":stage_rna,"peaks":stage_peaks,"cistopic":stage_cistopic,"topics":stage_topics,
          "regionsets":stage_regionsets,"scenicplus":stage_scenicplus,"extract":stage_extract}
if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv)>1 else ""
    if st not in STAGES:
        print("usage: run_scenicplus_pipeline.py [%s]" % "|".join(STAGES)); sys.exit(1)
    print(f"=== STAGE {st} START ==="); STAGES[st](); print(f"=== STAGE {st} DONE ===")
