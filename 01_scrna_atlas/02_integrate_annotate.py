# -*- coding: utf-8 -*-
"""
================================================================================
[SERVER, npc env (scanpy/decoupler/harmonypy)] -- run on the aliyun 16-32 vCPU /
128 GB instance. Stage-by-stage, checkpointed, idempotent.

STEP 02 of the POWERED / max-data SSc virtual-perturbation pipeline
--------------------------------------------------------------------------------
MAX-DATA ATLAS INTEGRATION + CLEAN ANNOTATION.

Loads the multiple fibroblast + immune-co-context scRNA cohorts (config list of
h5ad / 10x paths), does per-sample:
  (a) ambient-RNA decontamination          -> DecontX (R/celda) OR scAR (python)
  (b) doublet removal                       -> Scrublet (python) OR scDblFinder (R)
  (c) data-driven QC                        -> MAD-based mito / nGene thresholds
then INTEGRATES across (cohort x sample) batch with Harmony (harmonypy, CPU;
scVI GPU path optional & flagged), and ANNOTATES CELL TYPES FROM SCRATCH via a
rigorous multi-reference CONSENSUS (NO prior labels are trusted -- see below),
including a myofibroblast SIGNATURE score, and keeps the immune compartments
(myeloid / T / B / mast) for the immune layer.

OUTPUTS (to /data/ssc/powered/):
  full_integrated_annotated.h5ad   -- ALL cells, integrated + annotated (from scratch)
  fibroblast_integrated.h5ad       -- fibroblast subset (input to CellOracle/regulon)
  per_sample_qc_report.csv         -- per-sample QC / decont / doublet counts
  --- annotation AUDIT bundle (for human confirmation, STEP3) ---
  annotation_marker_scores.csv     -- cluster x marker-panel z-score matrix
  annotation_cluster_DE.csv        -- per-cluster top rank_genes_groups DE genes
  annotation_celltypist.csv        -- CellTypist per-cluster majority label(s)+conf
  annotation_consensus.csv         -- per-cluster consensus table + needs_review flag
  --- fibroblast subtype AUDIT bundle (STEP4) ---
  fibro_subtype_marker_scores.csv  -- fibro-cluster x subtype z-score matrix
  fibro_subtype_DE.csv             -- per-fibro-cluster top DE genes
  fibro_subtype_consensus.csv      -- fibro subtype consensus + needs_review flag
Each STEP writes an intermediate checkpoint and is SKIPPED if its output exists.

*** ANNOTATION IS REDONE FROM SCRATCH -- NO PRIOR LABELS TRUSTED ***
The user is skeptical of ALL pre-existing annotations, INCLUDING the project's own
Tabib GSE138669 labels (server_results/tabib/tabib_annotated.h5ad). This script
therefore ignores every incoming *_annotated.h5ad label and re-derives cell types
on the FRESHLY integrated object, like a careful human annotator:
  (1) Leiden clustering with a RESOLUTION SWEEP + a stability pick (not one res).
  (2) Each cluster is scored against MULTIPLE published marker panels (below), the
      score is z-scored across clusters, and argmax gives a MARKER call.
  (3) An independent REFERENCE-MAPPING call is produced by CellTypist
      (Immune_All_Low + Adult_Human_Skin models, majority_voting) -- verified real
      API (celltypist.models.download_models / Model.load / annotate).
  (4) The FINAL label = CONSENSUS of the marker argmax AND the reference-mapping
      call, mapped to a shared vocabulary; disagreements raise a needs_review flag.
  (5) Full AUDIT CSVs (DE genes, score matrix, CellTypist labels, consensus table)
      are written so a human can confirm every call. A per-cluster verdict is printed.

MARKER PANELS -- cross-compared, each with its literature source (verbatim genes,
verified against the cited papers; see MARKERS / FIB_SUBTYPES blocks below):
  * Tabib 2018  (PMID 29080679, J Invest Dermatol): dermal fibroblast axes
    SFRP2/DPP4 (secretory-reticular) vs FMO1/LSP1 (secretory-papillary).
  * Tabib 2021  (PMC8289865, Nat Commun): SFRP2hi progenitor -> myofibroblast
    trajectory; PRSS23/THBS1; SFRP4+ pro-fibrotic expansion in SSc.
  * Gur 2022    (PMC7612792, Cell): LGR5+ ScAF; COMP/PRSS23/WIF1 functional subsets.
  * ARD 2025    (PMID 40410053): fibroblast-macrophage niche (context markers).
  * canonical immune/stromal panels (T/myeloid/mast/NK/B/plasma/EC/pericyte-SMC/
    keratinocyte/melanocyte) from standard skin scRNA references.

WHY THIS DESIGN (reuse, not from scratch):
  * QC (DecontX/Scrublet/MAD) + Harmony integration are KEPT verbatim from the prior
    version of this script (STEP1/STEP2 unchanged).
  * The annotation stage is REBUILT: marker-score assignment now cross-compares
    multiple panels AND a reference-mapping model, instead of a single argmax.
  * MYO program + per-cell fibroblast-purity filter == rerun_clean/_qc_common.py.
  * scVI-optional integration path                  == library module 506
    (08_singlecell_spatial_trajectory/506_scvi_scanvi_integration): honest-baseline
    note -- learned embeddings do NOT automatically beat Harmony (Genome Biol 2025);
    Harmony is the CPU default here, scVI is opt-in and compared, not assumed better.
  * CellTypist reference-mapping           == verified real API (Teichlab); models
    Immune_All_Low.pkl + Adult_Human_Skin.pkl (per-model try/except, never hard-fails).

HONEST-FRAMING NOTES (baked into comments, per the iron laws):
  * Immune compartments are kept ONLY for association / co-context downstream;
    this script does NOT and MUST NOT assert any immune -> HES1 direction.
  * TREM2+ macrophages are PROTECTIVE in SSc skin -- annotation is descriptive only.
  * Nothing here is causal: cell-type calls are marker-correlation, not proof.

REPRODUCIBILITY: fixed SEED, sparse/backed reads, gzip-compressed writes, every
heavy step checkpointed + idempotent, all paths in the CONFIG block below.
Run:
  conda run -n npc python /data/ssc/server_powered_2026/02_integrate_annotate.py
================================================================================
"""
import os, sys, gc, json, glob, re, time, subprocess, tempfile, traceback
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp

sc.settings.verbosity = 1
SEED = 0
np.random.seed(SEED)

# ==============================================================================
# CONFIG  (edit paths here only -- no hardcoding below this block)
# ==============================================================================
RAW   = "/data/ssc/raw"          # all downloaded cohorts live under here
OUT   = "/data/ssc/powered"      # step-02 outputs
CKPT  = os.path.join(OUT, "_ckpt")   # per-step intermediate checkpoints
for d in (OUT, CKPT):
    os.makedirs(d, exist_ok=True)

# ---- ambient / doublet backend switches -------------------------------------
# DecontX lives in the R package 'celda' (Bioconductor); there is NO python port.
# Two honest choices, controlled here:
#   AMBIENT_METHOD = "decontx"  -> call R via Rscript (celda::decontX); most-cited,
#                                  matches what SSc atlas papers use, but needs the
#                                  'decont_env' R env with celda installed.
#   AMBIENT_METHOD = "scar"     -> pure-python (scar, scvi-tools ecosystem); GPU-
#                                  friendly, no R dependency.
#   AMBIENT_METHOD = "none"     -> skip ambient correction (fallback if neither env
#                                  is ready; QC + doublet still run). LOG-only.
AMBIENT_METHOD  = "decontx"      # "decontx" | "scar" | "none"
RSCRIPT_BIN     = os.environ.get("RSCRIPT", "/data/ssc/miniconda3/envs/r44/bin/Rscript")  # r44 有 celda/DropletUtils (2026-07-07 修)

# scDblFinder is R (Bioconductor); Scrublet is python. Python default = no extra env.
DOUBLET_METHOD  = "scrublet"     # "scrublet" | "scdblfinder" | "none"

# ---- integration backend -----------------------------------------------------
# Harmony (harmonypy) is the CPU default and always runs. scVI is an OPTIONAL GPU
# path (library module 506); set USE_SCVI=True only when a GPU + scvi-tools are
# present -- it is COMPARED against Harmony (dual metric), never assumed better.
USE_SCVI      = False            # optional GPU path; flagged off by default
SCVI_EPOCHS   = 200

# ---- QC (MAD-based, data-driven) --------------------------------------------
N_MADS        = 5.0              # outlier = median +/- N_MADS*MAD  (log1p scale)
MITO_N_MADS   = 3.0              # TIGHTER MAD for mito only (A3): sc-best-practices
                                 # tutorial uses 3 MADs for pct_mt (vs 5 for nGene/nUMI).
MITO_HARD_MAX = 20.0             # absolute ceiling on pct_counts_mt (%). KEPT AT 20 (NOT the
                                 # tutorial's 8%): our own M3 sensitivity analysis showed the 8%
                                 # cap is unsuitable for the MARS-seq Gur discovery cohort
                                 # (median mito ~16%; an 8% cut drops ~95% of cells) and the
                                 # leads are NOT high-mito artifacts -> tighten the MAD (3) but
                                 # keep the justified 20% ceiling.
MIN_GENES     = 200             # floor before per-cell metrics (matches project)
MIN_CELLS     = 10               # gene filter (matches project)

# emptyDrops (DropletUtils) cell-calling for RAW 10x matrices (unfiltered barcodes).
# Cohorts whose input is a CellRanger raw_feature_bc_matrix (full ~737k/1.8M barcode
# whitelist, e.g. Tabib GSE138669, GSE181957, GSE292979) get statistical empty-droplet
# removal BEFORE ambient/doublet/QC -- min_genes alone is not a valid substitute
# (Lun 2019 DropletUtils; sc-best-practices / Heumos 2023 / Luecken 2019 all require it).
EMPTYDROPS_FDR   = 0.001   # A13: DELIBERATE STRICT cut (the DropletUtils vignette uses
                           # 0.01; 0.001 trades a few real cells for near-zero empty-droplet
                           # contamination -- a conservative atlas-quality choice, not standard).
EMPTYDROPS_LOWER = 100     # total-count threshold below which barcodes seed the ambient pool

# ---- integration / clustering ------------------------------------------------
N_HVG         = 3000             # HVG for the shared embedding (atlas-scale)
N_PCS         = 40
# Annotation is done on a RESOLUTION SWEEP + stability pick (human-like), not one res.
# We cluster at each res in LEIDEN_RES_SWEEP and choose the one whose clustering is
# most STABLE across neighbouring resolutions (max mean ARI to its two neighbours).
LEIDEN_RES_SWEEP = [0.4, 0.6, 0.8, 1.0, 1.2, 1.5]   # global compartment sweep
LEIDEN_RES       = 1.0          # fallback/default if the sweep can't be scored
FIB_SUBCL_RES    = 0.6          # fibroblast-only subtype resolution

# ---- reference-mapping (CellTypist) -----------------------------------------
# Verified REAL API against the actual Teichlab celltypist source + readthedocs
# tutorial (checked 2026-07-02, celltypist 1.7.x):
#   from celltypist import models
#   models.download_models(model=[...])            # subset download by .pkl name
#   celltypist.annotate(adata, model='X.pkl', majority_voting=True,
#                       over_clustering=<obs col or array>) -> AnnotationResult
#     * annotate() signature (source annotate.py) has exactly these kwargs:
#       filename, model, transpose_input, gene_file, cell_file, mode, p_thres,
#       majority_voting, over_clustering, use_GPU, min_prop.
#     * We pass over_clustering="leiden" so majority voting is refined WITHIN our
#       own Leiden clusters (documented: over_clustering may be an obs column key),
#       which makes the per-cluster reference label well defined and comparable to
#       the marker call. (Default over_clustering is an internal res=5 clustering.)
#   result.predicted_labels -> DataFrame with columns 'predicted_labels',
#       'over_clustering', 'majority_voting'  (verified in classifier.py).
#   result.to_adata() -> writes those + 'conf_score' into .obs (verified to_adata()).
# INPUT REQUIREMENT (tutorial, verbatim): adata.X must be log1p-normalised to 1e4
#   counts/cell. Our integrated .X is exactly that (STEP2: normalize_total 1e4 ->
#   log1p; scaling is done only on a throwaway HVG copy). We still pass a fresh
#   AnnData built from adata.raw (full gene space, log1p-1e4) to be unambiguous.
# MODELS (both VERIFIED to exist in celltypist models.json, 2026-07-02):
#   Immune_All_Low.pkl   -- Pan-immune, 98 subtypes (Science 2022, abl5197). Default.
#   Adult_Human_Skin.pkl -- healthy adult skin, 34 types (Human_Skin_Reynolds/v1,
#                           Science 2021, aba6500) -> gives stromal/epithelial labels
#                           that Immune_All_Low cannot (fibroblast/KC/EC/melanocyte).
# Each model is loaded in its own try/except: a missing/renamed/undownloadable model
# is LOGGED and skipped, never crashes annotation (marker consensus still runs). Set
# USE_CELLTYPIST=False to disable reference-mapping entirely (marker-only; the audit
# consensus then flags every cluster ref_label='NA', method='marker_only').
USE_CELLTYPIST     = True
# Immune_All_Low (coarse immune) + Adult_Human_Skin (structural/stromal) + Immune_All_High
# (FINE immune subtypes, 98 types; Dominguez Conde 2022 Science abl5197) -- the High model
# is the literature-grounded reference that resolves the immune fine-type clusters (NK-vs-
# cytotoxic-T, B, pDC) that a coarse panel leaves marker-unsupported (needs_review).
CELLTYPIST_MODELS  = ["Immune_All_Low.pkl", "Immune_All_High.pkl", "Adult_Human_Skin.pkl"]
CELLTYPIST_CONF_MIN = 0.5   # N2: per-cell max class-probability below this = low-confidence
                            # (flagged 'Unknown' even inside a confident cluster; documented cut)
# needs_review is raised when marker-call and reference-mapping-call disagree after
# mapping both into the shared vocabulary (SHARED_VOCAB below).

# ---- COHORT MANIFEST ---------------------------------------------------------
# Each entry describes ONE cohort. Loader picks the right reader by 'kind':
#   kind="h5ad"   -> single pre-built .h5ad (path = file)
#   kind="10x_h5" -> directory of *.h5 (one per sample; filtered_feature_bc_matrix)
#   kind="10x_mtx"-> directory of per-sample 10x triplets (matrix/barcodes/features)
#   kind="csv"    -> genes x cells .csv(.gz) counts (discovery-style)
#   kind="mars"   -> Gur-style per-plate genes x cells *.txt.gz tables
# 'compartment' just documents the intended content (fibro atlas vs immune co-context);
# ALL cells are kept through integration -- annotation decides fibro vs immune.
# 'sample_regex' extracts a per-sample id from filenames (for the batch key).
# kind/glob/sample_regex VERIFIED against the actual GEO deposits (peeked 2026-07-04);
# the earlier guesses (all "10x_mtx") were WRONG for 5/7 cohorts -- fixed below.
COHORTS = [
    # ---- fibroblast atlas ----
    dict(label="Gur_GSE195452",   kind="mars",             # per-plate genes x cells txt
         path=os.path.join(RAW, "GSE195452"),      compartment="fibro_atlas",
         glob="*.txt.gz", sample_regex=r"(GSM\d+)"),
    dict(label="Tabib_GSE138669", kind="10x_h5",            # per-sample raw_feature_bc_matrix.h5
         path=os.path.join(RAW, "GSE138669"),      compartment="fibro_atlas",
         glob="*.h5",     sample_regex=r"_(SC\d+)raw", raw_matrix=True),   # 737280 barcodes = RAW
    dict(label="GSE292979",       kind="10x_mtx_prefixed",  # flat GSM_*_matrix.mtx.gz+barcodes+features
         path=os.path.join(RAW, "GSE292979"),      compartment="fibro_atlas",
         glob=None,       sample_regex=r"(GSM\d+)", raw_matrix=True),      # 1.83M barcodes = RAW
    dict(label="GSE249279",       kind="10x_h5_merged",     # single merged 10x h5 (aggr); sample=barcode suffix
         path=os.path.join(RAW, "GSE249279", "GSE249279_scRNA-seq_merged_raw_counts.h5"),
         compartment="fibro_atlas", glob=None,     sample_regex=None),
    # ---- immune co-context (kept for the immune-association layer ONLY) ----
    dict(label="GSE181957",       kind="10x_h5",            # per-sample raw_feature_bc_matrix.h5 (like Tabib)
         path=os.path.join(RAW, "GSE181957"),      compartment="immune_cocontext",
         glob="*.h5",     sample_regex=r"_(SC\d+)raw", raw_matrix=True),   # 737280 barcodes = RAW
    dict(label="GSE210395",       kind="long_triplet",      # (gene,cell,count) long table; sample=barcode suffix
         path=os.path.join(RAW, "GSE210395"),      compartment="immune_cocontext",
         glob="*countMatrix.tsv.gz", sample_regex=None, sep="\t"),
    dict(label="GSE236111",       kind="csv",               # per-sample genes x cells counts.tsv.gz (+cellname map)
         path=os.path.join(RAW, "GSE236111"),      compartment="immune_cocontext",
         glob="*counts.tsv.gz", sample_regex=r"(GSM\d+)", sep="\t"),
]

# ==============================================================================
# GENE PANELS  (reused verbatim from the project; extended for subtypes)
# ==============================================================================
# myofibroblast program -- IDENTICAL to _qc_common.MYO / srv03_regulon_core.MYO
MYO = ["ACTA2", "TAGLN", "POSTN", "COL1A1", "COL1A2", "COMP", "COL11A1", "CTHRC1", "COL3A1", "FN1"]

# TOP-LEVEL COMPARTMENT MARKERS -- cross-compared multi-literature panels.
# Each panel = canonical lineage genes verified against the cited references so a
# reviewer can trace every marker. Prompt-mandated genes are all present:
#   T CD3D/CD8A ; myeloid LYZ/CD68/SPP1/TREM2/FCN1/LYVE1 ; mast KIT/TPSAB1 ;
#   NK NCAM1/NKG7 ; B/plasma MS4A1/MZB1 ; endothelial PECAM1/VWF ;
#   pericyte/SMC RGS5/ACTA2/MYH11 ; keratinocyte KRT1/KRT14 ; melanocyte MLANA.
# Refs: Tabib 2018 (PMID 29080679) & 2021 (PMC8289865) skin fibroblast; Gur 2022
# (PMC7612792) ScAF; ARD 2025 (PMID 40410053) fibro-macrophage niche; canonical
# skin/immune atlas panels (Reynolds Science 2021 skin; Monaco Cell Rep 2019 immune).
MARKERS = {
    "Fibroblast":   ["COL1A1", "COL1A2", "LUM", "DCN", "PDGFRA", "PDGFRB"],   # Tabib18/21, Gur22
    "Tcell":        ["CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD8A", "CD4"],   # +CD8A (prompt)
    "Myeloid":      ["LYZ", "CD68", "CD14", "AIF1", "ITGAX", "CD163",
                     "FCN1", "SPP1", "TREM2", "LYVE1"],                       # +FCN1/SPP1/TREM2/LYVE1 (prompt, ARD25 niche)
    "Endothelial":  ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Keratinocyte": ["KRT14", "KRT5", "KRT1", "KRT10", "KRT15"],
    "SmoothMuscle": ["MYH11", "DES", "TAGLN", "ACTA2", "CNN1"],               # RGS5/ACTA2/MYH11 (prompt)
    "Pericyte":     ["RGS5", "PDGFRB", "NOTCH3", "MCAM", "KCNJ8"],            # RGS5 (prompt)
    "Bcell":        ["MS4A1", "CD79A", "CD79B", "IGHM", "BANK1"],             # MS4A1 (prompt)
    "Plasma":       ["MZB1", "IGHG1", "JCHAIN", "XBP1", "DERL3"],             # MZB1 (prompt) -- split from B
    "Mast":         ["TPSAB1", "CPA3", "MS4A2", "KIT"],                       # KIT/TPSAB1 (prompt)
    "Melanocyte":   ["MLANA", "PMEL", "DCT", "TYRP1"],                        # MLANA (prompt)
    "NK":           ["NKG7", "GNLY", "KLRD1", "NCAM1", "KLRF1"],              # NCAM1/NKG7 (prompt)
    "Lymphatic":    ["PROX1", "LYVE1", "CCL21", "MMRN1"],                     # endothelial-lymphatic
    "pDC":          ["LILRA4", "CLEC4C", "IRF7", "GZMB"],
}

# SHARED VOCABULARY: maps the many CellTypist label strings (fine immune + skin
# models emit dozens of leaf labels) onto the coarse compartment vocabulary above,
# so the marker-call and the reference-mapping-call are compared in ONE space.
# Matching is substring/keyword based (case-insensitive) and is applied to the
# CellTypist majority label; anything unmatched -> "Other" (never silently forced).
SHARED_VOCAB = {
    "Fibroblast":   ["fibroblast", "stromal", "myofibro", "sfrp", "scaf"],
    "Tcell":        ["t cell", "t-cell", " t ", "cd4", "cd8", "treg", "tcm", "tem",
                     "th1", "th2", "th17", "mait", "gamma", "ilc"],
    "NK":           ["nk", "natural killer", "cytotoxic"],
    "Bcell":        ["b cell", "b-cell", "memory b", "naive b", "germinal"],
    "Plasma":       ["plasma", "plasmablast"],
    "Myeloid":      ["macrophage", "monocyte", "dendritic", "dc", "myeloid",
                     "langerhans", "kupffer", "microglia", "mono", "mph"],
    "Mast":         ["mast"],
    "pDC":          ["plasmacytoid", "pdc"],
    "Endothelial":  ["endothelial", "vascular", "vein", "artery", "capillary",
                     "ve1", "ve2", "ve3", "ve4"],   # Reynolds skin-atlas VE1-4 codes
    "Lymphatic":    ["lymphatic", "lec", "le1", "le2"],   # skin-atlas LE1/LE2 codes
    "Pericyte":     ["pericyte", "mural"],
    "SmoothMuscle": ["smooth muscle", "myocyte", "vsmc"],
    "Keratinocyte": ["keratinocyte", "keratin", "epiderm", "epithelial", "_kc",
                     "spinous", "basal", "granular", "corneocyte", "suprabasal"],
    "Melanocyte":   ["melanocyte", "melano"],
}

# ---- fibroblast subtype reference states -------------------------------------
# Tabib 2018 (J Invest Dermatol) SFRP2/DPP4 + FMO1/LSP1 papillary/reticular axis;
# Tabib 2021 (Nat Commun) SFRP4+/PRSS23+ pro-fibrotic expansion in SSc;
# Gur 2022 (Cell) LGR5+/SFRP2+/... functional fibroblast states + myofibroblast.
# These are marker-score REFERENCE panels for descriptive subtyping (correlation,
# not a causal claim). Myofibroblast is scored separately by MYO (above).
FIB_SUBTYPES = {
    "SFRP2_DPP4":     ["SFRP2", "DPP4", "COL1A1", "PI16"],       # Tabib18 papillary/universal
    "FMO1_LSP1":      ["FMO1", "LSP1", "APOD", "CD34"],          # Tabib18 reticular
    "SFRP4_proFib":   ["SFRP4", "PRSS23", "COMP", "TNC"],        # Tabib21 SSc pro-fibrotic
    "Myofibroblast":  ["ACTA2", "CTHRC1", "COMP", "POSTN", "COL11A1", "FN1"],  # the fibrotic program
    "LGR5_Gur":       ["LGR5", "SFRP2", "WIF1"],                 # Gur22 functional state
    "Inflammatory":   ["CCL19", "CXCL12", "CD74", "HLA-DRA"],    # inflammatory/immune-interacting
    "Adipogenic":     ["APOE", "PPARG", "CFD", "GPX3"],          # adipose-derived
}

# myofibroblast SIGNATURE (the requested keep) -- explicit list from the prompt
MYOFIB_SIG = ["ACTA2", "CTHRC1", "COMP", "POSTN", "COL11A1", "FN1"]

# ==============================================================================
# LITERATURE-GROUNDED CLEANING PANELS (added 2026-07-04). Every gene set is quoted
# from a cited source -- NOT hand-curated -- so the QC is defensible at top-journal
# review. Used by rbc_filter() / score_dissoc_cellcycle() below.
# ==============================================================================
# (i) DISSOCIATION / immediate-early artifact signature -- the genes van den Brink
#     et al. 2017 (Nat Methods, nmeth.4437, Fig 1a) report as INDUCED purely by warm
#     enzymatic tissue dissociation (AP-1 IEGs + HSPs + SOCS3/ZFP36). ★FOSB is IN this
#     list and is a project lead -> this score is the mandatory control that shows the
#     FOSB signal is not a dissociation artifact (used as a covariate in layer 07).
DISSOC_SIG = ["FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3", "EGR1", "CEBPB",
              "CEBPD", "HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "SOCS3", "ZFP36"]

# (ii) HEMOGLOBIN / erythrocyte genes -- flag & remove RBC-contaminated barcodes
#     (erythrocytes express almost exclusively Hb; MAD-mito does not catch them).
HB_GENES   = ["HBA1", "HBA2", "HBB", "HBD", "HBM", "HBQ1", "HBZ", "ALAS2"]
HB_MAX_PCT = 20.0   # documented CFG: drop cells with >20% counts from Hb genes (clear RBC)

# (iii) CELL-CYCLE gene sets -- the canonical Tirosh et al. 2016 (Science aad0501)
#     "regev_lab" S-phase and G2/M lists shipped with Seurat/scanpy. Score-only (we do
#     NOT regress: fibroblast proliferation may be real biology, not a confound).
CC_S_GENES = ["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
              "CDCA7","DTL","PRIM1","UHRF1","CENPU","HELLS","RFC2","RPA2","NASP",
              "RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2",
              "RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2",
              "USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"]
CC_G2M_GENES = ["HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
                "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","PIMREG","SMC4","CCNB2",
                "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1",
                "KIF20B","HJURP","CDCA3","JPT1","CDC20","TTK","CDC25C","KIF2C",
                "RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR",
                "AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3",
                "GAS2L3","CBX5","CENPA"]


# ==============================================================================
# helpers
# ==============================================================================
def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _ck(name):
    return os.path.join(CKPT, name)


def mad_bounds(x, n_mads=N_MADS, log1p=True):
    """Data-driven outlier bounds via median +/- n_mads*MAD (optionally on log1p).
    Returns (lo, hi) on the ORIGINAL scale. This is the standard scanpy-tutorial /
    sc-best-practices MAD rule -- replaces the fixed 300-7000 gene / 15% mito cuts
    from the project scripts with a per-dataset data-driven threshold."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]   # ★2026-07-08: drop NaN/inf (e.g. pct_mt=0/0 from empty MARS-seq wells)
    if x.size == 0:         #   BEFORE the median — np.median(NaN)=NaN poisoned the bound -> mito<=NaN kept 0 cells
        return 0.0, float("inf")   # no finite data -> non-restrictive bounds; caller's MIN_GENES floor + QC still apply
    v = np.log1p(x) if log1p else x
    med = np.median(v)
    mad = np.median(np.abs(v - med)) * 1.4826  # scale to ~sd
    lo, hi = med - n_mads * mad, med + n_mads * mad
    if log1p:
        lo, hi = np.expm1(lo), np.expm1(hi)
    return float(lo), float(hi)


# ------------------------------------------------------------------ LOADERS ----
def _sample_id(fname, regex, fallback):
    if regex:
        m = re.search(regex, os.path.basename(fname))
        if m:
            return m.group(1)
    return fallback


def _call_cells(x, c):
    """Per-sample cell-calling. RAW matrices (c['raw_matrix']=True, i.e. CellRanger
    raw_feature_bc_matrix with the full barcode whitelist) get statistical empty-droplet
    removal via emptyDrops; all others get the light min_genes floor (they were already
    cell-called upstream). Returns the cell-called AnnData."""
    if c.get("raw_matrix"):
        sc.pp.filter_cells(x, min_counts=1)          # drop all-zero barcodes pre-emptyDrops
        with tempfile.TemporaryDirectory() as tmp:
            x, _, _ = empty_drops_filter(x, tmp)
        return x
    sc.pp.filter_cells(x, min_genes=MIN_GENES)
    return x


def load_cohort(c):
    """Return a per-cohort AnnData (cells x genes, raw counts in .X) with
    obs['sample'] (batch) and obs['cohort']. Reuses the exact readers from
    tabib_process.py / srv04 / srv06. Light min_genes filter only (full QC later)."""
    label, kind, path = c["label"], c["kind"], c["path"]
    ads = []
    if kind == "h5ad":
        a = sc.read_h5ad(path)
        a.var_names_make_unique()
        if "sample" not in a.obs:
            a.obs["sample"] = label
        ads = [a]

    elif kind == "10x_h5":
        files = sorted(glob.glob(os.path.join(path, "**", c["glob"] or "*.h5"), recursive=True))
        for f in files:
            s = _sample_id(f, c["sample_regex"], os.path.basename(f).split("_")[0])
            x = sc.read_10x_h5(f); x.var_names_make_unique()
            x = _call_cells(x, c)
            x.obs["sample"] = s
            x.obs_names = [f"{s}_{b}" for b in x.obs_names]
            ads.append(x)

    elif kind == "10x_mtx":
        # per-sample 10x triplets under subdirs; else *.h5 fallback
        mtx = sorted(glob.glob(os.path.join(path, "**", "matrix.mtx*"), recursive=True))
        if mtx:
            for m in mtx:
                d = os.path.dirname(m); s = os.path.basename(d)
                x = sc.read_10x_mtx(d); x.var_names_make_unique()
                x = _call_cells(x, c)
                x.obs["sample"] = s
                x.obs_names = [f"{s}_{b}" for b in x.obs_names]
                ads.append(x)
        else:
            for f in sorted(glob.glob(os.path.join(path, "**", "*.h5"), recursive=True)):
                s = _sample_id(f, c["sample_regex"], os.path.basename(f).split("_")[0])
                x = sc.read_10x_h5(f); x.var_names_make_unique()
                x = _call_cells(x, c)
                x.obs["sample"] = s
                x.obs_names = [f"{s}_{b}" for b in x.obs_names]
                ads.append(x)

    elif kind == "csv":
        files = sorted(glob.glob(os.path.join(path, c["glob"] or "*.csv*")))
        for f in files:
            # genes x cells (project convention); sep=c['sep'] supports .tsv (\t) too
            df = pd.read_csv(f, index_col=0, sep=c.get("sep", ","))
            s = _sample_id(f, c["sample_regex"], os.path.basename(f).split("_")[0])
            x = ad.AnnData(sp.csr_matrix(df.T.values.astype("float32")))
            x.obs_names = [f"{s}_{col}" for col in df.columns]; x.var_names = list(df.index)
            x.var_names_make_unique(); x.obs["sample"] = s
            x = _call_cells(x, c)
            ads.append(x)

    elif kind == "mars":  # Gur-style MARS-seq per-plate tables (srv06_gur.py)
        files = sorted(glob.glob(os.path.join(path, "**", c["glob"] or "*.txt.gz"), recursive=True))
        files = [f for f in files if "metadata" not in os.path.basename(f).lower()]  # 排除 GSE195452 Cell_metadata(非 plate; 2026-07-07 修)
        ref_genes, parts, cells, samps = None, [], [], []
        for i, f in enumerate(files):
            df = pd.read_csv(f, sep="\t", index_col=0)          # genes x cells
            if ref_genes is None:
                ref_genes = df.index
            elif not df.index.equals(ref_genes):
                df = df.reindex(ref_genes, fill_value=0)
            parts.append(sp.csr_matrix(df.values.T.astype("float32")))   # cells x genes
            s = _sample_id(f, c["sample_regex"], f"plate{i}")
            cells += [f"{s}_{col}" for col in df.columns]; samps += [s] * df.shape[1]
            if i % 25 == 0:
                log(f"    {label}: plate {i}/{len(files)}")
        X = sp.vstack(parts).tocsr()
        a = ad.AnnData(X); a.obs_names = cells; a.var_names = list(ref_genes)
        a.var_names_make_unique(); a.obs["sample"] = samps
        ads = [a]

    elif kind == "10x_mtx_prefixed":
        # flat dir of per-sample 10x triplets sharing a file PREFIX (GEO RAW.tar
        # convention), e.g. GSM..._SAMPLE_matrix.mtx.gz + _barcodes.tsv.gz +
        # _features.tsv.gz. Group by the matrix-file prefix; read each triplet.
        import gzip
        from scipy.io import mmread
        mtxs = sorted(glob.glob(os.path.join(path, "**", "*matrix.mtx*"), recursive=True))
        if not mtxs:
            raise FileNotFoundError(f"{label}: no *matrix.mtx* under {path}")
        for mp in mtxs:
            d = os.path.dirname(mp); base = os.path.basename(mp)
            prefix = base[:base.index("matrix.mtx")]      # '..._'
            s = _sample_id(prefix, c["sample_regex"], prefix.rstrip("_").split("_")[-1])

            def _sib(stem):
                g = glob.glob(os.path.join(d, prefix + stem + "*"))
                return g[0] if g else None
            bcp = _sib("barcodes.tsv")
            ftp = _sib("features.tsv") or _sib("genes.tsv")
            if not (bcp and ftp):
                raise FileNotFoundError(f"{label}: barcodes/features missing for prefix {prefix} in {d}")
            op = (lambda p: gzip.open(p, "rb")) if mp.endswith(".gz") else (lambda p: open(p, "rb"))
            M = mmread(op(mp)).T.tocsr().astype("float32")    # mtx = genes x cells -> cells x genes
            barcodes = pd.read_csv(bcp, header=None, sep="\t")[0].astype(str).tolist()
            feats = pd.read_csv(ftp, header=None, sep="\t")
            genes = (feats[1] if feats.shape[1] > 1 else feats[0]).astype(str).tolist()
            x = ad.AnnData(M); x.var_names = genes
            x.obs_names = [f"{s}_{b}" for b in barcodes]
            x.var_names_make_unique(); x.obs["sample"] = s
            x = _call_cells(x, c)
            ads.append(x)

    elif kind == "10x_h5_merged":
        # a SINGLE merged 10x-format .h5 (cellranger aggr): per-cell sample = the
        # barcode suffix '-N' (aggr library index). Split gives the per-sample batch
        # key needed for per-sample decontamination / QC downstream.
        x = sc.read_10x_h5(path); x.var_names_make_unique()
        suf = [bc.rsplit("-", 1)[1] if "-" in bc else "1" for bc in x.obs_names]
        x.obs["sample"] = [f"{label}_agg{si}" for si in suf]
        x.obs_names = [f"{label}_{b}" for b in x.obs_names]
        sc.pp.filter_cells(x, min_genes=MIN_GENES)
        ads = [x]

    elif kind == "long_triplet":
        # single LONG-format table of (gene, cell, count) triples (col0=gene,
        # col1=cell, col2=count). Pivot to a sparse cells x genes matrix. Per-cell
        # sample = the 10x barcode suffix after the first '-' (e.g. '-SI-GA-B5').
        f = sorted(glob.glob(os.path.join(path, c["glob"] or "*.tsv*")))[0]
        df = pd.read_csv(f, sep=c.get("sep", "\t"))
        gcol, ccol, ncol = df.columns[0], df.columns[1], df.columns[2]
        gcat = pd.Categorical(df[gcol].astype(str))
        ccat = pd.Categorical(df[ccol].astype(str))
        M = sp.coo_matrix((df[ncol].values.astype("float32"), (ccat.codes, gcat.codes)),
                          shape=(len(ccat.categories), len(gcat.categories))).tocsr()
        x = ad.AnnData(M); x.var_names = list(gcat.categories)
        x.obs["_bc"] = list(ccat.categories)
        x.obs["sample"] = [bc.split("-", 1)[1] if "-" in bc else label for bc in x.obs["_bc"]]
        x.obs_names = [f"{s}_{bc}" for s, bc in zip(x.obs["sample"], x.obs["_bc"])]
        x.obs.drop(columns=["_bc"], inplace=True)
        x.var_names_make_unique()
        sc.pp.filter_cells(x, min_genes=MIN_GENES)
        ads = [x]

    else:
        raise ValueError(f"unknown kind={kind} for {label}")

    if not ads:
        raise FileNotFoundError(f"{label}: no input files found under {path}")
    A = ad.concat(ads, join="outer", index_unique=None)
    A.var_names_make_unique()
    A.obs["cohort"] = label
    A.obs["compartment_src"] = c["compartment"]
    # densify NaN from outer-join to 0 (missing gene in a sample = 0 count)
    if sp.issparse(A.X):
        A.X.data = np.nan_to_num(A.X.data)
    A.X = A.X.astype("float32")
    return A


# ---------------------------------------------------- AMBIENT DECONTAMINATION --
def _quick_z(adata_sample):
    """Per-cell cluster labels handed to DecontX so it SKIPS its own internal UMAP+dbscan
    cell-type estimation -- that UMAP (uwot) CRASHES when a sample is smaller than its
    n_neighbors ('n_neighbors must be smaller than the dataset size'). z is a documented,
    supported decontX input. Robust BY CONSTRUCTION (no try/except): tiny samples get a
    deterministic 2-way total-count split; larger ones a quick Leiden, guaranteed >=2
    groups (decontX needs >=2)."""
    n = adata_sample.n_obs
    tot = np.asarray(adata_sample.X.sum(1)).ravel()
    if n < 30:
        return (tot > np.median(tot)).astype(int)
    a = adata_sample.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=min(2000, a.n_vars))
    a = a[:, a.var.highly_variable].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=max(2, min(30, a.n_obs - 1, a.n_vars - 1)), random_state=SEED)
    sc.pp.neighbors(a, n_neighbors=min(15, a.n_obs - 1), random_state=SEED)
    sc.tl.leiden(a, resolution=1.0, random_state=SEED, flavor="igraph",
                 n_iterations=2, directed=False)
    z = a.obs["leiden"].astype(int).values
    if len(np.unique(z)) < 2:
        z = (tot > np.median(tot)).astype(int)
    return z


def decontx_R(counts_csr, z, tmp):
    """Call celda::decontX (R) on a raw-count matrix; return the decontaminated
    counts as CSR. DecontX has NO python port -- this shells out to Rscript in an
    R env with 'celda' installed. z = per-cell cluster labels (see _quick_z): passing
    z skips decontX's dbscan/kmeans cell-type ESTIMATION, but celda 1.26.0 STILL runs
    scater::calculateUMAP unconditionally, which crashes on tiny samples (<15 cells:
    n_neighbors >= n_cells). So we guard tiny samples below. Otherwise NO raw fallback:
    if Rscript/celda errors, this RAISES (set AMBIENT_METHOD='none' to skip honestly)."""
    # ★2026-07-07 修: <20 细胞样本 decontX 内部 UMAP 必崩且去污无意义 -> 返回原始 counts, 不 raise
    if counts_csr.shape[0] < 20:
        return counts_csr.astype("float32"), "none_tiny_sample"
    mtx_in  = os.path.join(tmp, "raw.mtx")
    mtx_out = os.path.join(tmp, "decont.mtx")
    z_in    = os.path.join(tmp, "z.txt")
    from scipy.io import mmwrite, mmread
    mmwrite(mtx_in, counts_csr.T.tocoo())   # R/celda expects genes x cells
    np.savetxt(z_in, np.asarray(z).astype(int), fmt="%d")
    # embed paths with forward slashes: Windows backslashes in an R string become
    # invalid escapes ('\\U...' -> R error). R accepts '/' on Windows and Linux.
    mtx_in_r, mtx_out_r = mtx_in.replace("\\", "/"), mtx_out.replace("\\", "/")
    z_in_r = z_in.replace("\\", "/")
    rscript = f"""
suppressMessages({{library(Matrix); library(celda)}})
set.seed(0)
m <- readMM("{mtx_in_r}")                    # genes x cells
z <- scan("{z_in_r}", quiet=TRUE)
res <- decontX(as(m, "CsparseMatrix"), z=z, seed=0)
writeMM(as(res$decontXcounts, "CsparseMatrix"), "{mtx_out_r}")
cat("DECONTX_OK\\n")
"""
    rf = os.path.join(tmp, "decontx.R")
    with open(rf, "w") as fh:
        fh.write(rscript)
    r = subprocess.run([RSCRIPT_BIN, rf], capture_output=True, text=True, timeout=7200)
    if "DECONTX_OK" in r.stdout and os.path.exists(mtx_out):
        dec = mmread(mtx_out).T.tocsr()   # back to cells x genes
        return dec.astype("float32"), "decontx"
    raise RuntimeError(
        f"decontX did NOT run (celda/Rscript env broken -- fix it, do not skip).\n"
        f"  RSCRIPT_BIN={RSCRIPT_BIN}\n  stdout tail: {r.stdout[-300:]}\n"
        f"  stderr tail: {r.stderr[-600:]}")


def empty_drops_filter(adata_sample, tmp, log=log):
    """Statistical empty-droplet removal on ONE raw sample via DropletUtils::emptyDrops
    (Lun et al. 2019). Input = a RAW barcode matrix (unfiltered; only all-zero barcodes
    pre-dropped). Barcodes with total counts <= EMPTYDROPS_LOWER seed the ambient profile;
    the rest are tested against it. Keeps barcodes with FDR <= EMPTYDROPS_FDR (real cells).
    Returns (filtered_adata, n_in, n_cells). NO fallback: if DropletUtils is missing or
    emptyDrops errors, RAISES (min_genes is NOT a valid substitute for a raw matrix)."""
    from scipy.io import mmwrite
    n_in = adata_sample.n_obs
    X = adata_sample.X if sp.issparse(adata_sample.X) else sp.csr_matrix(adata_sample.X)
    mtx_in = os.path.join(tmp, "raw_ed.mtx")
    mmwrite(mtx_in, X.T.tocoo())   # genes x cells
    fdr_out = os.path.join(tmp, "ed_fdr.txt")
    # forward slashes: Windows backslashes ('\\U...' in an R string) are invalid escapes
    mtx_in_r, fdr_out_r = mtx_in.replace("\\", "/"), fdr_out.replace("\\", "/")
    # A4: also capture res$Limited. A barcode is 'Limited' when its Monte-Carlo p-value
    # is capped by the number of permutations (niters); if any such barcode is still
    # NON-significant at our FDR cut, the permutation count was insufficient and niters
    # should be raised. We WARN (do NOT auto-loop) and pass the count to the Python side.
    lim_out = fdr_out.replace("ed_fdr.txt", "ed_limited.txt")
    lim_out_r = lim_out.replace("\\", "/")
    rscript = f"""
suppressMessages({{library(Matrix); library(DropletUtils)}})
set.seed(0)
m <- readMM("{mtx_in_r}")                     # genes x barcodes
res <- emptyDrops(as(m, "CsparseMatrix"), lower={EMPTYDROPS_LOWER})
fdr <- res$FDR; fdr[is.na(fdr)] <- 1          # NA (ambient/below-lower) -> not a cell
lim <- res$Limited; lim[is.na(lim)] <- FALSE  # TRUE = p-value capped by niters
# barcodes that are 'Limited' AND still non-significant -> permutations too few
n_lim_ns <- sum(lim & (fdr > {EMPTYDROPS_FDR}), na.rm=TRUE)
if (n_lim_ns > 0) cat(sprintf("EMPTYDROPS_LIMITED_WARN %d\\n", n_lim_ns))
write.table(fdr, "{fdr_out_r}", row.names=FALSE, col.names=FALSE)
write.table(as.integer(n_lim_ns), "{lim_out_r}", row.names=FALSE, col.names=FALSE)
cat("EMPTYDROPS_OK\\n")
"""
    rf = os.path.join(tmp, "emptydrops.R")
    with open(rf, "w") as fh:
        fh.write(rscript)
    r = subprocess.run([RSCRIPT_BIN, rf], capture_output=True, text=True, timeout=10800)
    if "EMPTYDROPS_OK" not in r.stdout or not os.path.exists(fdr_out):
        raise RuntimeError(
            f"emptyDrops did NOT run (DropletUtils/Rscript env broken -- fix it, do not "
            f"fall back to min_genes on a RAW matrix).\n  stdout: {r.stdout[-300:]}\n"
            f"  stderr: {r.stderr[-800:]}")
    fdr = np.loadtxt(fdr_out)
    keep = np.asarray(fdr) <= EMPTYDROPS_FDR
    n_cells = int(keep.sum())
    # A4: surface the emptyDrops 'Limited' diagnostic -- if any Limited barcode is still
    # non-significant, the permutation count (niters) was too low (WARN only, no auto-loop).
    lim_file = fdr_out.replace("ed_fdr.txt", "ed_limited.txt")
    if os.path.exists(lim_file):
        try:
            n_lim_ns = int(np.loadtxt(lim_file))
        except Exception:
            n_lim_ns = 0
        if n_lim_ns > 0:
            log(f"    [WARN] emptyDrops: {n_lim_ns} barcodes are 'Limited' (Monte-Carlo "
                f"p capped by niters) yet non-significant at FDR<={EMPTYDROPS_FDR} -- "
                f"permutation count may be insufficient; consider raising emptyDrops niters.")
    log(f"    emptyDrops: {n_in} raw barcodes -> {n_cells} cells (FDR<={EMPTYDROPS_FDR})")
    return adata_sample[keep].copy(), n_in, n_cells


def scar_ambient(adata):
    """Pure-python ambient removal via scar (scvi-tools ecosystem). GPU-friendly.
    Returns denoised counts; falls back to raw if scar not importable."""
    try:
        from scar import model as scar_model  # noqa
        import torch  # noqa
        m = scar_model(raw_count=adata.X, ambient_profile=None, feature_type="mRNA",
                       device="cuda" if _has_cuda() else "cpu")
        m.train(epochs=150, verbose=False)
        m.inference()
        return sp.csr_matrix(np.asarray(m.native_counts)).astype("float32"), "scar"
    except Exception as e:
        log(f"    scar unavailable ({e}); raw fallback")
        return adata.X, "raw_fallback"


def _has_cuda():
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False


def apply_ambient(adata_sample, tmp):
    """Dispatch ambient correction per the AMBIENT_METHOD switch. Operates on ONE
    sample's raw counts. Returns (new_counts_csr, method_used)."""
    counts = adata_sample.X if sp.issparse(adata_sample.X) else sp.csr_matrix(adata_sample.X)
    if AMBIENT_METHOD == "decontx":
        z = _quick_z(adata_sample)   # pass clusters so decontX skips its fragile UMAP
        return decontx_R(counts.tocsr(), z, tmp)
    if AMBIENT_METHOD == "scar":
        return scar_ambient(adata_sample)
    return counts.tocsr(), "none"


# ----------------------------------------------------------- DOUBLET REMOVAL --
def flag_doublets(adata_sample):
    """Per-sample doublet call. Scrublet (python) default; scDblFinder (R) optional.
    Returns a boolean 'is_doublet' array aligned to adata_sample.obs_names."""
    n = adata_sample.n_obs
    if DOUBLET_METHOD == "none" or n < 30:
        return np.zeros(n, dtype=bool)
    if DOUBLET_METHOD == "scrublet":
        # scanpy's BUILT-IN Scrublet (sc.pp.scrublet) -- identical Scrublet algorithm,
        # actively maintained, and needs NO standalone 'scrublet' package (which is
        # unmaintained and absent in modern envs). NO except-swallow: a crash is an
        # env/input BUG to fix, not a reason to silently skip doublet removal
        # (skipping = doublets contaminate the atlas while we report "removed").
        As = adata_sample.copy()   # keep caller's counts untouched; scrublet writes obs
        try:
            sc.pp.scrublet(As, random_state=SEED)
        except (ValueError, ZeroDivisionError) as e:
            # ★2026-07-08: scrublet's internal gene/cell filtering (min_cells=3) can empty a very
            # sparse plate -> log1p on a (0,0) array (ValueError: 0 samples). This is a DATA condition
            # (a near-empty MARS-seq plate), NOT an env bug: LOUDLY skip doublet flagging for this
            # sample (same spirit as the n<30 guard); QC downstream drops it if nothing survives.
            log(f"    scrublet could not run (too-sparse sample, n={n}): {str(e)[:70]} "
                f"-> 0 doublets flagged (LOGGED, not silently skipped)")
            return np.zeros(n, dtype=bool)
        if "predicted_doublet" not in As.obs:
            raise RuntimeError("sc.pp.scrublet produced no predicted_doublet; investigate")
        pred = As.obs["predicted_doublet"].values
        if pd.isna(pred).any():
            # scrublet genuinely could NOT auto-threshold (bimodality absent in a
            # small/uniform sample). Real scrublet outcome, not a hidden error: apply
            # its own score at the standard 0.25 cut, LOGGED -- never silently keep all.
            sco = As.obs.get("doublet_score")
            pred = (np.asarray(sco.values) >= 0.25) if sco is not None else np.zeros(n, bool)
            log(f"    scrublet auto-threshold failed; fixed 0.25 cut "
                f"-> {int(np.sum(pred))}/{n} doublets (LOGGED, not skipped)")
        return np.asarray(pred, dtype=bool)
    if DOUBLET_METHOD == "scdblfinder":
        # NOT wired: raise instead of a silent no-op that would report "doublets removed"
        # while removing none (the forbidden reported-but-not-done trap). Wire an Rscript
        # scDblFinder hook here before selecting this backend, or use DOUBLET_METHOD=scrublet.
        raise NotImplementedError(
            "DOUBLET_METHOD='scdblfinder' is not wired in-loop -- use 'scrublet' or add the "
            "Rscript scDblFinder hook. Refusing to silently skip doublet removal.")
    raise ValueError(f"unknown DOUBLET_METHOD={DOUBLET_METHOD}")


# ---------------------------------------------- LITERATURE-GROUNDED CLEANING --
def rbc_filter(adata, log=log, hb_max_pct=HB_MAX_PCT, counts_layer=None):
    """Remove RBC/erythrocyte-contaminated barcodes by hemoglobin fraction. The Hb
    fraction is computed on COUNTS: pass counts_layer='counts' when adata.X is already
    log-normalised, else it uses .X (assumed counts). Erythrocytes express almost only
    Hb (HB_GENES); a cell with >hb_max_pct% of its counts from Hb is RBC contamination
    the MAD-mito filter misses. Returns (filtered_adata, n_removed). If Hb genes are
    absent it logs and returns unchanged (an honest, logged no-op, not a fallback)."""
    hb = [g for g in HB_GENES if g in adata.var_names]
    if not hb:
        log("    rbc_filter: no Hb genes present -> nothing to filter (logged)")
        return adata, 0
    X = adata.layers[counts_layer] if (counts_layer and counts_layer in adata.layers) else adata.X
    hb_idx = [adata.var_names.get_loc(g) for g in hb]
    hb_sum = np.asarray(X[:, hb_idx].sum(1)).ravel()
    tot = np.asarray(X.sum(1)).ravel()
    pct_hb = 100.0 * hb_sum / np.maximum(tot, 1.0)
    keep = pct_hb <= hb_max_pct
    n_removed = int((~keep).sum())
    log(f"    rbc_filter: removed {n_removed}/{adata.n_obs} cells with >{hb_max_pct:.0f}% Hb (RBC)")
    return adata[keep].copy(), n_removed


def score_dissoc_cellcycle(adata, log=log, use_raw=True):
    """Add literature-grounded per-cell scores to adata.obs (needs LOG-NORM data; uses
    .raw when use_raw). NO cells are removed and NOTHING is regressed -- these are
    reported/used as COVARIATES (esp. the dissociation score for the FOSB control in
    layer 07). Adds: dissoc_score (van den Brink 2017), S_score/G2M_score/phase
    (Tirosh 2016 cell cycle)."""
    src = adata.raw if (use_raw and adata.raw is not None) else adata
    dg = [g for g in DISSOC_SIG if g in src.var_names]
    sc.tl.score_genes(adata, dg, score_name="dissoc_score", use_raw=use_raw,
                      random_state=SEED)
    log(f"    dissociation score: {len(dg)}/{len(DISSOC_SIG)} genes "
        f"(median={float(np.median(adata.obs['dissoc_score'])):.3f}); FOSB present={'FOSB' in dg}")
    sg = [g for g in CC_S_GENES if g in src.var_names]
    g2 = [g for g in CC_G2M_GENES if g in src.var_names]
    # score_genes_cell_cycle scores on .X (log-norm); it needs the genes IN .X. If .X is
    # the HVG-subset we score on a full-gene copy so no cycle gene is silently dropped.
    if use_raw and adata.raw is not None and adata.raw.shape[1] != adata.n_vars:
        tmp = ad.AnnData(X=adata.raw.X, obs=adata.obs[[]].copy(),
                         var=pd.DataFrame(index=adata.raw.var_names))
        sc.tl.score_genes_cell_cycle(tmp, s_genes=sg, g2m_genes=g2)
        for cc in ("S_score", "G2M_score", "phase"):
            adata.obs[cc] = tmp.obs[cc].values
    else:
        sc.tl.score_genes_cell_cycle(adata, s_genes=sg, g2m_genes=g2)
    log(f"    cell-cycle: S={len(sg)} G2M={len(g2)} genes; "
        f"phase mix={adata.obs['phase'].value_counts().to_dict()}")
    return adata


# ==============================================================================
# STEP 1 -- per-cohort / per-sample: ambient + doublet + MAD QC -> clean h5ad
# ==============================================================================
def step1_qc_per_cohort():
    qc_rows = []
    for c in COHORTS:
        out_h5 = _ck(f"clean__{c['label']}.h5ad")
        if os.path.exists(out_h5):
            log(f"STEP1 skip {c['label']} (checkpoint exists)")
            # recover the qc summary from the sidecar if present
            side = _ck(f"qc__{c['label']}.json")
            if os.path.exists(side):
                qc_rows += json.load(open(side))
            continue
        log(f"STEP1 load {c['label']} ({c['kind']})")
        # NO try/except-skip: a load failure is a loader/format BUG to FIX, not to
        # silently skip (skipping a cohort = running on less data while reporting
        # success = the forbidden silent degradation). Let it raise with traceback.
        A = load_cohort(c)
        log(f"  raw: {A.shape}, samples={A.obs['sample'].nunique()}")
        A.var["mt"] = A.var_names.str.upper().str.startswith("MT-")
        A.var["ribo"] = A.var_names.str.upper().str.startswith(("RPS", "RPL"))

        cleaned_samples = []
        rows_this = []
        for s in A.obs["sample"].astype(str).unique():
            As = A[A.obs["sample"].astype(str) == s].copy()
            n0 = As.n_obs
            # skip degenerate tiny samples: DecontX's internal UMAP needs > n_neighbors
            # cells or it crashes ('n_neighbors must be smaller than the dataset size');
            # a <100-cell library cannot be reliably decontaminated/clustered/annotated.
            # LOGGED QC exclusion (reportable), not a silent fallback. (MARS plates exempt.)
            if c["kind"] != "mars" and n0 < 100:
                log(f"    {s}: n={n0} < 100 -> excluded (too small for DecontX/annotation, logged)")
                rows_this.append(dict(cohort=c["label"], sample=s, compartment=c["compartment"],
                                      n_raw=int(n0), n_after_doublet=0, n_doublet=0, n_final=0,
                                      ambient_method="skip_too_small", doublet_method=DOUBLET_METHOD,
                                      gene_lo=0, gene_hi=0, umi_lo=0, umi_hi=0, mito_hi=0))
                del As; gc.collect(); continue
            # A13: the ambient(a) -> doublet(b) -> QC(c) ORDER is BY DESIGN, not incidental.
            # DecontX/scAR estimate the ambient profile from the FULL set of (already
            # cell-called) barcodes, so ambient decontamination must run FIRST, on the
            # complete ambient-containing barcode set -- before doublets and QC prune cells
            # and distort that profile. Doublet detection then runs on decontaminated counts,
            # and data-driven MAD QC runs last on the cleaned matrix.
            # (a) ambient -- PLATFORM-AWARE. DecontX/scAR model DROPLET ambient RNA
            # (soup from cell-free transcripts co-encapsulated in GEMs). MARS-seq is
            # plate-based FACS-sorted (one cell sorted per well): there is no shared
            # droplet soup, so droplet-ambient decontamination is inappropriate and
            # its assumptions don't hold. For MARS cohorts we skip it (multiplets are
            # instead removed upstream via the plate metadata Number_of_cells==1), an
            # honest LOGGED choice -- not a silent fallback.
            if c["kind"] == "mars":
                dec_counts, amb_used = (As.X if sp.issparse(As.X)
                                        else sp.csr_matrix(As.X)), "skip_mars_plate"
            else:
                with tempfile.TemporaryDirectory() as tmp:
                    dec_counts, amb_used = apply_ambient(As, tmp)
            As.layers["counts_raw"] = (As.X.copy() if sp.issparse(As.X) else sp.csr_matrix(As.X))
            As.X = dec_counts if sp.issparse(dec_counts) else sp.csr_matrix(dec_counts)
            As.layers["counts"] = As.X.copy()   # decontaminated counts = canonical

            # (b) doublets (on decontaminated counts)
            is_dbl = flag_doublets(As)
            As = As[~is_dbl].copy()
            n_after_dbl = As.n_obs

            # (c) data-driven QC (MAD on nGene + nUMI; MAD-upper on mito, hard ceiling)
            sc.pp.calculate_qc_metrics(As, qc_vars=["mt", "ribo"], inplace=True, percent_top=None)
            g_lo, g_hi = mad_bounds(As.obs["n_genes_by_counts"].values)
            u_lo, u_hi = mad_bounds(As.obs["total_counts"].values)
            _, mt_hi   = mad_bounds(As.obs["pct_counts_mt"].values, n_mads=MITO_N_MADS, log1p=False)  # A3: tighter 3-MAD for mito
            mt_hi = MITO_HARD_MAX if not np.isfinite(mt_hi) else min(max(mt_hi, 5.0), MITO_HARD_MAX)  # ★2026-07-08 nan-safe; never below 5%, never above ceiling
            g_lo = max(g_lo, MIN_GENES)
            keep = ((As.obs["n_genes_by_counts"] >= g_lo) & (As.obs["n_genes_by_counts"] <= g_hi) &
                    (As.obs["total_counts"] >= u_lo) & (As.obs["total_counts"] <= u_hi) &
                    (As.obs["pct_counts_mt"] <= mt_hi)).values
            As = As[keep].copy()

            rows_this.append(dict(cohort=c["label"], sample=s, compartment=c["compartment"],
                                  n_raw=int(n0), n_after_doublet=int(n_after_dbl),
                                  n_doublet=int(n0 - n_after_dbl), n_final=int(As.n_obs),
                                  ambient_method=amb_used, doublet_method=DOUBLET_METHOD,
                                  gene_lo=round(g_lo, 1), gene_hi=round(g_hi, 1),
                                  umi_lo=round(u_lo, 1), umi_hi=round(u_hi, 1),
                                  mito_hi=round(mt_hi, 2)))
            if As.n_obs >= 20:
                cleaned_samples.append(As)
            log(f"    {s}: {n0} -> dbl {n0-n_after_dbl} -> QC keep {As.n_obs} "
                f"[amb={amb_used}, gene {g_lo:.0f}-{g_hi:.0f}, mito<={mt_hi:.1f}]")
            del As; gc.collect()

        if not cleaned_samples:
            log(f"  !! {c['label']}: nothing survived QC -- skip"); continue
        Ac = ad.concat(cleaned_samples, join="outer", index_unique=None)
        Ac.var_names_make_unique()
        sc.pp.filter_genes(Ac, min_cells=MIN_CELLS)
        Ac.obs["cohort"] = c["label"]; Ac.obs["compartment_src"] = c["compartment"]
        Ac.write(out_h5, compression="gzip")
        json.dump(rows_this, open(_ck(f"qc__{c['label']}.json"), "w"))
        qc_rows += rows_this
        log(f"  STEP1 wrote {out_h5}: {Ac.shape}")
        del Ac, cleaned_samples; gc.collect()

    pd.DataFrame(qc_rows).to_csv(os.path.join(OUT, "per_sample_qc_report.csv"), index=False)
    log("STEP1 done -> per_sample_qc_report.csv")


# ==============================================================================
# STEP 2 -- merge all cleaned cohorts, normalize, HVG, PCA, Harmony (scVI optional)
# ==============================================================================
def step2_integrate():
    out_h5 = _ck("integrated_embedding.h5ad")
    if os.path.exists(out_h5):
        log("STEP2 skip (integrated_embedding.h5ad exists)")
        return out_h5
    parts = sorted(glob.glob(_ck("clean__*.h5ad")))
    if not parts:
        raise FileNotFoundError("STEP2: no clean__*.h5ad checkpoints from STEP1")
    log(f"STEP2 merge {len(parts)} cleaned cohorts")
    ads = [sc.read_h5ad(p) for p in parts]
    adata = ad.concat(ads, join="outer", index_unique=None)
    adata.var_names_make_unique()
    del ads; gc.collect()
    # outer-join can introduce NaN for genes absent in a cohort -> 0 counts
    if sp.issparse(adata.X):
        adata.X.data = np.nan_to_num(adata.X.data)
    log(f"  merged: {adata.shape}, cohorts={adata.obs['cohort'].nunique()}, "
        f"samples={adata.obs['sample'].nunique()}")

    # keep raw decontaminated counts before normalization
    adata.layers["counts"] = (adata.X.copy() if sp.issparse(adata.X) else sp.csr_matrix(adata.X))
    # A13: target_sum=1e4 is NOT arbitrary -- it is a HARD CellTypist INPUT CONSTRAINT
    # (models are trained on log1p counts-per-10k; a different target_sum silently mis-
    # scales the reference-mapping call). adata.raw (set below) becomes the CellTypist input.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata   # log-norm full gene space, needed by score_genes(use_raw=True)

    # HVG computed batch-aware (per cohort) so no single cohort dominates
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat",
                                batch_key="cohort")
    h = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(h, max_value=10)
    sc.tl.pca(h, n_comps=N_PCS, random_state=SEED)
    adata.obsm["X_pca"] = h.obsm["X_pca"]

    # ---- Harmony (harmonypy, CPU) -- default embedding -----------------------
    # BATCH KEY = cohort (study) only, NOT cohort x sample. Rationale:
    #   * cohort (7 studies) is the dominant technical/platform batch axis in a
    #     multi-study atlas -- Harmony on study is the standard integration level.
    #   * "sample" gives 776 groups, but ONE cohort (Gur_GSE195452) contributes 665 of
    #     them, each tiny (median ~77 cells) -- these are almost certainly genetically
    #     de-multiplexed INDIVIDUALS from pooled runs, i.e. biological units, NOT separate
    #     technical batches. Correcting them as 776 batches would (a) over-integrate and
    #     erase inter-patient biology we NEED downstream (05 regulon-severity), (b) rely on
    #     unreliable per-batch statistics from <100-cell groups, and (c) cost ~10 min/Harmony
    #     iteration (~1.7 h). cohort-level runs in seconds and preserves patient-level signal.
    #   (If sample-level residual batch is later a concern, revisit with scVI/scanorama on
    #    size-filtered samples -- do not naively split multiplexed cohorts into 665 batches.)
    adata.obs["batch"] = adata.obs["cohort"].astype("category")   # integration batch = study
    log("STEP2 Harmony integration (harmonypy, CPU) on cohort (study) batch")
    # Direct harmonypy call (NOT sc.external.pp.harmony_integrate): the scanpy wrapper
    # crashed opaquely at 424K with a degenerate Z_corr shape (it fed a float32/non-contiguous
    # obsm view straight in). Here we hand harmonypy a contiguous float64 copy, log the exact
    # input/output, and handle both Z_corr orientations so any failure is diagnosable.
    import harmonypy
    _xpca = np.ascontiguousarray(np.asarray(adata.obsm["X_pca"], dtype=np.float64))
    log(f"  [harmony-in] X_pca shape={_xpca.shape} dtype={_xpca.dtype} "
        f"nan={int(np.isnan(_xpca).sum())} inf={int(np.isinf(_xpca).sum())} "
        f"cohort={adata.obs['cohort'].nunique()} sample={adata.obs['sample'].nunique()} "
        f"(integrating on cohort)")
    _ho = harmonypy.run_harmony(_xpca, adata.obs, ["cohort"], random_state=SEED)
    _Z = np.asarray(_ho.Z_corr)
    log(f"  [harmony-out] Z_corr shape={_Z.shape} (n_obs={adata.n_obs}, n_pcs={N_PCS})")
    if _Z.shape == (adata.n_obs, N_PCS):
        adata.obsm["X_pca_harmony"] = _Z
    elif _Z.shape == (N_PCS, adata.n_obs):
        adata.obsm["X_pca_harmony"] = _Z.T
    else:
        raise RuntimeError(
            f"Harmony returned degenerate Z_corr {_Z.shape}; X_pca in was {_xpca.shape}. "
            "Neither orientation matches (n_obs,n_pcs) -- harmonypy/data problem, not the wrapper.")
    INTEG_REP = "X_pca_harmony"

    # ---- scVI (OPTIONAL, GPU) -- library module 506; compared, not assumed ----
    # Honest-baseline note (Genome Biol 2025): learned embeddings do NOT auto-beat
    # Harmony. If USE_SCVI, we ALSO compute X_scVI and store both; downstream can
    # pick by the batch-mixing / biology-conservation metric pair, not the prettier UMAP.
    if USE_SCVI:
        try:
            import scvi
            scvi.settings.seed = SEED
            scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
            m = scvi.model.SCVI(adata, n_latent=30, n_layers=2, gene_likelihood="nb")
            m.train(max_epochs=SCVI_EPOCHS,
                    accelerator="gpu" if _has_cuda() else "cpu",
                    enable_progress_bar=False)
            adata.obsm["X_scVI"] = m.get_latent_representation()
            log("STEP2 scVI latent computed (X_scVI) -- optional GPU path")
            # let scVI drive clustering only if a GPU actually ran it
            INTEG_REP = "X_scVI" if _has_cuda() else INTEG_REP
        except Exception as e:
            log(f"STEP2 scVI optional path skipped ({e}); staying on Harmony")

    adata.uns["integration_rep"] = INTEG_REP
    sc.pp.neighbors(adata, use_rep=INTEG_REP, n_neighbors=15, random_state=SEED)
    # provisional clustering only (the from-scratch STEP3 re-does the resolution
    # sweep + stability pick and OVERWRITES obs['leiden']); igraph flavor = the
    # scanpy-recommended default, avoids the leidenalg FutureWarning.
    sc.tl.leiden(adata, resolution=LEIDEN_RES, random_state=SEED, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=SEED)
    log(f"  clusters: {adata.obs['leiden'].nunique()} (rep={INTEG_REP})")
    adata.write(out_h5, compression="gzip")
    log(f"STEP2 wrote {out_h5}")
    del h; gc.collect()
    return out_h5


# ==============================================================================
# ANNOTATION HELPERS -- shared by STEP3 (compartments) and STEP4 (fibro subtypes).
# All from-scratch: no pre-existing label is ever read. Every call below is a real,
# verified scanpy / celltypist API (signatures checked 2026-07-02; no invented args).
# ==============================================================================
def leiden_stability_pick(adata, res_sweep, key_prefix="leiden_r", nn_ready=True):
    """Human-like RESOLUTION SWEEP + STABILITY pick (not a single hard-coded res).

    Cluster at every resolution in `res_sweep` on the ALREADY-BUILT neighbour graph,
    then pick the resolution whose partition is most STABLE, measured as the mean
    Adjusted Rand Index (ARI) to its two neighbouring resolutions in the sweep. The
    intuition (standard practice, e.g. sc-best-practices clustering chapter): a
    resolution sitting on a stable plateau -- where nudging the parameter barely
    changes the partition -- is a more defensible clustering than a knife-edge one.
    Ties / un-scorable sweeps fall back to LEIDEN_RES. Returns (best_res, best_key).
    Stores every partition in obs[key_prefix+res] so the audit can show them all.
    """
    from sklearn.metrics import adjusted_rand_score
    keys = []
    for r in res_sweep:
        k = f"{key_prefix}{r}"
        sc.tl.leiden(adata, resolution=r, random_state=SEED, key_added=k,
                     flavor="igraph", n_iterations=2, directed=False)
        keys.append((r, k))
        log(f"    leiden res={r}: {adata.obs[k].nunique()} clusters")
    if len(keys) < 3:
        best = LEIDEN_RES if any(r == LEIDEN_RES for r, _ in keys) else keys[0][0]
        return best, f"{key_prefix}{best}"
    # mean ARI of each interior resolution to its two neighbours
    scores = {}
    for i in range(1, len(keys) - 1):
        r, k = keys[i]
        a_prev = adjusted_rand_score(adata.obs[keys[i - 1][1]], adata.obs[k])
        a_next = adjusted_rand_score(adata.obs[keys[i + 1][1]], adata.obs[k])
        scores[r] = float(np.mean([a_prev, a_next]))
    best_res = max(scores, key=scores.get)
    log("    stability (mean ARI to neighbours): " +
        json.dumps({str(k): round(v, 3) for k, v in scores.items()}))
    log(f"    -> stability-picked resolution = {best_res}")
    adata.uns["leiden_stability_scores"] = {str(k): v for k, v in scores.items()}
    return best_res, f"{key_prefix}{best_res}"


def marker_score_matrix(adata, panels, cluster_key, prefix):
    """Score each panel with sc.tl.score_genes (control-gene-normalised, use_raw),
    average per cluster, then Z-SCORE each panel ACROSS clusters so panels with
    different baseline magnitudes are comparable. Returns (zscored_df, argmax_call,
    present_panels). argmax over z-scored panels = the MARKER call per cluster."""
    present = []
    for name, gs in panels.items():
        gg = [g for g in gs if g in adata.raw.var_names]
        if len(gg) >= 2:                       # need >=2 present genes to trust a panel
            sc.tl.score_genes(adata, gg, score_name=f"{prefix}{name}", use_raw=True,
                              random_state=SEED)
            present.append(name)
        else:
            log(f"    [marker] panel '{name}' skipped: <2 of its genes present")
    raw_mean = (adata.obs[[cluster_key] + [f"{prefix}{n}" for n in present]]
                .groupby(cluster_key, observed=True).mean())
    raw_mean.columns = present
    # z-score each panel column across clusters (population std; guard zero-var)
    z = raw_mean.copy()
    for c in present:
        col = raw_mean[c].values.astype(float)
        sd = col.std()
        z[c] = (col - col.mean()) / sd if sd > 1e-9 else 0.0
    argmax_call = {cl: z.loc[cl].idxmax() for cl in z.index}
    return z, argmax_call, present


def celltypist_reference_labels(adata, cluster_key):
    """Independent REFERENCE-MAPPING call via CellTypist (verified real API).

    Runs each model in CELLTYPIST_MODELS with majority_voting refined WITHIN our own
    Leiden clusters (over_clustering=cluster_key). Returns a dict:
        per_cell   -> DataFrame (index=cells) with one majority_voting column/model
        per_cluster-> DataFrame (index=clusters) with, per model, the modal
                      majority-voting label and its within-cluster fraction (conf).
    Input is a FRESH log1p-1e4 AnnData built from adata.raw (the required format).
    Each model is isolated in try/except: a failure is LOGGED and that model column
    is dropped -- annotation never hard-fails on a missing/renamed model.
    """
    if not USE_CELLTYPIST:
        log("    CellTypist disabled (USE_CELLTYPIST=False) -> marker-only annotation")
        return None
    try:
        import celltypist
        from celltypist import models as ct_models
    except Exception as e:
        log(f"    CellTypist not importable ({e}); marker-only annotation")
        return None

    # Build the exact input CellTypist wants: log1p-normalised-to-1e4 in .X.
    # adata.raw was set in STEP2 to the full log1p(normalize_total 1e4) matrix.
    ref = ad.AnnData(X=adata.raw.X.copy(),
                     obs=adata.obs[[cluster_key]].copy(),
                     var=pd.DataFrame(index=adata.raw.var_names))
    per_cell = pd.DataFrame(index=adata.obs_names)
    per_cluster = pd.DataFrame(index=sorted(adata.obs[cluster_key].astype(str).unique()))
    for mdl in CELLTYPIST_MODELS:
        try:
            try:
                ct_models.download_models(model=[mdl], force_update=False)
            except Exception as de:
                log(f"    [celltypist] download check for {mdl}: {de} (using cached if present)")
            res = celltypist.annotate(ref, model=mdl, majority_voting=True,
                                      over_clustering=cluster_key)
            pl = res.predicted_labels
            col = "majority_voting" if "majority_voting" in pl.columns else "predicted_labels"
            tag = mdl.replace(".pkl", "")
            lab = pl[col].astype(str)
            per_cell[f"ct_{tag}"] = lab.reindex(adata.obs_names).values
            # N2: per-cell confidence = max class probability from the SAME annotate call
            # (no extra compute). The majority_voting `frac` below is ~1.0 by construction
            # and is NOT a real per-cell confidence; this max-probability is. A capture
            # failure only drops the diagnostic column, never the model's label.
            try:
                pm = res.probability_matrix
                conf = np.asarray(pm.to_numpy().max(axis=1)).ravel()
                per_cell[f"conf_{tag}"] = pd.Series(conf, index=pm.index).reindex(adata.obs_names).values
            except Exception as pe:
                log(f"    [celltypist] {mdl} no per-cell probability_matrix ({pe})")
            # modal label per Leiden cluster + its within-cluster fraction (confidence)
            grp = pd.DataFrame({"cl": adata.obs[cluster_key].astype(str).values,
                                "lab": lab.reindex(adata.obs_names).values})
            modal, frac = {}, {}
            for cl, sub in grp.groupby("cl"):
                vc = sub["lab"].value_counts()
                modal[cl] = vc.index[0]
                frac[cl] = round(float(vc.iloc[0]) / float(len(sub)), 3)
            per_cluster[f"ct_{tag}"] = pd.Series(modal)
            per_cluster[f"ct_{tag}_frac"] = pd.Series(frac)
            log(f"    [celltypist] {mdl}: mapped {lab.nunique()} distinct leaf labels")
        except Exception as e:
            log(f"    [celltypist] model {mdl} FAILED ({e}); skipped\n{traceback.format_exc()}")
    if per_cell.shape[1] == 0:
        log("    CellTypist produced no usable model -> marker-only annotation")
        return None
    return dict(per_cell=per_cell, per_cluster=per_cluster)


def map_to_vocab(label):
    """Collapse a free-text CellTypist leaf label onto the coarse SHARED_VOCAB.
    Case-insensitive substring/keyword match; unmatched -> 'Other' (never forced)."""
    if not isinstance(label, str) or label.strip() == "":
        return "NA"
    s = f" {label.lower()} "
    for vocab, keys in SHARED_VOCAB.items():
        for kw in keys:
            if kw in s:
                return vocab
    return "Other"


# ---------------------------------------------- M1: INTEGRATION BENCHMARK ------
def integration_benchmark(adata):
    """M1 (top-journal must-have): QUANTITATIVE integration benchmark via scib-metrics
    on the stored embeddings -- unintegrated (X_pca) vs Harmony (X_pca_harmony) vs scVI
    (X_scVI, if computed) -- keyed on the from-scratch celltype labels. Reports batch-
    mixing (iLISI/kBET/graph-connectivity) + bio-conservation (cLISI/ARI/NMI/ASW) so the
    integration choice is MEASURED, not asserted (Luecken 2019 / Heumos 2023). Writes
    integration_benchmark.csv. If scib-metrics is not installed it logs LOUDLY and skips
    (this is a reporting deliverable, not a correctness gate -- install scib-metrics on
    the server to produce it); NOT a silent fallback."""
    try:
        from scib_metrics.benchmark import Benchmarker
    except Exception as e:
        log(f"  [M1] scib-metrics NOT installed ({e}) -> integration benchmark SKIPPED. "
            f"Run `pip install scib-metrics` on the server to produce integration_benchmark.csv")
        return
    embs = [k for k in ("X_pca", "X_pca_harmony", "X_scVI") if k in adata.obsm]
    if len(embs) < 2 or "batch" not in adata.obs or "celltype" not in adata.obs:
        log(f"  [M1] need >=2 embeddings + batch + celltype; have {embs} -> skip")
        return
    # scib-metrics silhouette/LISI are O(n^2); on the full 424K atlas the benchmark runs for
    # HOURS. STANDARD PRACTICE (Luecken 2022 scib) is to evaluate integration on a representative
    # SUBSAMPLE. Stratify by celltype so rare labels stay represented for the isolated-labels /
    # label-ASW metrics; the resulting scores are near-identical to full-data but computed in minutes.
    N_BENCH = 40000
    if adata.n_obs > N_BENCH:
        rng = np.random.default_rng(SEED)
        cts = adata.obs["celltype"].astype(str).values
        keep = []
        for ct in np.unique(cts):
            pos = np.where(cts == ct)[0]
            take = min(len(pos), max(20, int(round(len(pos) * N_BENCH / adata.n_obs))))
            keep.append(rng.choice(pos, size=take, replace=False))
        sub = np.sort(np.concatenate(keep))
        adata_bench = adata[sub].copy()
        log(f"  [M1] scib benchmark on STRATIFIED SUBSAMPLE {adata_bench.n_obs}/{adata.n_obs} cells "
            f"(O(n^2) metrics -> subsample is scib best-practice; seed={SEED})")
    else:
        adata_bench = adata
    log(f"  [M1] scib benchmark on {embs} (batch_key=batch, label_key=celltype)...")
    bm = Benchmarker(adata_bench, batch_key="batch", label_key="celltype",
                     embedding_obsm_keys=embs, n_jobs=-1)
    bm.benchmark()
    res = bm.get_results(min_max_scale=False)
    res.to_csv(os.path.join(OUT, "integration_benchmark.csv"))
    log(f"  [M1] wrote integration_benchmark.csv (n={adata_bench.n_obs}):\n" + res.to_string())


# ==============================================================================
# STEP 3 -- FROM-SCRATCH compartment annotation by MULTI-PANEL + REF-MAP CONSENSUS.
# NO prior label is read. Marker argmax (multi-panel, z-scored) AND CellTypist
# reference-mapping are computed independently, mapped to a shared vocabulary, and
# reconciled into a consensus with a needs_review flag. Full audit CSVs emitted.
# ==============================================================================
def step3_annotate(integ_h5):
    out_full = os.path.join(OUT, "full_integrated_annotated.h5ad")
    if os.path.exists(out_full):
        log("STEP3 skip (full_integrated_annotated.h5ad exists)")
        return out_full
    adata = sc.read_h5ad(integ_h5)
    log(f"STEP3 annotate FROM SCRATCH {adata.shape} (no prior labels trusted)")

    # ---- (0) LITERATURE-GROUNDED CLEANING (RBC removal + covariate scores) ----
    # RBC/erythrocyte removal (Hb fraction on the counts layer) BEFORE clustering so
    # contamination does not seed a spurious cluster; then per-cell dissociation
    # (van den Brink 2017) + cell-cycle (Tirosh 2016) scores kept as covariates
    # (the dissociation score is the FOSB-artifact control consumed by layer 07).
    adata, n_rbc = rbc_filter(adata, counts_layer="counts")
    adata = score_dissoc_cellcycle(adata, use_raw=False)   # .X here is already log-norm
    adata.uns["n_rbc_removed"] = int(n_rbc)

    # ---- (1) RESOLUTION SWEEP + STABILITY PICK -------------------------------
    # STEP2 already built neighbours on the integration rep. Sweep resolutions and
    # pick the most stable partition; use THAT as the annotation clustering.
    best_res, best_key = leiden_stability_pick(adata, LEIDEN_RES_SWEEP,
                                               key_prefix="leiden_r")
    adata.obs["leiden"] = adata.obs[best_key].astype("category")
    adata.uns["annotation_resolution"] = float(best_res)
    log(f"  annotation clustering: res={best_res} -> {adata.obs['leiden'].nunique()} clusters")

    # ---- (2a) MARKER CALL: multi-panel, z-scored across clusters, argmax --------
    zmat, marker_call, present = marker_score_matrix(adata, MARKERS, "leiden", "sc_")
    zmat.round(3).to_csv(os.path.join(OUT, "annotation_marker_scores.csv"))
    log(f"  marker panels used ({len(present)}): {present}")

    # ---- (2b) REFERENCE-MAPPING CALL: CellTypist (independent) ------------------
    ctres = celltypist_reference_labels(adata, "leiden")

    # ---- (2c) per-cluster DE audit (rank_genes_groups) --------------------------
    # Wilcoxon on the log-norm data so a human can read what each cluster actually is.
    # DEFERRABLE: this full-atlas Wilcoxon (~47min on 424K) ONLY feeds the descriptive
    # annotation_cluster_DE.csv -- nothing downstream reads uns['rgg'] (consensus annotation
    # below uses marker-argmax + celltypist). SSC_SKIP_RANKGENES=1 skips it so the atlas
    # write / STEP4 / 03 unblock fast; rank_genes_standalone.py backfills it in parallel.
    if os.environ.get("SSC_SKIP_RANKGENES") == "1":
        log("  [2c] rank_genes_groups(leiden) SKIPPED (SSC_SKIP_RANKGENES=1) -> backfill via rank_genes_standalone.py")
        de_names = pd.DataFrame()   # keep DEFINED (used later for consensus top_DE); blank until standalone backfills
    else:
        sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True,
                                n_genes=25, key_added="rgg")
        de_names = pd.DataFrame(adata.uns["rgg"]["names"])
        de_rows = []
        for cl in de_names.columns:
            de_rows.append(dict(cluster=cl, top_DE_genes=", ".join(de_names[cl].tolist()[:25])))
        pd.DataFrame(de_rows).to_csv(os.path.join(OUT, "annotation_cluster_DE.csv"), index=False)

    # ---- (3) CONSENSUS: reconcile marker argmax vs reference-mapping ------------
    clusters = list(zmat.index)
    rows = []
    # Collect EACH reference model's vocab-mapped call per cluster. The models are
    # COMPLEMENTARY, not competing: Adult_Human_Skin resolves structural/stromal
    # lineages (fibro/pericyte/SMC/endo/melanocyte/KC), Immune_All_Low resolves
    # immune fine-types. A marker call is CONFIRMED if ANY reference model agrees;
    # a wrong call from the model that does NOT cover this lineage (e.g. the
    # immune-only model labelling a smooth-muscle cluster "T cell") therefore no
    # longer spuriously flags a correct marker call. (Previously the code kept only
    # the single highest-fraction model; because majority_voting makes every model
    # frac=1.0, the FIRST-listed model always won, so the immune model over-wrote
    # every non-immune structural call -> 10/45 false needs_review flags.)
    ct_vocabs_by_cluster, ct_detail_by_cluster = {}, {}
    if ctres is not None:
        pc = ctres["per_cluster"]
        pc.to_csv(os.path.join(OUT, "annotation_celltypist.csv"))
        model_cols = [c for c in pc.columns if not c.endswith("_frac")]
        for cl in clusters:
            vocs, detail = [], []
            for mc in model_cols:
                if cl in pc.index and isinstance(pc.loc[cl, mc], str):
                    raw = pc.loc[cl, mc]
                    fr = pc.loc[cl, mc + "_frac"] if (mc + "_frac") in pc.columns else np.nan
                    v = map_to_vocab(raw)
                    detail.append(f"{mc}={raw}({v},{fr})")
                    if v not in ("NA", "Other"):
                        vocs.append(v)
            ct_vocabs_by_cluster[cl] = vocs
            ct_detail_by_cluster[cl] = "; ".join(detail)
    else:
        for cl in clusters:
            ct_vocabs_by_cluster[cl] = []
            ct_detail_by_cluster[cl] = "celltypist_unavailable"

    for cl in clusters:
        mk = marker_call[cl]
        vocs = ct_vocabs_by_cluster[cl]
        # consensus rule (human-like): if ANY reference model confirms the marker
        # call -> confident. If no model produced a comparable call -> take marker
        # call, flag for review. If reference model(s) gave real calls but NONE
        # match -> take marker call, flag REVIEW (genuine ambiguity to inspect).
        if not vocs:
            rf, final, needs_review, why = "NA", mk, True, "no comparable reference call"
        elif mk in vocs:
            rf, final, needs_review, why = mk, mk, False, "marker confirmed by reference"
        else:
            rf = "/".join(sorted(set(vocs)))
            final, needs_review, why = mk, True, f"marker({mk}) != reference({rf})"
        rows.append(dict(cluster=cl, marker_call=mk, reference_call=rf,
                         final_celltype=final, needs_review=needs_review, reason=why,
                         marker_top_z=round(float(zmat.loc[cl].max()), 3),
                         celltypist_detail=ct_detail_by_cluster[cl],
                         top_DE=", ".join(de_names[cl].tolist()[:10]) if cl in de_names.columns else ""))
    consensus = pd.DataFrame(rows).set_index("cluster")
    consensus.to_csv(os.path.join(OUT, "annotation_consensus.csv"))

    assign = consensus["final_celltype"].to_dict()
    adata.obs["celltype"] = adata.obs["leiden"].map(assign).astype(str)
    adata.obs["celltype_needs_review"] = adata.obs["leiden"].map(
        consensus["needs_review"].to_dict()).astype(bool)

    # immune FINE subtype (literature-grounded refinement): for immune-compartment
    # clusters, attach the Immune_All_High per-cluster label (Dominguez Conde 2022
    # Science) -- the published fine reference that resolves the marker-unsupported
    # immune fine-types (NK-vs-cytotoxic-T, B, pDC) instead of a coarse marker call.
    IMMUNE_COMPARTMENTS = {"Tcell", "NK", "Bcell", "Plasma", "Myeloid", "Mast", "pDC"}
    adata.obs["immune_subtype"] = ""
    if ctres is not None and "ct_Immune_All_High" in ctres["per_cluster"].columns:
        hi = ctres["per_cluster"]["ct_Immune_All_High"].astype(str).to_dict()
        im = adata.obs["celltype"].isin(IMMUNE_COMPARTMENTS).values
        adata.obs.loc[im, "immune_subtype"] = (
            adata.obs.loc[im, "leiden"].astype(str).map(hi).values)

    # N2: per-cell CellTypist confidence (max class probability across models) + a
    # low-confidence 'Unknown' gate. A cell can be low-confidence even inside a confident
    # cluster; the majority_voting frac (~1.0) is NOT this. celltype is KEPT; low_confidence
    # is a flag so downstream/figures can mark or exclude these cells.
    adata.obs["celltypist_conf"] = np.nan
    adata.obs["low_confidence"] = False
    if ctres is not None:
        conf_cols = [c for c in ctres["per_cell"].columns if c.startswith("conf_")]
        if conf_cols:
            conf = ctres["per_cell"][conf_cols].reindex(adata.obs_names).max(axis=1)
            adata.obs["celltypist_conf"] = conf.values
            adata.obs["low_confidence"] = (conf < CELLTYPIST_CONF_MIN).fillna(True).values
            log(f"  low-confidence cells (CellTypist max prob < {CELLTYPIST_CONF_MIN}): "
                f"{int(adata.obs['low_confidence'].sum())}/{adata.n_obs} (flagged, celltype kept)")

    # ---- per-cluster VERDICT (printed, human-readable) -------------------------
    log("  ===== PER-CLUSTER ANNOTATION VERDICT (confirm the flagged ones) =====")
    for cl in clusters:
        r = consensus.loc[cl]
        flag = "  <-- NEEDS_REVIEW" if r["needs_review"] else ""
        log(f"    cluster {cl:>3}: {r['final_celltype']:<13} "
            f"[marker={r['marker_call']} | ref={r['reference_call']}] "
            f"top_DE: {r['top_DE']}{flag}")
    n_flag = int(consensus["needs_review"].sum())
    log(f"  consensus: {len(clusters)} clusters, {n_flag} flagged needs_review")
    log("  celltype counts: " + json.dumps(adata.obs.celltype.value_counts().to_dict()))

    # ---- myofibroblast SIGNATURE + MYO program (unchanged contract downstream) --
    mg = [g for g in MYOFIB_SIG if g in adata.raw.var_names]
    sc.tl.score_genes(adata, mg, score_name="myofib_signature", use_raw=True,
                      random_state=SEED)
    pg = [g for g in MYO if g in adata.raw.var_names]
    adata.obs["myo"] = np.asarray(adata.raw[:, pg].X.mean(1)).ravel()

    # attach per-cell CellTypist raw labels for the audit trail (if available)
    if ctres is not None:
        for c in ctres["per_cell"].columns:
            adata.obs[c] = ctres["per_cell"][c].reindex(adata.obs_names).astype(str).values

    # M1: quantitative integration benchmark (scib-metrics) on the annotated atlas
    # (embeddings + batch + celltype now all present). Server-only; skips w/ loud log locally.
    # NON-FATAL: benchmark is a reporting deliverable, NOT a correctness gate -- a failure here
    # must NOT discard the (expensive) annotation + rank_genes or block the atlas write/STEP4.
    try:
        integration_benchmark(adata)
    except Exception as _bench_e:
        log(f"  [M1] integration benchmark FAILED ({type(_bench_e).__name__}: {_bench_e}) -> SKIPPED "
            f"(non-gate deliverable; atlas write proceeds so STEP4/downstream is NOT blocked). "
            f"Traceback:\n{traceback.format_exc()}")

    adata.write(out_full, compression="gzip")
    log(f"STEP3 wrote {out_full}  (+ 4 audit CSVs)")
    return out_full


# ==============================================================================
# STEP 4 -- fibroblast subset + subtype annotation (Tabib/Gur states) + purity
# ==============================================================================
def _purity_filter(adata, mito_max=0.15):
    """Per-cell fibroblast-purity filter: keep cells whose fibroblast-module score
    >= epidermal-module score. (★2026-07-07: the extra hard-coded mito<15% cut was
    REMOVED -- cells already passed the main QC's per-sample 3-MAD + 20% ceiling; a
    third inconsistent mito threshold silently changed the myofibroblast count fed to 03.)"""
    EPI = ["KRT1", "KRT10", "KRT14", "KRT5", "KRT15", "DMKN", "KRTDAP", "SFN",
           "LGALS7B", "KRT6A", "CALML5"]
    FIB = ["DCN", "LUM", "COL1A1", "COL1A2", "PDGFRA", "COL3A1"]
    src = adata.raw if adata.raw is not None else adata
    epi = [g for g in EPI if g in src.var_names]
    fib = [g for g in FIB if g in src.var_names]
    sc.tl.score_genes(adata, epi, score_name="_epi", use_raw=adata.raw is not None)
    sc.tl.score_genes(adata, fib, score_name="_fib", use_raw=adata.raw is not None)
    keep = adata.obs["_fib"].values >= adata.obs["_epi"].values
    # ★2026-07-07 修: 去掉此处第二道硬编码 mito<15% 切——细胞已在主 QC 过逐样本 3-MAD + 20%
    # 天花板; 再套无依据的 15% 是第三条互不一致的线粒体阈值, 会静默改变喂给 03 的肌成纤维数。
    # 纯度过滤只保留 fib>=epi 判据(mito_max 参数保留仅为向后兼容, 不再施用)。
    rep = dict(n_in=int(adata.n_obs), n_keep=int(keep.sum()),
               pct_drop=round(100 * (~keep).mean(), 2))
    out = adata[keep].copy()
    # drop the scratch score columns; DataFrame.pop() has NO default arg (dict-style
    # .pop(cc, None) is a TypeError), so use drop(errors="ignore").
    out.obs = out.obs.drop(columns=["_epi", "_fib"], errors="ignore")
    return out, rep


def step4_fibroblast(full_h5):
    out_fib = os.path.join(OUT, "fibroblast_integrated.h5ad")
    if os.path.exists(out_fib):
        log("STEP4 skip (fibroblast_integrated.h5ad exists)")
        return out_fib
    adata = sc.read_h5ad(full_h5)
    fib = adata[adata.obs["celltype"] == "Fibroblast"].copy()
    log(f"STEP4 fibroblast subset: {fib.shape}")
    fib, prep = _purity_filter(fib)
    log(f"  purity filter: keep {prep['n_keep']}/{prep['n_in']} (drop {prep['pct_drop']}%)")

    # re-embed fibroblasts on the integrated rep then subcluster for subtypes
    rep = adata.uns.get("integration_rep", "X_pca_harmony")
    rep = rep if rep in fib.obsm else "X_pca"
    sc.pp.neighbors(fib, use_rep=rep, n_neighbors=15, random_state=SEED)
    # subtype clustering also uses a resolution sweep + stability pick (human-like),
    # centred on FIB_SUBCL_RES; keeps the fibroblast subtyping defensible, not 1 res.
    fib_sweep = sorted({round(FIB_SUBCL_RES * m, 2) for m in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)})
    fbest_res, fbest_key = leiden_stability_pick(fib, fib_sweep, key_prefix="fib_leiden_r")
    fib.obs["fib_leiden"] = fib.obs[fbest_key].astype("category")
    fib.uns["fib_subtype_resolution"] = float(fbest_res)
    sc.tl.umap(fib, random_state=SEED)
    log(f"  fibroblast subclustering: res={fbest_res} -> {fib.obs['fib_leiden'].nunique()} clusters")

    # ---- CONSENSUS subtype call (same design as STEP3) -------------------------
    # (a) MARKER call: z-scored multi-state score matrix, argmax over reference states
    fz, fmarker_call, fpres = marker_score_matrix(fib, FIB_SUBTYPES, "fib_leiden", "fs_")
    fz.round(3).to_csv(os.path.join(OUT, "fibro_subtype_marker_scores.csv"))
    log(f"  fibro subtype panels used ({len(fpres)}): {fpres}")

    # (b) per-cluster DE audit -- DEFERRABLE (only feeds descriptive fibro_subtype_DE.csv;
    # the subtype call below uses marker-argmax + myofib signature, NOT uns['frgg']).
    if os.environ.get("SSC_SKIP_RANKGENES") == "1":
        log("  [subtype] rank_genes_groups(fib_leiden) SKIPPED (SSC_SKIP_RANKGENES=1) -> backfill via rank_genes_standalone.py")
        fde = pd.DataFrame()   # keep DEFINED (used later for subtype top_DE); blank until standalone backfills
    else:
        sc.tl.rank_genes_groups(fib, "fib_leiden", method="wilcoxon", use_raw=True,
                                n_genes=25, key_added="frgg")
        fde = pd.DataFrame(fib.uns["frgg"]["names"])
        pd.DataFrame([dict(cluster=cl, top_DE_genes=", ".join(fde[cl].tolist()[:25]))
                      for cl in fde.columns]).to_csv(
            os.path.join(OUT, "fibro_subtype_DE.csv"), index=False)

    # (c) myofibroblast = a DISCRETE subtype called by the SAME marker-argmax rule as
    #     every other fibroblast state (Tabib18/21 + Gur22 panels), NOT a self-created
    #     top-tertile signature cut. The Tabib21/Gur22 Myofibroblast panel (ACTA2/
    #     CTHRC1/COMP/POSTN/COL11A1/FN1) is one of the FIB_SUBTYPES scored above, so a
    #     subcluster is Myofibroblast iff that panel is its argmax. This keeps the
    #     subtyping uniform and literature-anchored, and avoids over-calling (the old
    #     2/3-quantile rule mechanically labelled 1/3 of clusters ~39% of fibroblasts
    #     myofibroblast, sweeping in clusters barely above the median signature). The
    #     KO / gate-relevant myofibroblast identity is the trajectory TERMINAL macro-
    #     state (STEP 04 CellRank2), i.e. the fibrotic end-point the virtual KO shifts
    #     cells toward -- the static label here just marks the atlas subtype.
    mg = [g for g in MYOFIB_SIG if g in fib.raw.var_names]
    sc.tl.score_genes(fib, mg, score_name="myofib_signature", use_raw=True,
                      random_state=SEED)
    myo_by_cl = fib.obs.groupby("fib_leiden", observed=True)["myofib_signature"].mean()
    # data-driven sanity reference (median of cluster means) -- used ONLY to flag a
    # cluster for human review when its argmax and fibrotic signature disagree; it
    # never overrides the label.
    myo_med = float(myo_by_cl.median())

    # (d) subtype = marker argmax (uniform); flag for review where the fibrotic
    #     signature and the argmax call disagree (argmax=Myofib but low signature, or
    #     argmax!=Myofib but signature above the cluster-median) so a human can look.
    frows, fassign = [], {}
    for cl in fz.index.astype(str):
        mk = fmarker_call[cl]
        sig_hi = float(myo_by_cl.loc[cl]) >= myo_med
        final = mk
        needs_review = (mk == "Myofibroblast") != sig_hi
        why = ("argmax==signature" if not needs_review
               else f"argmax({mk}) vs fibrotic-signature {'high' if sig_hi else 'low'} -- inspect")
        fassign[cl] = final
        frows.append(dict(cluster=cl, marker_call=mk, final_subtype=final,
                          myofib_signature_mean=round(float(myo_by_cl.loc[cl]), 3),
                          fibrotic_signature_high=sig_hi,
                          needs_review=needs_review, reason=why,
                          top_DE=", ".join(fde[cl].tolist()[:10]) if cl in fde.columns else ""))
    fcons = pd.DataFrame(frows).set_index("cluster")
    fcons.to_csv(os.path.join(OUT, "fibro_subtype_consensus.csv"))
    fib.obs["fib_subtype"] = fib.obs["fib_leiden"].astype(str).map(fassign).astype(str)
    fib.obs["fib_subtype_needs_review"] = fib.obs["fib_leiden"].astype(str).map(
        fcons["needs_review"].to_dict()).astype(bool)

    log("  ===== PER-CLUSTER FIBRO SUBTYPE VERDICT =====")
    for cl in fz.index.astype(str):
        r = fcons.loc[cl]
        flag = "  <-- NEEDS_REVIEW" if r["needs_review"] else ""
        log(f"    fib cluster {cl:>3}: {r['final_subtype']:<15} "
            f"[marker={r['marker_call']}, myo_sig={r['myofib_signature_mean']}, "
            f"sig_hi={r['fibrotic_signature_high']}] {r['top_DE']}{flag}")
    log("  subtype counts: " + json.dumps(fib.obs.fib_subtype.value_counts().to_dict()))

    # MYO program mean (identical genes/definition to _qc_common.MYO) for continuity
    if "myo" not in fib.obs:
        pg = [g for g in MYO if g in fib.raw.var_names]
        fib.obs["myo"] = np.asarray(fib.raw[:, pg].X.mean(1)).ravel()

    fib.write(out_fib, compression="gzip")
    log(f"STEP4 wrote {out_fib}: {fib.shape}")
    return out_fib


# ==============================================================================
def main():
    t0 = time.time()
    log("=== STEP 02: max-data integration + annotation ===")
    log(f"CONFIG ambient={AMBIENT_METHOD} doublet={DOUBLET_METHOD} "
        f"scVI={USE_SCVI} N_MADS={N_MADS} HVG={N_HVG}")
    step1_qc_per_cohort()
    integ = step2_integrate()
    full  = step3_annotate(integ)
    fibh  = step4_fibroblast(full)
    log(f"=== DONE in {(time.time()-t0)/60:.1f} min ===")
    log(f"OUTPUTS: {full}\n         {fibh}")
    log("INTEGRATE_ANNOTATE_DONE")


if __name__ == "__main__":
    main()
