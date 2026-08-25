# Environments

Several tools have conflicting dependencies — use isolated environments.

## 1. single-cell core (conda/mamba, Python 3.10)
scanpy 1.12, anndata 0.12, celloracle, decoupler 2.1, cellrank 2, squidpy, scikit-learn, networkx

## 2. SCENIC+ (separate env, Python 3.11)
scenicplus, pycisTopic 2.0a0, pycistarget, snakemake, Mallet (LDA)

## 3. R (>= 4.3)
Seurat, nichenetr, liana, IOBR, spacexr (RCTD), Banksy, circlize, ComplexHeatmap, ggraph, igraph

## 4. molecular dynamics (separate env)
GROMACS (CUDA build), gmx_MMPBSA 1.6.5 (conda mpi4py + pip), CHARMM-GUI inputs, AutoDock Vina

Genome: GRCh38. Priors: CollecTRI (via decoupler), NicheNet prior model.
Exact versions are pinned in each stage's install script where present.
