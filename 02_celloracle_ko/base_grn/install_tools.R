options(Ncpus=6, repos="https://cloud.r-project.org")
ok <- function(p) requireNamespace(p, quietly=TRUE)
if(!ok("BiocManager")) install.packages("BiocManager")
if(!ok("remotes"))     install.packages("remotes")
if(!ok("Signac"))      install.packages("Signac")
BiocManager::install(c("BiocGenerics","DelayedArray","DelayedMatrixStats","limma","lme4",
  "S4Vectors","SingleCellExperiment","SummarizedExperiment","batchelor","HDF5Array",
  "terra","ggrastr","GenomicRanges","Rsamtools","rtracklayer"), update=FALSE, ask=FALSE)
if(!ok("monocle3")) remotes::install_github("cole-trapnell-lab/monocle3", upgrade="never")
if(!ok("cicero"))   remotes::install_github("cole-trapnell-lab/cicero-release", ref="monocle3", upgrade="never")
cat("INSTALL_DONE Signac=",ok("Signac")," monocle3=",ok("monocle3")," cicero=",ok("cicero"),"\n")
