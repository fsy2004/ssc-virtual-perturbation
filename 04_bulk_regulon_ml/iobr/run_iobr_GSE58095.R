# auto-generated IOBR deconvolution (module 492) for GSE58095  [SEED 123]
suppressMessages({library(IOBR); library(tidyverse)})
set.seed(123)
eset <- as.matrix(read.csv("./out_05_bulk/stageB/GSE58095_expr_for_deconv.csv", row.names=1, check.names=FALSE))
# IOBR expects genes(rows, symbol) x samples(cols); ★线性 TPM/FPKM(已在 _deconvolve 反 log)。
methods <- c("mcpcounter","xcell","epic","quantiseq","timer","cibersort","estimate")
decon <- lapply(methods, function(m)
  tryCatch(deconvo_tme(eset=eset, method=m, arrays=TRUE, perm=1000),
           error=function(e){message(m," FAILED: ",conditionMessage(e)); NULL}))
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
write.csv(tme, "./out_05_bulk/stageB/GSE58095_deconv.csv")
cat("IOBR done:", ncol(tme), "features x", nrow(tme), "samples\n")
