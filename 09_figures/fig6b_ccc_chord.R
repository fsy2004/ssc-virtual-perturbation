# Fig 6 (advanced) — Notch niche -> myofibroblast chord diagram (circlize).
# Perivascular/niche senders signal via Notch ligands to myofibroblast NOTCH receptors.
# Data: LIANA SSc Notch niche->myofibroblast (verified). No bar charts. Vector PDF.
suppressMessages({library(circlize); library(grid)})
DL <- "C:/Users/fsy/Desktop/SSc_virtualKO_MR/04_manuscript/plot_data_local/ccc_nichenet/liana"
OUT <- "C:/Users/fsy/Desktop/SSc_virtualKO_MR/03_results_figures/figures/Fig6b_ccc_chord"

d <- read.csv(file.path(DL, "liana_SSc_notch_niche2myo.csv"), stringsAsFactors = FALSE)
d <- d[order(-d$lr_means), ]
# keep informative senders; build sender -> receptor edges weighted by lr_means, coloured by ligand
d$sender  <- d$source
d$recept  <- paste0("Myofib:", d$receptor_complex)   # NOTCH2 / NOTCH3 on myofibroblast
d$ligand  <- d$ligand_complex

senders <- unique(d$sender)
recepts <- unique(d$recept)
ligs    <- unique(d$ligand)

# colour-blind-safe palette
sender_col <- setNames(c("#4C72B0","#DD8452","#55A868","#8172B3","#937860","#DA8BC3")[seq_along(senders)], senders)
recept_col <- setNames(c("#C44E52","#B2182B")[seq_along(recepts)], recepts)
grid_col   <- c(sender_col, recept_col)
lig_col    <- setNames(c("#2166AC","#4C72B0","#8C8C8C","#CCB974")[seq_along(ligs)], ligs)  # DLL4/JAG1 emphasised

# edge list: sender -> receptor, width = lr_means, ribbon colour = ligand
edges <- data.frame(from = d$sender, to = d$recept, value = d$lr_means,
                    col = lig_col[d$ligand], stringsAsFactors = FALSE)

for (dev in c("pdf","png")) {
  if (dev=="pdf") pdf(paste0(OUT,".pdf"), width=7.4, height=7.4, useDingbats=FALSE)
  else png(paste0(OUT,".png"), width=1750, height=1750, res=230)
  par(mar = c(1,1,2,1))
  circos.clear()
  circos.par(canvas.xlim = c(-1.45, 1.45), canvas.ylim = c(-1.45, 1.45),
             gap.after = c(rep(3, length(senders)-1), 12, rep(3, length(recepts)-1), 12),
             start.degree = 90)
  chordDiagram(edges[,c("from","to","value")],
               grid.col = grid_col, col = edges$col,
               directional = 1, direction.type = c("diffHeight","arrows"),
               link.arr.type = "big.arrow", diffHeight = mm_h(2),
               annotationTrack = "grid", preAllocateTracks = list(track.height = 0.09),
               link.sort = TRUE, link.decreasing = TRUE)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    s <- get.cell.meta.data("sector.index")
    circos.text(mean(get.cell.meta.data("xlim")), get.cell.meta.data("ylim")[1] + 0.9,
                gsub("Myofib:","", s), facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.85, font = 2)
  }, bg.border = NA)
  # ligand legend
  legend("topright", legend = names(lig_col), col = lig_col, lwd = 4, bty = "n",
         cex = 0.8, title = "Notch ligand", title.font = 2)
  title(main = "Niche -> myofibroblast Notch signalling (SSc)", cex.main = 1.05)
  dev.off()
}
cat("chord saved:", OUT, "\nsenders:", paste(senders, collapse=", "),
    "\nreceptors:", paste(recepts, collapse=", "),
    "\nligands:", paste(ligs, collapse=", "), "\n")
