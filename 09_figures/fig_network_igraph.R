# ==========================================================================
# fig_network_igraph.R
# 用途: HES1/FOSB 调控网络图 — 基于库 047_RcisTarget igraph 方案改编
# 布局: Kamada-Kawai (layout_with_kk), fibrotic effector 高亮标注
# 输出: Fig_regulon_network.pdf + Fig_regulon_network.png
# ==========================================================================
library(igraph)

# ---------- 数据 ----------
net <- read.csv("C:/Users/fsy/Desktop/SSc_virtualKO_MR/server_results/data/collectri_net.csv",
                stringsAsFactors = FALSE)

hubs    <- c("HES1", "FOSB")
net_sub <- net[net$source %in% hubs, ]

FIBRO <- c(
  "COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL11A1","COL6A3",
  "ACTA2","TAGLN","FN1","POSTN","CTHRC1","COMP","TIMP1","CCN2","CTGF",
  "SERPINE1","LOX","LOXL1","LOXL2","SPARC","THBS1","THBS2","FBN1",
  "ELN","MMP2","MMP9","MMP14","VIM","SNAI1","SNAI2","TWIST1","ZEB1",
  "TNC","PDGFRB","TGFB1","TGFBR1","SMAD7","VCAN","BGN","DCN","LUM",
  "PLOD2","P4HA1","CDKN1A","VEGFA","HEY1","MYC","JUN","FOS","IL6",
  "MMP1","MMP3"
)

# ---------- 构图 ----------
all_nodes <- unique(c(net_sub$source, net_sub$target))
g <- graph_from_data_frame(
  d        = data.frame(from=net_sub$source, to=net_sub$target,
                         weight=net_sub$weight),
  directed = TRUE,
  vertices = data.frame(name=all_nodes)
)

is_hub   <- V(g)$name %in% hubs
is_fibro <- V(g)$name %in% FIBRO & !is_hub

# 共享靶基因: 被 HES1 和 FOSB 同时靶向
tgt_counts <- table(net_sub$target)
is_shared  <- V(g)$name %in% names(tgt_counts[tgt_counts >= 2]) & is_fibro

# 边类型
src <- ends(g, E(g))[, 1]
tgt <- ends(g, E(g))[, 2]
edge_to_fibro <- tgt %in% FIBRO

# ---------- 节点样式 ----------
V(g)$size <- ifelse(is_hub, 38,
             ifelse(is_fibro, 13, 4))

V(g)$color <- ifelse(V(g)$name == "HES1", "#2166AC",
              ifelse(V(g)$name == "FOSB", "#D6604D",
              ifelse(is_shared,  "#8E44AD",
              ifelse(is_fibro,   "#B2182B",
                                 "#D9D9D9"))))

V(g)$frame.color <- ifelse(is_hub, "#1A1A1A",
                    ifelse(is_fibro, "#7B0000", "#AAAAAA"))
V(g)$frame.width <- ifelse(is_hub, 2.5, ifelse(is_fibro, 1.2, 0.3))

V(g)$label <- ifelse(is_hub | is_fibro, V(g)$name, "")
V(g)$label.cex   <- ifelse(is_hub, 1.15, 0.72)
V(g)$label.color <- ifelse(is_hub, "white", "#111111")
V(g)$label.font  <- ifelse(is_hub, 2L, 2L)  # bold

# ---------- 边样式 ----------
E(g)$color <- ifelse(E(g)$weight < 0, "#F4A582", "#BDD7E7")
E(g)$width <- ifelse(edge_to_fibro,
                     ifelse(E(g)$weight < 0, 1.8, 1.5),
                     0.25)
E(g)$arrow.size <- 0

# ---------- 布局 ----------
set.seed(42)
# Pin hub positions: HES1 left, FOSB right; let targets spread freely via FR
n  <- vcount(g)
minx <- maxx <- miny <- maxy <- rep(NA_real_, n)
minx[V(g)$name == "HES1"] <- maxx[V(g)$name == "HES1"] <- -5.0
miny[V(g)$name == "HES1"] <- maxy[V(g)$name == "HES1"] <-  0.0
minx[V(g)$name == "FOSB"] <- maxx[V(g)$name == "FOSB"] <-  5.0
miny[V(g)$name == "FOSB"] <- maxy[V(g)$name == "FOSB"] <-  0.0
l <- layout_with_fr(g,
                    weights = abs(E(g)$weight) + 0.1,
                    niter   = 2000,
                    minx = minx, maxx = maxx,
                    miny = miny, maxy = maxy)

# ---------- 输出目录 ----------
outdir <- "C:/Users/fsy/Desktop/SSc_virtualKO_MR/figures"

# ---------- 绘图 ----------
draw_net <- function() {
  par(bg = "white", mar = c(4, 1, 3, 1))
  plot(g,
       layout            = l,
       vertex.size       = V(g)$size,
       vertex.color      = V(g)$color,
       vertex.frame.color= V(g)$frame.color,
       vertex.label      = V(g)$label,
       vertex.label.cex  = V(g)$label.cex,
       vertex.label.color= V(g)$label.color,
       vertex.label.font = V(g)$label.font,
       edge.color        = E(g)$color,
       edge.width        = E(g)$width,
       edge.arrow.size   = 0,
       main              = "")
  title(main = "HES1 and FOSB regulons converge on fibrotic effector genes",
        cex.main = 1.25, font.main = 2, line = 0.5)

  # Legend
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.82,
         pt.cex = c(2.2, 2.2, 1.5, 1.0),
         pch    = 21,
         pt.bg  = c("#2166AC", "#D6604D", "#B2182B", "#D9D9D9"),
         col    = c("#1A1A1A", "#1A1A1A", "#7B0000", "#AAAAAA"),
         legend = c("HES1 (Notch hub)", "FOSB (AP-1 hub)",
                    "fibrotic effector target", "other regulon target"),
         inset  = 0.01)

  # Edge legend
  legend("bottomleft", bty = "n", cex = 0.78,
         lty = 1, lwd = c(2, 2),
         col = c("#BDD7E7", "#F4A582"),
         legend = c("activation", "repression"),
         inset  = 0.01)
}

# PDF
pdf(file.path(outdir, "Fig_regulon_network.pdf"), width = 14, height = 9.5)
draw_net()
dev.off()

png(file.path(outdir, "Fig_regulon_network.png"),
    width = 14 * 400, height = 9.5 * 400, res = 400)
draw_net()
dev.off()

cat("saved Fig_regulon_network.pdf / .png\n")
cat(sprintf("nodes=%d  hub=%d  fibro=%d  other=%d\n",
            vcount(g), sum(is_hub), sum(is_fibro), sum(!is_hub & !is_fibro)))
