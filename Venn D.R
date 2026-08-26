pkgs <- c("readxl","dplyr","tidyr","ggplot2","ggforce")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(ggforce)

# =====================================================================
# PART 1: (unchanged) Venn diagram function - keep using as before
# =====================================================================
make_venn3 <- function(xlsx_path, groups, group_colors, title, out_pdf, out_png) {
  df <- read_excel(xlsx_path, sheet = 1)
  names(df) <- trimws(names(df))
  
  build_set <- function(cols) {
    sub <- df %>% select(NAME, all_of(cols)) %>%
      mutate(across(all_of(cols), as.numeric))
    present <- sub$NAME[rowSums(sub[cols] > 0, na.rm=TRUE) > 0]
    present[!is.na(present)]
  }
  
  set_list <- lapply(groups, build_set)
  names(set_list) <- names(groups)
  
  gA <- names(groups)[1]; gB <- names(groups)[2]; gC <- names(groups)[3]
  A <- set_list[[gA]]; B <- set_list[[gB]]; C <- set_list[[gC]]
  
  onlyA   <- setdiff(A, union(B, C))
  onlyB   <- setdiff(B, union(A, C))
  onlyC   <- setdiff(C, union(A, B))
  AB_only <- setdiff(intersect(A, B), C)
  AC_only <- setdiff(intersect(A, C), B)
  BC_only <- setdiff(intersect(B, C), A)
  ABC     <- intersect(intersect(A, B), C)
  
  wrap_list <- function(v) {
    if (length(v)==0) return("")
    if (length(v) <= 8) return(paste(v, collapse="\n"))
    ncol <- 2
    nrow <- ceiling(length(v)/ncol)
    v_pad <- c(v, rep("", nrow*ncol - length(v)))
    m <- matrix(v_pad, nrow=nrow, ncol=ncol)
    paste(apply(m, 1, function(row) paste(format(row, width=9), collapse="  ")), collapse="\n")
  }
  
  circles <- data.frame(
    x0 = c(-1.2, 1.2, 0),
    y0 = c(0.7, 0.7, -1.1),
    r  = c(2, 2, 2),
    group = c(gA, gB, gC)
  )
  
  p <- ggplot() +
    geom_circle(data=circles, aes(x0=x0, y0=y0, r=r, color=group, fill=group),
                alpha=0.25, linewidth=1.1) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    coord_fixed() +
    theme_void() +
    theme(legend.position="none",
          plot.title = element_text(hjust=0.5, face="bold", size=14))
  
  p <- p +
    annotate("text", x=-1.2, y=3.35, label=gA, fontface="bold", size=5, color=group_colors[[gA]]) +
    annotate("text", x=1.2,  y=3.35, label=gB, fontface="bold", size=5, color=group_colors[[gB]]) +
    annotate("text", x=0,    y=-3.75, label=gC, fontface="bold", size=5, color=group_colors[[gC]])
  
  p <- p +
    annotate("text", x=-1.2, y=3.05,  label=paste(groups[[gA]], collapse="-"), size=3.5) +
    annotate("text", x=1.2,  y=3.05,  label=paste(groups[[gB]], collapse="-"), size=3.5) +
    annotate("text", x=0,    y=-3.45, label=paste(groups[[gC]], collapse="-"), size=3.5)
  
  p <- p +
    annotate("text", x=-2.6, y=1.3,  label=length(onlyA),   color="red", fontface="bold", size=6) +
    annotate("text", x=2.6,  y=1.3,  label=length(onlyB),   color="red", fontface="bold", size=6) +
    annotate("text", x=0,    y=-2.6, label=length(onlyC),   color="red", fontface="bold", size=6) +
    annotate("text", x=0,    y=1.9,  label=length(AB_only), color="red", fontface="bold", size=6) +
    annotate("text", x=-1.5, y=-0.8, label=length(AC_only), color="red", fontface="bold", size=6) +
    annotate("text", x=1.5,  y=-0.8, label=length(BC_only), color="red", fontface="bold", size=6) +
    annotate("text", x=0,    y=0.2,  label=length(ABC),     color="red", fontface="bold", size=6)
  
  p <- p +
    annotate("text", x=-2.6, y=1.05, label=wrap_list(onlyA),   size=4, fontface="bold", hjust=0.5, vjust=1, lineheight=0.9, family="mono") +
    annotate("text", x=2.7,  y=1.05, label=wrap_list(onlyB),   size=4, fontface="bold", hjust=0.5, vjust=1, lineheight=0.9, family="mono") +
    annotate("text", x=0,    y=-2.85,label=wrap_list(onlyC),   size=4, fontface="bold", hjust=0.5, vjust=1, lineheight=0.9, family="mono") +
    annotate("text", x=0,    y=1.65, label=wrap_list(AB_only), size=4, fontface="bold", hjust=0.5, vjust=1, lineheight=0.9, family="mono") +
    annotate("text", x=-1.5, y=-0.95,label=wrap_list(AC_only), size=4, fontface="bold", hjust=0.5, vjust=1, lineheight=0.9, family="mono") +
    annotate("text", x=1.5,  y=-0.95,label=wrap_list(BC_only), size=4, fontface="bold", hjust=0.5, vjust=1, lineheight=0.9, family="mono") +
    annotate("text", x=0,    y=0.05, label=wrap_list(ABC),     size=4, fontface="bold", hjust=0.5, vjust=1, lineheight=0.9, family="mono") +
    labs(title=title) +
    xlim(-4.5, 4.5) + ylim(-4.1, 3.7)
  
  ggsave(out_pdf, p, width=11, height=11)
  ggsave(out_png, p, width=11, height=11, dpi=300)
  
  return(p)
}


# ── Hive 1 ─────────────────────────────────────────────────────
hive1_groups <- list(Honeybee=c("C1","C2","C3"), Varroa=c("A1","A2","A3"), `Queen Bee`=c("E1"))
hive1_all_cols  <- unlist(hive1_groups, use.names = FALSE)                                  # Venn + group-count order (unchanged)
hive1_bar_order <- c(hive1_groups$Honeybee, hive1_groups$`Queen Bee`, hive1_groups$Varroa)  # bar-plot only: Queen Bee in the middle

make_venn3(
  xlsx_path = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive1_TPM_matrix.xlsx",
  groups = hive1_groups,
  group_colors = group_colors,
  title = "Virus Distribution across Host Groups (Hive 1)",
  out_pdf = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive1_Venn.pdf",
  out_png = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive1_Venn.png"
)


# ── Hive 2 ─────────────────────────────────────────────────────
hive2_groups <- list(Honeybee=c("D1","D2","D3"), Varroa=c("B1","B2","B3"), `Queen Bee`=c("F2"))
hive2_all_cols  <- unlist(hive2_groups, use.names = FALSE)                                  # Venn + group-count order (unchanged)
hive2_bar_order <- c(hive2_groups$Honeybee, hive2_groups$`Queen Bee`, hive2_groups$Varroa)  # bar-plot only: Queen Bee in the middle

make_venn3(
  xlsx_path = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive2_TPM_matrix.xlsx",
  groups = hive2_groups,
  group_colors = group_colors,
  title = "Virus Distribution across Host Groups (Hive 2)",
  out_pdf = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive2_Venn.pdf",
  out_png = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive2_Venn.png"
)









