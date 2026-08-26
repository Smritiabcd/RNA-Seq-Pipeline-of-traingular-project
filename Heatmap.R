pkgs <- c("readxl","ggplot2","dplyr","tidyr","scales","patchwork")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(readxl); library(ggplot2); library(dplyr)
library(tidyr);  library(scales);  library(patchwork)

make_virus_heatmap <- function(xlsx_path, out_pdf, out_png, plot_title) {
  
  # ── 1. Read ──────────────────────────────────────────────────
  df_raw <- read_excel(xlsx_path, sheet = 1)
  cat("Columns:", paste(names(df_raw), collapse=", "), "\n")
  
  cn <- tolower(trimws(names(df_raw)))
  get_col <- function(pat){
    h <- names(df_raw)[grepl(pat, cn, ignore.case=TRUE)]
    if (!length(h)) stop(paste("No column matching:", pat))
    h[1]
  }
  id_col     <- get_col("id|accession")
  name_col   <- get_col("^name|virus.?name")
  family_col <- get_col("famil")
  meta_cols  <- c(id_col, name_col, family_col)
  tpm_cols   <- setdiff(names(df_raw), meta_cols)
  tpm_cols   <- tpm_cols[!grepl("coverage|length|reads|mapped", tpm_cols, ignore.case=TRUE)]
  
  # ── 2. Sample rename ─────────────────────────────────────────
  sample_map <- c(
    "C1"="Honeybee 1","C2"="Honeybee 2","C3"="Honeybee 3",
    "D1"="Honeybee 1","D2"="Honeybee 2","D3"="Honeybee 3",
    "E1"="Queen Bee",  "F2"="Queen Bee",
    "A1"="Varroa Mite 1","A2"="Varroa Mite 2","A3"="Varroa Mite 3",
    "B1"="Varroa Mite 1","B2"="Varroa Mite 2","B3"="Varroa Mite 3"
  )
  rs <- function(x){ k<-trimws(x); if(k %in% names(sample_map)) sample_map[[k]] else k }
  
  # ── 3. Tidy ───────────────────────────────────────────────────
  df <- df_raw %>%
    rename(Accession=all_of(id_col), Virus_Name=all_of(name_col), Family=all_of(family_col)) %>%
    mutate(
      # Manual fix for MH243376.1 which is missing NAME/FAMILY in the source file
      Virus_Name = case_when(
        Accession == "MH243376.1" ~ "AMFV",
        TRUE ~ Virus_Name
      ),
      Family = case_when(
        Accession == "MH243376.1" ~ "Others",   # change if AMFV belongs to a different family
        TRUE ~ Family
      )
    ) %>%
    select(Accession, Virus_Name, Family, all_of(tpm_cols)) %>%
    mutate(across(all_of(tpm_cols), as.numeric),
           Family = ifelse(is.na(Family)|Family=="","Others",Family),
           # Guard against any other blank/NA virus names causing row misalignment
           Virus_Name = ifelse(is.na(Virus_Name) | trimws(Virus_Name)=="",
                               Accession, trimws(Virus_Name)),
           # Guard against duplicate names doing the same
           Virus_Name = make.unique(Virus_Name))
  
  virus_order_top_to_bottom <- df$Virus_Name
  df$Virus_Name <- factor(df$Virus_Name, levels = rev(virus_order_top_to_bottom))
  
  df_long <- df %>%
    pivot_longer(all_of(tpm_cols), names_to="Sample", values_to="TPM") %>%
    mutate(Sample  = sapply(Sample, rs),
           TPM_log = ifelse(is.na(TPM)|TPM<=0, NA_real_, log10(TPM)))
  
  # ── 4. Sample order ───────────────────────────────────────────
  so  <- unique(df_long$Sample)
  hb  <- sort(so[grepl("Honeybee",  so, ignore.case=TRUE)])
  qb  <- sort(so[grepl("Queen",     so, ignore.case=TRUE)])
  vm  <- sort(so[grepl("Varroa",    so, ignore.case=TRUE)])
  sample_order <- c(hb, qb, vm, setdiff(so, c(hb,qb,vm)))
  df_long$Sample <- factor(df_long$Sample, levels=sample_order)
  n_samples <- length(sample_order)
  n_viruses <- nlevels(df_long$Virus_Name)
  
  # ── 5. Family info ────────────────────────────────────────────
  family_df <- df %>%
    select(Virus_Name, Family) %>%
    mutate(y_int = as.integer(Virus_Name))
  
  fam_groups <- family_df %>%
    group_by(Family) %>%
    summarise(y_min=min(y_int), y_max=max(y_int),
              y_mid=(min(y_int)+max(y_int))/2, n=n(), .groups="drop") %>%
    arrange(y_mid)
  
  # ── 6. Color palette ─────────────────────────────────────────
  pal <- c("#C0392B","#2471A3","#1E8449","#7D3C98","#CA6F1E",
           "#117A65","#943126","#1A5276","#196F3D","#6C3483",
           "#784212","#0E6655","#884EA0","#1F618D","#B7950B")
  fams_u     <- fam_groups$Family
  fam_colors <- setNames(pal[(seq_along(fams_u)-1) %% length(pal)+1], fams_u)
  
  # ── 7. Gradient ───────────────────────────────────────────────
  grad_cols  <- c("#FFFAEE","#FDE8C8","#FBCB8A","#F5A623","#E2611A","#B22222","#6B0000")
  tpm_breaks <- 0:6
  
  # ── 8. MERGED PANEL: family brackets + virus names together (no gap) ──
  bracket_data <- fam_groups %>% mutate(col = fam_colors[Family])
  
  seg_vert <- bracket_data %>% transmute(x=1, xend=1,    y=y_min-0.4, yend=y_max+0.4, col=col)
  seg_top  <- bracket_data %>% transmute(x=1, xend=1.15, y=y_max+0.4, yend=y_max+0.4, col=col)
  seg_bot  <- bracket_data %>% transmute(x=1, xend=1.15, y=y_min-0.4, yend=y_min-0.4, col=col)
  lab_data <- bracket_data %>% transmute(x=0.97, y=y_mid, label=Family, col=col)
  
  name_data <- data.frame(
    y     = 1:n_viruses,
    label = levels(df$Virus_Name)
  )
  
  p_left <- ggplot() +
    geom_segment(data=seg_vert, aes(x=x,xend=xend,y=y,yend=yend,color=I(col)), linewidth=0.9) +
    geom_segment(data=seg_top,  aes(x=x,xend=xend,y=y,yend=yend,color=I(col)), linewidth=0.9) +
    geom_segment(data=seg_bot,  aes(x=x,xend=xend,y=y,yend=yend,color=I(col)), linewidth=0.9) +
    geom_text(data=lab_data, aes(x=x,y=y,label=label,color=I(col)),
              hjust=1, size=3.3, fontface="bold") +
    geom_text(data=name_data, aes(x=1.25, y=y, label=label),
              hjust=0, size=3.8, fontface="bold", color="black") +
    scale_y_continuous(limits=c(0.5, n_viruses+0.5), expand=c(0,0)) +
    scale_x_continuous(limits=c(0, 1.6), expand=c(0,0)) +
    theme_void() +
    theme(plot.margin = margin(t=50, r=0, b=10, l=2))
  
  # ── 9. RIGHT PANEL: heatmap with top axis only ───────────────
  p_heat <- ggplot(df_long, aes(x=Sample, y=Virus_Name, fill=TPM_log)) +
    geom_tile(color="white", linewidth=0.35) +
    scale_fill_gradientn(
      colors   = grad_cols,
      values   = rescale(tpm_breaks),
      limits   = c(0,6),
      na.value = "#F5F5F0",
      name     = "TPM",
      breaks   = tpm_breaks,
      labels   = expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6),
      guide    = guide_colorbar(barwidth=1, barheight=8)
    ) +
    scale_x_discrete(expand=c(0,0), position="top") +
    scale_y_discrete(expand=c(0,0)) +
    labs(x=NULL, y=NULL) +
    theme_minimal(base_size=10) +
    theme(
      axis.text.x      = element_text(face="bold", size=11, angle=0, hjust=0.5, vjust=0),
      axis.text.x.top  = element_text(face="bold", size=11, angle=0, hjust=0.5, vjust=0),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      panel.grid       = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(face="bold", size=9),
      plot.margin      = margin(t=40, r=10, b=10, l=0)
    )
  
  # ── 10. Combine ───────────────────────────────────────────────
  p_final <- p_left + p_heat +
    plot_layout(widths = c(1.5, 6)) +
    plot_annotation(
      title = plot_title,
      theme = theme(plot.title = element_text(face="bold", size=13, hjust=0.5))
    )
  
  # ── 11. Save ─────────────────────────────────────────────────
  fig_h <- max(12, n_viruses * 0.30 + 4)
  fig_w <- max(12, n_samples * 1.5  + 5)
  
  ggsave(out_pdf, p_final, width=fig_w, height=fig_h, device="pdf")
  ggsave(out_png, p_final, width=fig_w, height=fig_h, dpi=300)
  
  cat("\n✅ Done:", plot_title, "\n")
  cat("PDF ->", out_pdf, "\n")
  cat("PNG ->", out_png, "\n")
  
  return(p_final)
}

# ── Run for Hive 1 ─────────────────────────────────────────────
make_virus_heatmap(
  xlsx_path = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive1_TPM_matrix.xlsx",
  out_pdf   = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive1_Virus_Heatmap.pdf",
  out_png   = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive1_Virus_Heatmap.png",
  plot_title = "Hive 1 – Virus Abundance Heatmap"
)

# ── Run for Hive 2 ─────────────────────────────────────────────
make_virus_heatmap(
  xlsx_path = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive2_TPM_matrix.xlsx",
  out_pdf   = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive2_Virus_Heatmap.pdf",
  out_png   = "D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive2_Virus_Heatmap.png",
  plot_title = "Hive 2 – Virus Abundance Heatmap"
)

