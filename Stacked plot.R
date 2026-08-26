
# ============================================================
# Hive A vs Hive B - Virome stacked bar plot (Top 10 viruses)
# ============================================================
# Reads Hive1_TPM_matrix.xlsx and Hive2_TPM_matrix.xlsx, identifies
# the top 10 viruses by combined TPM across BOTH hives, and plots
# stacked bar charts with ONE shared color scheme -- so the same
# virus always gets the same color in both panels. Any virus
# outside the top 10 is grouped into "Other".

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)   # for side-by-side panels with a shared legend

# ---- 1. Load data ----
hive1 <- read_excel("D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive1_TPM_matrix.xlsx", sheet = "TPM_matrix")
hive2 <- read_excel("D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive2_TPM_matrix.xlsx", sheet = "TPM_matrix")

# Strip any accidental whitespace in virus names (e.g. "HPLV 34 " -> "HPLV 34")
hive1$NAME <- trimws(hive1$NAME)
hive2$NAME <- trimws(hive2$NAME)

# Sample columns for each hive (edit these if your column names differ)
hive1_samples <- c("C1", "C2", "C3", "E1", "A1", "A2", "A3")
hive2_samples <- c("D1", "D2", "D3", "F2", "B1", "B2", "B3")

# Nicer display labels for the x-axis
hive1_labels <- c("Honeybee-C-1", "Honeybee-C-2", "Honeybee-C-3", "Queenbee-E-1",
                  "Varroa-A-1", "Varroa-A-2", "Varroa-A-3")
hive2_labels <- c("Honeybee-D-1", "Honeybee-D-2", "Honeybee-D-3", "Queenbee-F-2",
                  "Varroa-B-1", "Varroa-B-2", "Varroa-B-3")

# ---- 2. Determine top 10 viruses by COMBINED total TPM across both hives ----
h1_totals <- hive1 %>%
  select(NAME, all_of(hive1_samples)) %>%
  mutate(total = rowSums(select(., all_of(hive1_samples)))) %>%
  select(NAME, total)

h2_totals <- hive2 %>%
  select(NAME, all_of(hive2_samples)) %>%
  mutate(total = rowSums(select(., all_of(hive2_samples)))) %>%
  select(NAME, total)

combined_totals <- full_join(h1_totals, h2_totals, by = "NAME", suffix = c("_h1", "_h2")) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
  mutate(grand_total = total_h1 + total_h2) %>%
  arrange(desc(grand_total))

top10 <- combined_totals$NAME[1:10]
print(top10)

# ---- 3. Reshape each hive to long format, collapsing non-top10 into "Other" ----
prep_long <- function(df, samples, labels) {
  df_long <- df %>%
    select(NAME, all_of(samples)) %>%
    pivot_longer(cols = all_of(samples), names_to = "sample", values_to = "TPM") %>%
    mutate(virus = if_else(NAME %in% top10, NAME, "Other")) %>%
    group_by(sample, virus) %>%
    summarise(TPM = sum(TPM), .groups = "drop")
  
  # relabel sample codes with nice display names, preserving order
  sample_map <- setNames(labels, samples)
  df_long$sample <- factor(sample_map[df_long$sample], levels = labels)
  df_long
}

hive1_long <- prep_long(hive1, hive1_samples, hive1_labels)
hive2_long <- prep_long(hive2, hive2_samples, hive2_labels)

# ---- 4. Shared color mapping (same virus = same color in both plots) ----
virus_levels <- c(top10, "Other")   # order controls stacking + legend order

color_map <- c(
  "DWV A"     = "#0b3d67",  # dark navy
  "LSV4"    = "#2e8b57",  # green
  "BeeML"   = "#f4a442",  # orange
  "VD9"     = "#d1352b",  # red
  "VD3"     = "#7d5ba6",  # purple
  "VD5"     = "#8c6d1f",  # brown/olive
  "LSV8"    = "#f2a6c1",  # pink
  "HPLV 34" = "#f4d35e",  # yellow
  "VD2"     = "#a8d8ea",  # light blue
  "IAPV"    = "#4fb3a9",  # teal
  "Other"   = "#d3d3d3"   # light gray
)
# NOTE: if your actual top10 differs from the list above, R will simply be
# missing a color for any new name -- add it to color_map before plotting.

# Stack order: lowest-abundance virus at the bottom, DWV (dominant) at the top
stack_order <- rev(virus_levels)
hive1_long$virus <- factor(hive1_long$virus, levels = stack_order)
hive2_long$virus <- factor(hive2_long$virus, levels = stack_order)

# ---- 5. Build each panel ----
make_panel <- function(df, title) {
  ggplot(df, aes(x = sample, y = TPM, fill = virus)) +
    geom_bar(stat = "identity", width = 0.65, color = "white", linewidth = 0.15) +
    scale_fill_manual(values = color_map, breaks = virus_levels, name = "Virus") +
    scale_y_continuous(labels = scales::comma) +
    coord_cartesian(ylim = c(0, 1e6)) +  # visually caps the axis WITHOUT dropping/NA-ing
    # any bar segments -- some sample totals in the
    # source data are 1,000,000.01 due to floating-point
    # rounding in the original TPM calculation; using
    # scale_y_continuous(limits=...) here would silently
    # convert those to NA and drop them with a warning
    labs(title = title, x = NULL, y = "TPM") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}

p1 <- make_panel(hive1_long, "Hive A")
p2 <- make_panel(hive2_long, "Hive B")

# ---- 6. Combine panels side by side with ONE shared legend ----
combined_plot <- (p1 + p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(combined_plot)

# ---- 7. Save ----
ggsave("D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive_A_B_virome_stacked_barplot.png", combined_plot, width = 14, height = 7, dpi = 300)
ggsave("D:/bee_immune_analysis/qc/Align/aligned/newvirus/Hive_A_B_virome_stacked_barplot.pdf", combined_plot, width = 14, height = 7)
