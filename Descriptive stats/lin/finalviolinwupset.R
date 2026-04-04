library(dplyr)
library(ggplot2)
library(readr)
library(wesanderson)
library(patchwork)

# install.packages("ComplexUpset")  # once
library(ComplexUpset)


snps <- read_csv("SignificantGenomewideSNPs_withEffects.csv", show_col_types = FALSE)

metabolic_traits <- c(
  "max_DeltaE_5","auc_5","slope_5","r_squared_5",
  "max_b","tp_max_b","auc_b","slope_b_5","r2_b_5","range_b"
)

snps <- snps %>%
  mutate(
    Trait_class = if_else(Trait %in% metabolic_traits, "Metabolic", "Structural"),
    Trait_class = factor(Trait_class, levels = c("Structural","Metabolic")),
    abs_Efsize  = abs(Efsize)
  )

# ---- Colors ----
pal <- wes_palette("IsleofDogs1", 2)
fill_vals <- c(Structural = pal[1], Metabolic = pal[2])

# -----------------------------
# A) Violin (all associations/hits)
# -----------------------------
lab_df <- snps %>%
  group_by(Trait_class) %>%
  summarise(hits = n(), snps = n_distinct(ID), .groups = "drop")

x_labs <- setNames(
  paste0(lab_df$Trait_class, "\n(hits=", lab_df$hits, ", SNPs=", lab_df$snps, ")"),
  lab_df$Trait_class
)

p_violin <- ggplot(snps, aes(x = Trait_class, y = abs_Efsize, fill = Trait_class)) +
  geom_violin(trim = TRUE, bw = 0.01, color = "black", linewidth = 0.6) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.6, fill = "white") +
  stat_summary(fun = median, geom = "point", size = 3, shape = 21,
               fill = "white", color = "black") +
  scale_fill_manual(values = fill_vals) +
  scale_x_discrete(labels = x_labs) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Standardized SNP effect sizes", x = NULL, y = "Standardized SNP effect size") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(8, 8, 8, 8)
  )

# -----------------------------
# B) UpSet (UNIQUE SNPs membership by class)
# -----------------------------
snp_sets <- snps %>%
  distinct(ID, Trait_class) %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(
    names_from  = Trait_class,
    values_from = present,
    values_fill = FALSE
  )

# Make sure columns exist
if (!("Structural" %in% names(snp_sets))) snp_sets$Structural <- FALSE
if (!("Metabolic"  %in% names(snp_sets))) snp_sets$Metabolic  <- FALSE
p_upset <- ComplexUpset::upset(
  snp_sets,
  intersect = c("Structural", "Metabolic"),
  name = "SNPs",
  width_ratio = 0.22,
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Intersection size" = intersection_size(
      counts = TRUE,
      text = list(size = 4)
    )
  ),
  set_sizes = (
    upset_set_size() +
      aes(fill = after_stat(set)) +
      scale_fill_manual(values = fill_vals) +
      guides(fill = "none")
  )
) +
  labs(title = "Significant SNPs by trait class") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )
# -----------------------------
# Combine: wrap UpSet so it behaves like ONE panel (prevents C/D tags)
# -----------------------------
p_combo <- (p_violin | patchwork::wrap_elements(full = p_upset)) +
  plot_layout(widths = c(1.05, 1.35)) +
  plot_annotation(tag_levels = "A")  # now you should only get A and B

p_combo

# -----------------------------
# Export (Illustrator-friendly PDF)
# -----------------------------
ggsave("Violin_plus_UpSet_clean.pdf", p_combo,
       width = 10.5, height = 5.2, units = "in",
       device = cairo_pdf, limitsize = FALSE)
