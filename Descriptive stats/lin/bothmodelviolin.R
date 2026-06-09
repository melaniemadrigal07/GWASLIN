 library(dplyr)
library(ggplot2)
library(readr)
library(wesanderson)
library(patchwork)
library(ComplexUpset)
library(tidyr)

snps <- read_csv("SignificantGenomewideSNPS.csv", show_col_types = FALSE)

metabolic_traits <- c(
  "max_DeltaE_5","auc_5","slope_5","r_squared_5",
  "max_b","tp_max_b","auc_b","slope_b_5","r2_b_5","range_b"
)

snps <- snps %>%
  mutate(
    Trait_class = if_else(Trait %in% metabolic_traits, "Metabolic", "Structural"),
    Trait_class = factor(Trait_class, levels = c("Structural", "Metabolic")),
    Model = factor(Model, levels = c("K", "PK")),
    abs_Efsize = abs(EFsize)
  )

pal <- wes_palette("IsleofDogs1", 2)

fill_vals <- c(
  Structural = pal[1],
  Metabolic  = pal[2]
)

make_model_panel <- function(df, model_name) {
  
  df_model <- df %>%
    filter(Model == model_name)
  
  lab_df <- df_model %>%
    group_by(Trait_class) %>%
    summarise(
      hits = n(),
      snps = n_distinct(ID),
      .groups = "drop"
    )
  
  x_labs <- setNames(
    paste0(
      lab_df$Trait_class,
      "\n(hits=", lab_df$hits,
      ", SNPs=", lab_df$snps,
      ")"
    ),
    lab_df$Trait_class
  )
  
  p_violin <- ggplot(df_model, aes(x = Trait_class, y = abs_Efsize, fill = Trait_class)) +
    geom_violin(trim = TRUE, bw = 0.01, color = "black", linewidth = 0.6) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.6, fill = "white") +
    stat_summary(
      fun = median,
      geom = "point",
      size = 3,
      shape = 21,
      fill = "white",
      color = "black"
    ) +
    scale_fill_manual(values = fill_vals) +
    scale_x_discrete(labels = x_labs) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0(model_name, " model: standardized SNP effect sizes"),
      x = NULL,
      y = "Standardized SNP effect size"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )
  
  snp_sets <- df_model %>%
    distinct(ID, Trait_class) %>%
    mutate(present = TRUE) %>%
    pivot_wider(
      names_from = Trait_class,
      values_from = present,
      values_fill = FALSE
    )
  
  if (!("Structural" %in% names(snp_sets))) snp_sets$Structural <- FALSE
  if (!("Metabolic" %in% names(snp_sets))) snp_sets$Metabolic <- FALSE
  
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
    set_sizes = upset_set_size()
  ) +
    labs(title = paste0(model_name, " model: SNPs by trait class")) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
  
  p_violin | patchwork::wrap_elements(full = p_upset)
}

p_K <- make_model_panel(snps, "K")
p_PK <- make_model_panel(snps, "PK")

p_combo_models <- (p_K / p_PK) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

p_combo_models

ggsave(
  "Violin_plus_UpSet_by_model.pdf",
  p_combo_models,
  width = 11,
  height = 10,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)


# ---------- Median effect sizes ----------

median_summary <- snps %>%
  group_by(Model, Trait_class) %>%
  summarise(
    n_hits = n(),
    n_unique_snps = n_distinct(ID),
    median_effect = median(abs_Efsize, na.rm = TRUE),
    mean_effect = mean(abs_Efsize, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Median standardized SNP effect sizes ===\n")
print(median_summary)

# ---------- Wilcoxon rank-sum tests ----------

wilcox_results <- snps %>%
  group_by(Model) %>%
  summarise(
    p_value = wilcox.test(
      abs_Efsize ~ Trait_class
    )$p.value,
    
    structural_median = median(
      abs_Efsize[Trait_class == "Structural"],
      na.rm = TRUE
    ),
    
    metabolic_median = median(
      abs_Efsize[Trait_class == "Metabolic"],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

cat("\n=== Wilcoxon rank-sum tests ===\n")
print(wilcox_results)
