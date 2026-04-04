library(dplyr)
library(ggplot2)
library(readr)

snps <- read_csv("SignificantGenomewideSNPs_withEffects.csv",
                 show_col_types = FALSE)

metabolic_traits <- c(
  "max_DeltaE_5","auc_5","slope_5",
  "max_b","tp_max_b","auc_b","slope_b_5","range_b"
)

snps <- snps %>%
  mutate(
    Trait_class = if_else(Trait %in% metabolic_traits,
                          "Metabolic", "Structural")
  )

snps <- snps %>%
  mutate(abs_Efsize = abs(Efsize))

ggplot(snps, aes(Trait_class, abs_Efsize, fill = Trait_class)) +
  geom_violin(
    trim = TRUE,
    bw = 0.01,
    color = "black",
    linewidth = 0.6
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    linewidth = 0.6,
    fill = "white"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 3,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Structural" = "#0072B2",  # Okabe–Ito blue
      "Metabolic"  = "#D55E00"   # Okabe–Ito vermillion
    )
  ) +
  labs(
    x = NULL,
    y = "Standardized SNP effect size",
    #title = "Distribution of SNP effect sizes",
    #subtitle = "Structural vs metabolic traits"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

#standard deviation of snps
snps %>%
  group_by(Trait_class) %>%
  summarize(
    mean_effect = mean(abs_Efsize),
    sd_effect   = sd(abs_Efsize),
    median_effect = median(abs_Efsize),
    n = n(),
    .groups = "drop"
  )

# 1) Number of significant SNP–trait associations (row count)
snps %>% count(Trait_class, name = "n_assoc")

# 2) Number of unique significant SNP IDs per class
snps %>%
  distinct(ID, Trait_class) %>%
  count(Trait_class, name = "n_unique_snps")

# 3) Number of traits contributing signals per class (optional)
snps %>%
  group_by(Trait_class) %>%
  summarize(n_traits = n_distinct(Trait), .groups = "drop")

snps %>%
  group_by(Trait_class) %>%
  summarize(
    median_abs_effect = median(abs_Efsize, na.rm = TRUE),
    mean_abs_effect   = mean(abs_Efsize, na.rm = TRUE),
    sd_abs_effect     = sd(abs_Efsize, na.rm = TRUE),
    n_assoc           = n(),
    .groups = "drop"
  )
