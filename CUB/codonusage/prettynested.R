library(dplyr)
library(ggplot2)
library(scales)

# Okabe–Ito / your palette
cb <- c(
  non_coding  = "#999999",
  coding_total = "#000000",
  synonymous  = "#0072B2",
  missense    = "#D55E00",
  nonsense    = "#E69F00",
  stop_loss   = "#009E73",
  other       = "#56B4E9",
  unknown     = "#CC79A7"
)

# -----------------------------
# Build "nested" plotting data
# -----------------------------
# Broad (All SNPs)
broad_plot <- gwas_snp_summary %>%
  filter(Category_Level == "broad") %>%
  mutate(
    class = if_else(SNP_Class == "coding", "coding", "non_coding"),
    group = "All significant SNPs"
  ) %>%
  group_by(group) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

# Coding subclasses
coding_plot <- gwas_snp_summary %>%
  filter(Category_Level == "coding") %>%
  mutate(
    class = SNP_Class,
    group = "Coding SNPs only"
  ) %>%
  group_by(group) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

plot_df <- bind_rows(broad_plot, coding_plot)

# Define a cleaner hierarchical palette
cb2 <- c(
  non_coding = "#999999",
  coding     = "#1B3A4B",     # dark blue (parent)
  synonymous = "#2C7FB8",
  missense   = "#41B6C4",
  nonsense   = "#7FCDBB",
  stop_loss  = "#C7E9B4",
  other      = "#EDF8B1",
  unknown    = "#CCCCCC"
)

p_nested <- ggplot(plot_df, aes(x = group, y = n, fill = class)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = cb2) +
  theme_classic(base_size = 13) +
  labs(
    title = "Genome-Wide Significant SNPs",
    x = NULL,
    y = "Count",
    fill = "Annotation class"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 11)
  )

p_nested
ggsave("Fig_SNP_annotation_nested.pdf", plot = p_nested,
       width = 6.5, height = 4.5, units = "in", dpi = 600)
