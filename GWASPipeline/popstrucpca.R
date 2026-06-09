library(data.table)
library(ggplot2)
library(patchwork)

# Read PCA + metadata
pca_raw <- fread(
  "Software/binMAF20NA10_PCA.eigenvec",
  header = FALSE
)

metadata <- fread(
  "data/Isolate_Host_origin_CORRECTED.txt"
)


# Format PCA dataframe
pca_plot_df <- data.frame(
  Sample = gsub("\\.", "_", pca_raw$V1),
  PC1 = pca_raw$V3,
  PC2 = pca_raw$V4,
  PC3 = pca_raw$V5
)

# Fix leading zeros
pca_plot_df$Sample <- sub("^1_", "01_", pca_plot_df$Sample)
pca_plot_df$Sample <- sub("^2_", "02_", pca_plot_df$Sample)

# Merge metadata
pca_plot_df2 <- merge(
  pca_plot_df,
  metadata,
  by.x = "Sample",
  by.y = "Isolate"
)

cat("Matching samples:", nrow(pca_plot_df2), "\n")

# Geography categories
pca_plot_df2$Geography <- ifelse(
  pca_plot_df2$Origin == "California",
  "California",
  "Other"
)

pca_plot_df2$Geography <- ifelse(
  pca_plot_df2$Origin %in%
    c("Netherlands", "Switzerland", "UnitedKingdom"),
  "Europe",
  pca_plot_df2$Geography
)

# Variance explained
eigenval <- fread(
  "Software/binMAF20NA10_PCA.eigenval",
  header = FALSE
)

pc_var <- eigenval$V1 / sum(eigenval$V1) * 100


# Scree plot dataframe
scree_df <- data.frame(
  PC_num = seq_along(pc_var),
  Variance = pc_var
)

# Shared PCA theme
pca_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )


# Panel A: PC1 vs PC2
p_pc12 <- ggplot(
  pca_plot_df2,
  aes(
    PC1,
    PC2,
    color = Host2,
    shape = Geography
  )
) +
  geom_point(size = 3.5, alpha = 0.85) +
  pca_theme +
  labs(
    x = paste0("PC1 (", round(pc_var[1], 2), "%)"),
    y = paste0("PC2 (", round(pc_var[2], 2), "%)"),
    color = "Host",
    shape = "Geography"
  )
-
# Panel B: PC1 vs PC3
p_pc13 <- ggplot(
  pca_plot_df2,
  aes(
    PC1,
    PC3,
    color = Host2,
    shape = Geography
  )
) +
  geom_point(size = 3.5, alpha = 0.85) +
  pca_theme +
  labs(
    x = paste0("PC1 (", round(pc_var[1], 2), "%)"),
    y = paste0("PC3 (", round(pc_var[3], 2), "%)"),
    color = "Host",
    shape = "Geography"
  )
-
# Panel C: PC2 vs PC3
p_pc23 <- ggplot(
  pca_plot_df2,
  aes(
    PC2,
    PC3,
    color = Host2,
    shape = Geography
  )
) +
  geom_point(size = 3.5, alpha = 0.85) +
  pca_theme +
  labs(
    x = paste0("PC2 (", round(pc_var[2], 2), "%)"),
    y = paste0("PC3 (", round(pc_var[3], 2), "%)"),
    color = "Host",
    shape = "Geography"
  )

# Panel D: Scree plot
p_scree <- ggplot(
  scree_df,
  aes(PC_num, Variance)
) +
  geom_point(size = 2.5) +
  geom_line() +
  theme_bw(base_size = 12) +
  labs(
    x = "Principal component",
    y = "Variance explained (%)"
  )

# turn figure into a Multi-panel figure
p_multi <- (
  p_pc12 | p_pc13
) / (
  p_pc23 | p_scree
) +
  plot_annotation(tag_levels = "A")

print(p_multi)

# Save figure
ggsave(
  "Plots/PCA/PCA_population_structure_multipanel.png",
  p_multi,
  width = 13,
  height = 10
)
