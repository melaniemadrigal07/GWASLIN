    library(ggtext)
library(wesanderson)
library(dplyr)
library(ggplot2)

#  Match violin/UpSet colors (Isle of Dogs) for trait class ----
pal <- wes_palette("IsleofDogs1", 2)
class_cols <- c(
  Structural = pal[1],
  Metabolic  = pal[2]
)

# ---- Direction-of-change colors (keep your current ones) ----
cb_palette <- c(
  toward_preferred = "#0072B2",
  toward_rare      = "#D55E00",
  no_change        = "#B0B0B0"
)

syn_with_usage_plot <- syn_with_usage %>%
  mutate(
    Gene_class = if_else(`Gene ID` == "Bcin07g05680", "Metabolic", "Structural"),
    Gene_class = factor(Gene_class, levels = c("Structural","Metabolic")),
    gene_lab = paste0(
      "<b><span style='color:", class_cols[as.character(Gene_class)], ";'>",
      `Gene ID`,
      "</span></b>"
    )
  )

p <- ggplot(syn_with_usage_plot, aes(
  y = reorder(gene_lab, DeltaUsage),
  x = DeltaUsage,
  fill = TowardPreferred
)) +
  geom_col(width = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, color = "black") +
  
  # invisible points only to generate a "Gene class" legend
  geom_point(aes(color = Gene_class), alpha = 0) +
  
  scale_fill_manual(values = cb_palette, name = "Direction of change") +
  scale_color_manual(values = class_cols, name = "Gene class") +
  
  guides(
    fill  = guide_legend(order = 1),
    color = guide_legend(order = 2, override.aes = list(alpha = 1, size = 4))
  ) +
  
  theme_classic(base_size = 13) +
  labs(
    x = expression(Delta~Codon~Usage~"(" * Alt - Ref * ", % )"),
    y = "Gene"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 11),
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

p
#save plot
ggsave("Fig_CUB.pdf", plot = p,
       width = 6.5, height = 4.5, units = "in", dpi = 600)
