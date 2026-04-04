library(dplyr)
library(stringr)
library(ggplot2)
library(readr)

# read your FungiDB annotation table
ann <- read_csv("lfungidb.csv", show_col_types = FALSE)

# assign broad functional categories from PFam descriptions
ann2 <- ann %>%
  mutate(
    PFam_Desc = coalesce(`PFam Description`, ""),
    Category = case_when(
      str_detect(PFam_Desc, regex("major facilitator|transporter|ABC", ignore_case = TRUE)) ~ "Transporter",
      str_detect(PFam_Desc, regex("cytochrome P450|oxidoreductase|dehydrogenase|NAD\\(P\\)-binding", ignore_case = TRUE)) ~ "Oxidoreductase",
      str_detect(PFam_Desc, regex("histidine kinase|signal transduction|receiver domain|kinase", ignore_case = TRUE)) ~ "Signal transduction",
      str_detect(PFam_Desc, regex("clathrin|ARID|vesicle|trafficking", ignore_case = TRUE)) ~ "Vesicle trafficking",
      str_detect(PFam_Desc, regex("glycosyl transferase|trehalose|hydrolase|neuraminidase|epimerase|metabolism", ignore_case = TRUE)) ~ "Metabolism",
      str_detect(PFam_Desc, regex("heterokaryon incompatibility|defense|compatibility", ignore_case = TRUE)) ~ "Cell compatibility / defense",
      PFam_Desc == "" | PFam_Desc == "N/A" ~ "Unknown",
      TRUE ~ "Other"
    )
  )

# make a clean annotation table that reflects the new categorization
annotation_table <- ann2 %>%
  distinct(`Gene ID`, `Interpro ID`, `Interpro Description`, PFam_Desc, Category) %>%
  rename(
    `InterPro ID` = `Interpro ID`,
    `InterPro Description` = `Interpro Description`,
    `PFam Description` = PFam_Desc,
    `Assigned Category` = Category
  )

# view the updated table
print(annotation_table)

# save the updated table
write_csv(annotation_table, "functional_annotation_table.csv")

# keep one row per gene for counting categories
counts <- ann2 %>%
  distinct(`Gene ID`, Category) %>%
  count(Category, sort = TRUE)

# set plotting order
counts <- counts %>%
  mutate(Category = factor(
    Category,
    levels = rev(c(
      "Metabolism",
      "Oxidoreductase",
      "Transporter",
      "Signal transduction",
      "Vesicle trafficking",
      "Cell compatibility / defense",
      "Unknown",
      "Other"
    ))
  )) %>%
  arrange(Category)

# make the bar plot
p_fun <- ggplot(counts, aes(x = Category, y = n)) +
  geom_col(fill = "#5C8D89", width = 0.8) +
  coord_flip() +
  labs(
    x = "Gene functional category",
    y = "Number of genes with significant SNPs"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 13)
  )

print(p_fun)

# save the plot
ggsave(
  "functional_gene_categories.pdf",
  plot = p_fun,
  width = 9,
  height = 5.5,
  units = "in",
  device = cairo_pdf
)
