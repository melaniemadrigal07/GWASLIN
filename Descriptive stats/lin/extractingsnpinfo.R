library(dplyr)
library(readr)

gwas <- read_csv("SignificantGenomewideSNPS_withEffects.csv", show_col_types = FALSE)

ids <- gwas %>%
  filter(!is.na(CHROM), !is.na(POS)) %>%
  transmute(SNP_ID = paste0("NGS_SNP.bcin_chr_", CHROM, ".", POS)) %>%
  distinct()

# WRITE AS PLAIN TEXT (one ID per line, no header, no quotes)
writeLines(ids$SNP_ID, "fungidb_snp_ids.txt")

library(readr)
library(dplyr)

ann <- read_csv("NgsSnpBySourceId_Summary.csv", show_col_types = FALSE)

gene_ids <- ann %>%
  filter(!is.na(`Gene ID`)) %>%
  distinct(`Gene ID`) %>%
  arrange(`Gene ID`)

write_csv(gene_ids, "unique_gene_ids_from_significant_snps.csv")
