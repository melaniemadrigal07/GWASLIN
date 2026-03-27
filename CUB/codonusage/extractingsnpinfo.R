library(dplyr)
library(readr)

gwas <- read_csv("SignificantGenomewideSNPS.csv", show_col_types = FALSE)

ids <- gwas %>%
  filter(!is.na(CHROM), !is.na(POS)) %>%
  transmute(SNP_ID = paste0("NGS_SNP.bcin_chr_", CHROM, ".", POS)) %>%
  distinct()

# WRITE AS PLAIN TEXT (one ID per line, no header, no quotes)
writeLines(ids$SNP_ID, "fungidb_snp_ids.txt")
