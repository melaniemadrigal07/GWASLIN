# GWAS SNP annotation summary with Model comparisons
# Separates:
#   1) SNP-trait associations = rows
#   2) Unique SNPs = distinct SNP Id per model

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)

# ---------- 1) Helper functions ----------

clean_context_bracket <- function(x) {
  x %>%
    str_replace_all("\\\\>", ">") %>%
    str_replace_all("<span[^>]*>\\s*([ACGT])\\s*</span>", "[\\1]") %>%
    str_replace_all("<[^>]+>", "") %>%
    str_replace_all("&lt;", "<") %>%
    str_replace_all("&gt;", ">") %>%
    str_replace_all("&amp;", "&") %>%
    str_replace_all("\\s+", "")
}

rev_comp <- function(seq) {
  if (is.na(seq)) return(NA_character_)
  s <- chartr("ACGT", "TGCA", toupper(seq))
  paste(rev(strsplit(s, "")[[1]]), collapse = "")
}

normalize_strand <- function(x) {
  x2 <- tolower(trimws(as.character(x)))
  case_when(
    x2 %in% c("+", "forward", "fwd", "plus", "sense") ~ "+",
    x2 %in% c("-", "reverse", "rev", "minus", "antisense") ~ "-",
    TRUE ~ NA_character_
  )
}

get_codons_from_ctx_strand <- function(ctx, cds_pos, alt_base, strand = "+") {
  if (is.na(ctx) || ctx == "N/A" || is.na(cds_pos) || is.na(alt_base) || is.na(strand)) return(NULL)
  
  alt_base <- toupper(alt_base)
  if (!alt_base %in% c("A","C","G","T")) return(NULL)
  
  pos_in_codon_coding <- ((cds_pos - 1) %% 3) + 1
  pos_in_codon_genomic <- if (strand == "-") 4 - pos_in_codon_coding else pos_in_codon_coding
  
  codon_raw <- case_when(
    pos_in_codon_genomic == 1 ~ str_extract(ctx, "\\[[ACGT]\\][ACGT]{2}"),
    pos_in_codon_genomic == 2 ~ str_extract(ctx, "[ACGT]\\[[ACGT]\\][ACGT]"),
    pos_in_codon_genomic == 3 ~ str_extract(ctx, "[ACGT]{2}\\[[ACGT]\\]"),
    TRUE ~ NA_character_
  )
  
  if (is.na(codon_raw)) return(NULL)
  
  refCodon_gen <- str_replace_all(codon_raw, "\\[|\\]", "")
  if (nchar(refCodon_gen) != 3) return(NULL)
  
  v <- strsplit(refCodon_gen, "")[[1]]
  v[pos_in_codon_genomic] <- alt_base
  altCodon_gen <- paste(v, collapse = "")
  
  if (strand == "-") {
    refCodon <- rev_comp(refCodon_gen)
    altCodon <- rev_comp(altCodon_gen)
  } else {
    refCodon <- refCodon_gen
    altCodon <- altCodon_gen
  }
  
  tibble(
    RefCodon = refCodon,
    AltCodon = altCodon,
    pos_in_codon = pos_in_codon_coding
  )
}

codon_table <- c(
  TTT="F", TTC="F", TTA="L", TTG="L", TCT="S", TCC="S", TCA="S", TCG="S",
  TAT="Y", TAC="Y", TAA="*", TAG="*", TGT="C", TGC="C", TGA="*", TGG="W",
  CTT="L", CTC="L", CTA="L", CTG="L", CCT="P", CCC="P", CCA="P", CCG="P",
  CAT="H", CAC="H", CAA="Q", CAG="Q", CGT="R", CGC="R", CGA="R", CGG="R",
  ATT="I", ATC="I", ATA="I", ATG="M", ACT="T", ACC="T", ACA="T", ACG="T",
  AAT="N", AAC="N", AAA="K", AAG="K", AGT="S", AGC="S", AGA="R", AGG="R",
  GTT="V", GTC="V", GTA="V", GTG="V", GCT="A", GCC="A", GCA="A", GCG="A",
  GAT="D", GAC="D", GAA="E", GAG="E", GGT="G", GGC="G", GGA="G", GGG="G"
)

codon_to_aa <- function(codon) {
  codon <- toupper(codon)
  aa <- unname(codon_table[codon])
  ifelse(is.na(aa), NA_character_, aa)
}

# ---------- 2) Load files ----------

snps_raw <- read_csv("NgsSnpBySourceId_Summary.csv", show_col_types = FALSE)
gwas <- read_csv("SignificantGenomewideSNPS.csv", show_col_types = FALSE)

# ---------- 3) Attach GWAS model/trait info ----------

gwas_ids <- gwas %>%
  mutate(
    SNP_ID = paste0("NGS_SNP.bcin_chr_", CHROM, ".", POS)
  ) %>%
  select(SNP_ID, CHROM, POS, ID, Trait, Environment, Model, P, BETA, AF) %>%
  distinct()

snps_joined <- snps_raw %>%
  rename(
    MajorAllele = `Major Allele`,
    MinorAllele = `Minor Allele`
  ) %>%
  left_join(
    gwas_ids,
    by = c("SNP Id" = "SNP_ID")
  ) %>%
  mutate(
    Coding = tolower(Coding)
  )

cat("\n=== Model matching check ===\n")
print(table(snps_joined$Model, useNA = "ifany"))

# ---------- 4) Association-level summaries ----------

assoc_by_model <- snps_joined %>%
  filter(!is.na(Model)) %>%
  count(Model, name = "n_snp_trait_associations")

cat("\n=== SNP-trait associations by model ===\n")
print(assoc_by_model)

assoc_by_chr_model <- snps_joined %>%
  filter(!is.na(Model)) %>%
  count(Model, CHROM, name = "n_associations") %>%
  arrange(Model, desc(n_associations))

cat("\n=== SNP-trait associations by chromosome and model ===\n")
print(assoc_by_chr_model, n = Inf)

coding_assoc_by_model <- snps_joined %>%
  filter(!is.na(Model)) %>%
  count(Model, Coding, name = "n_associations") %>%
  group_by(Model) %>%
  mutate(prop = n_associations / sum(n_associations)) %>%
  ungroup()

cat("\n=== Coding vs non-coding associations by model ===\n")
print(coding_assoc_by_model)

# ---------- 5) Unique SNP tables ----------

unique_snps_by_model_full <- snps_joined %>%
  filter(!is.na(Model)) %>%
  distinct(Model, `SNP Id`, .keep_all = TRUE)

unique_snps_by_model <- unique_snps_by_model_full %>%
  count(Model, name = "n_unique_snps")

cat("\n=== Unique SNPs by model ===\n")
print(unique_snps_by_model)

unique_chr_by_model <- unique_snps_by_model_full %>%
  count(Model, CHROM, name = "n_unique_snps") %>%
  arrange(Model, desc(n_unique_snps))

cat("\n=== Unique SNPs by chromosome and model ===\n")
print(unique_chr_by_model, n = Inf)

unique_coding_by_model <- unique_snps_by_model_full %>%
  count(Model, Coding, name = "n_unique_snps") %>%
  group_by(Model) %>%
  mutate(prop = n_unique_snps / sum(n_unique_snps)) %>%
  ungroup()

cat("\n=== Unique coding vs non-coding SNPs by model ===\n")
print(unique_coding_by_model)

# ---------- 6) Build codons and SNP classes ----------

snps_codons <- snps_joined %>%
  mutate(
    cds_pos = parse_integer(as.character(na_if(`Position in CDS`, "N/A"))),
    strand  = normalize_strand(`Gene strand`),
    ctx     = clean_context_bracket(`SNP context`),
    alt     = str_extract(toupper(as.character(MinorAllele)), "[ACGT]")
  ) %>%
  filter(Coding == "coding", !is.na(cds_pos), !is.na(strand), !is.na(alt)) %>%
  rowwise() %>%
  mutate(tmp = list(get_codons_from_ctx_strand(ctx, cds_pos, alt, strand))) %>%
  ungroup() %>%
  unnest_wider(tmp) %>%
  mutate(
    RefAA = codon_to_aa(RefCodon),
    AltAA = codon_to_aa(AltCodon),
    SNP_Class = case_when(
      is.na(RefAA) | is.na(AltAA) ~ "unknown",
      RefAA == AltAA              ~ "synonymous",
      RefAA != "*" & AltAA != "*" ~ "missense",
      RefAA != "*" & AltAA == "*" ~ "nonsense",
      RefAA == "*" & AltAA != "*" ~ "stop_loss",
      TRUE                        ~ "other"
    )
  )

cat("\n=== Coding SNP associations with codons extracted ===\n")
cat("n =", nrow(snps_codons), "\n")

# ---------- 7) Coding classes: association-level ----------

coding_class_assoc_by_model <- snps_codons %>%
  filter(!is.na(Model)) %>%
  count(Model, SNP_Class, name = "n_associations") %>%
  group_by(Model) %>%
  mutate(prop = n_associations / sum(n_associations)) %>%
  ungroup()

cat("\n=== Coding SNP classes by model, association-level ===\n")
print(coding_class_assoc_by_model)

# ---------- 8) Coding classes: unique SNP-level ----------

unique_codons_by_model <- snps_codons %>%
  filter(!is.na(Model)) %>%
  distinct(Model, `SNP Id`, .keep_all = TRUE)

coding_class_unique_by_model <- unique_codons_by_model %>%
  count(Model, SNP_Class, name = "n_unique_snps") %>%
  group_by(Model) %>%
  mutate(prop = n_unique_snps / sum(n_unique_snps)) %>%
  ungroup()

cat("\n=== Coding SNP classes by model, unique SNP-level ===\n")
print(coding_class_unique_by_model)

# ---------- 9) Three-category unique SNP summary ----------

noncoding_unique <- unique_snps_by_model_full %>%
  filter(Coding != "coding") %>%
  transmute(
    Model,
    `SNP Id`,
    CHROM,
    POS,
    Category3 = "noncoding"
  )

coding_unique <- unique_codons_by_model %>%
  transmute(
    Model,
    `SNP Id`,
    CHROM,
    POS,
    Category3 = case_when(
      SNP_Class == "synonymous" ~ "coding-synonymous",
      SNP_Class %in% c("missense", "nonsense", "stop_loss", "other", "unknown") ~ "coding-missense_or_other",
      TRUE ~ "coding-unknown"
    )
  )

three_cat_unique_by_model <- bind_rows(noncoding_unique, coding_unique) %>%
  count(Model, Category3, name = "n_unique_snps") %>%
  group_by(Model) %>%
  mutate(prop = n_unique_snps / sum(n_unique_snps)) %>%
  ungroup()

cat("\n=== Three-category unique SNP summary by model ===\n")
print(three_cat_unique_by_model)

# K-model manuscript checks
k_unique_chr <- unique_snps_by_model_full %>%
  filter(Model == "K") %>%
  count(CHROM, name = "n_unique_snps") %>%
  arrange(desc(n_unique_snps))

k_assoc_chr <- snps_joined %>%
  filter(Model == "K") %>%
  count(CHROM, name = "n_associations") %>%
  arrange(desc(n_associations))

k_three_cat <- three_cat_unique_by_model %>%
  filter(Model == "K")

cat("\n=== K model: unique SNPs by chromosome ===\n")
print(k_unique_chr, n = Inf)

cat("\n=== K model: SNP-trait associations by chromosome ===\n")
print(k_assoc_chr, n = Inf)

cat("\n=== K model: three-category unique SNP annotation ===\n")
print(k_three_cat)

cat("\n=== K model manuscript totals ===\n")
cat("K associations =", sum(k_assoc_chr$n_associations), "\n")
cat("K unique SNPs =", sum(k_unique_chr$n_unique_snps), "\n")

# ---------- 10) Plots ----------

cb <- c(
  non_coding  = "#999999",
  `non-coding` = "#999999",
  coding      = "#000000",
  synonymous  = "#0072B2",
  missense    = "#D55E00",
  nonsense    = "#E69F00",
  stop_loss   = "#009E73",
  other       = "#56B4E9",
  unknown     = "#CC79A7",
  noncoding   = "#999999",
  `coding-synonymous` = "#0072B2",
  coding_missense_or_other = "#D55E00",
  `coding-missense_or_other` = "#D55E00"
)

p_assoc_chr_k <- k_assoc_chr %>%
  ggplot(aes(x = reorder(as.factor(CHROM), -n_associations), y = n_associations)) +
  geom_col(width = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    x = "Chromosome",
    y = "SNP-trait associations",
    title = "K model SNP-trait associations by chromosome"
  )

print(p_assoc_chr_k)

p_unique_chr_k <- k_unique_chr %>%
  ggplot(aes(x = reorder(as.factor(CHROM), -n_unique_snps), y = n_unique_snps)) +
  geom_col(width = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    x = "Chromosome",
    y = "Unique SNPs",
    title = "K model unique SNPs by chromosome"
  )

print(p_unique_chr_k)

p_unique_coding_model <- unique_coding_by_model %>%
  ggplot(aes(x = Coding, y = n_unique_snps, fill = Coding)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Model) +
  theme_classic(base_size = 13) +
  labs(
    x = "Annotation class",
    y = "Unique SNP count",
    title = "Unique coding vs non-coding SNPs by model"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(p_unique_coding_model)

p_unique_class_model <- coding_class_unique_by_model %>%
  mutate(
    SNP_Class = factor(
      SNP_Class,
      levels = c("synonymous", "missense", "nonsense", "stop_loss", "other", "unknown")
    )
  ) %>%
  ggplot(aes(x = SNP_Class, y = n_unique_snps, fill = SNP_Class)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Model) +
  scale_fill_manual(values = cb, drop = FALSE) +
  theme_classic(base_size = 13) +
  labs(
    x = "Coding SNP class",
    y = "Unique SNP count",
    title = "Unique coding SNP classes by model"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(p_unique_class_model)

p_three_cat <- three_cat_unique_by_model %>%
  ggplot(aes(x = Category3, y = n_unique_snps, fill = Category3)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Model) +
  scale_fill_manual(values = cb, drop = FALSE) +
  theme_classic(base_size = 13) +
  labs(
    x = "Unique SNP annotation class",
    y = "Unique SNP count",
    title = "Three-category unique SNP summary by model"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(p_three_cat)

# -----------------------------
# Nested annotation plot by model: UNIQUE SNPs only
# -----------------------------

snps2_unique <- snps2 %>%
  filter(!is.na(Model)) %>%
  distinct(Model, `SNP Id`, .keep_all = TRUE) %>%
  mutate(
    SNP_Class = case_when(
      Coding != "coding" ~ "non_coding",
      FungiDB_Effect == "synonymous" ~ "synonymous",
      FungiDB_Effect == "nonsynonymous" ~ "missense",
      TRUE ~ "unknown"
    ),
    Broad_Class = case_when(
      SNP_Class == "non_coding" ~ "non_coding",
      TRUE ~ "coding"
    )
  )

# Broad coding vs non-coding
broad_plot_unique <- snps2_unique %>%
  count(Model, Broad_Class, name = "n") %>%
  mutate(
    SNP_Class = Broad_Class,
    group = "All unique SNPs"
  ) %>%
  select(Model, SNP_Class, group, n)

# Coding subclasses only
coding_plot_unique <- snps2_unique %>%
  filter(Broad_Class == "coding") %>%
  count(Model, SNP_Class, name = "n") %>%
  mutate(
    group = "Coding SNPs only"
  ) %>%
  select(Model, SNP_Class, group, n)

plot_df_unique <- bind_rows(broad_plot_unique, coding_plot_unique) %>%
  group_by(Model, group) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

cb2 <- c(
  non_coding = "#999999",
  coding     = "#1B3A4B",
  synonymous = "#2C7FB8",
  missense   = "#41B6C4",
  unknown    = "#CCCCCC"
)

p_nested_model_unique <- ggplot(plot_df_unique, aes(x = group, y = n, fill = SNP_Class)) +
  geom_col(width = 0.65) +
  facet_wrap(~ Model) +
  scale_fill_manual(values = cb2, drop = FALSE) +
  theme_classic(base_size = 13) +
  labs(
    #title = "Unique Genome-Wide Significant SNP Annotation by GWAS Model",
    x = NULL,
    y = "Unique SNP count",
    fill = "Annotation class"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 11),
    strip.text = element_text(face = "bold")
  )

p_nested_model_unique

ggsave(
  "Fig_UNIQUE_SNP_annotation_nested_by_model.pdf",
  plot = p_nested_model_unique,
  width = 10,
  height = 4.5,
  units = "in",
  dpi = 600
)

# ---------- 11) Missense details ----------

aa_props <- tibble::tribble(
  ~AA, ~AA_Group,
  "A", "nonpolar", "V", "nonpolar", "L", "nonpolar", "I", "nonpolar",
  "M", "nonpolar", "P", "nonpolar", "G", "special",
  "F", "aromatic", "W", "aromatic", "Y", "aromatic",
  "S", "polar", "T", "polar", "N", "polar", "Q", "polar", "C", "polar",
  "K", "basic", "R", "basic", "H", "basic",
  "D", "acidic", "E", "acidic", "*", "stop"
)

missense_snps <- snps_codons %>%
  filter(SNP_Class == "missense") %>%
  mutate(
    RefBase = str_match(ctx, "\\[([ACGT])\\]")[,2],
    AltBase = alt,
    NT_Change = paste0(RefBase, ">", AltBase),
    NT_Change_Type = case_when(
      RefBase %in% c("A","G") & AltBase %in% c("A","G") & RefBase != AltBase ~ "transition",
      RefBase %in% c("C","T") & AltBase %in% c("C","T") & RefBase != AltBase ~ "transition",
      RefBase %in% c("A","C","G","T") & AltBase %in% c("A","C","G","T") & RefBase != AltBase ~ "transversion",
      TRUE ~ NA_character_
    ),
    AA_Change = paste0(RefAA, ">", AltAA),
    Codon_Change = paste0(RefCodon, ">", AltCodon)
  ) %>%
  left_join(aa_props %>% rename(RefAA = AA, Ref_AA_Group = AA_Group), by = "RefAA") %>%
  left_join(aa_props %>% rename(AltAA = AA, Alt_AA_Group = AA_Group), by = "AltAA") %>%
  mutate(
    Missense_Type = case_when(
      is.na(Ref_AA_Group) | is.na(Alt_AA_Group) ~ "unknown",
      Ref_AA_Group == Alt_AA_Group ~ "conservative_like",
      TRUE ~ "nonconservative_like"
    )
  )

missense_unique_by_model <- missense_snps %>%
  filter(!is.na(Model)) %>%
  distinct(Model, `SNP Id`, .keep_all = TRUE) %>%
  count(Model, Missense_Type, name = "n_unique_missense_snps") %>%
  group_by(Model) %>%
  mutate(prop = n_unique_missense_snps / sum(n_unique_missense_snps)) %>%
  ungroup()

cat("\n=== Unique missense SNP type by model ===\n")
print(missense_unique_by_model)

# ---------- 12) Save outputs ----------

write_csv(assoc_by_model, "gwas_ASSOCIATIONS_by_model.csv")
write_csv(assoc_by_chr_model, "gwas_ASSOCIATIONS_by_chromosome_and_model.csv")
write_csv(coding_assoc_by_model, "gwas_ASSOCIATIONS_coding_vs_noncoding_by_model.csv")
write_csv(coding_class_assoc_by_model, "gwas_ASSOCIATIONS_coding_class_by_model.csv")

write_csv(unique_snps_by_model, "gwas_UNIQUE_snps_by_model.csv")
write_csv(unique_chr_by_model, "gwas_UNIQUE_snps_by_chromosome_and_model.csv")
write_csv(unique_coding_by_model, "gwas_UNIQUE_coding_vs_noncoding_by_model.csv")
write_csv(coding_class_unique_by_model, "gwas_UNIQUE_coding_snp_class_by_model.csv")
write_csv(three_cat_unique_by_model, "gwas_UNIQUE_three_category_summary_by_model.csv")

write_csv(k_unique_chr, "K_model_UNIQUE_snps_by_chromosome.csv")
write_csv(k_assoc_chr, "K_model_ASSOCIATIONS_by_chromosome.csv")
write_csv(k_three_cat, "K_model_UNIQUE_three_category_summary.csv")

write_csv(snps_joined, "gwas_snps_joined_with_model_trait_info.csv")
write_csv(snps_codons, "gwas_snps_with_codons_and_snp_class.csv")
write_csv(missense_snps, "gwas_missense_snp_details.csv")
write_csv(missense_unique_by_model, "gwas_UNIQUE_missense_by_model.csv")

ggsave("plot_K_model_associations_by_chromosome.png", p_assoc_chr_k, width = 7, height = 5, dpi = 300)
ggsave("plot_K_model_unique_snps_by_chromosome.png", p_unique_chr_k, width = 7, height = 5, dpi = 300)
ggsave("plot_unique_coding_vs_noncoding_by_model.png", p_unique_coding_model, width = 8, height = 5, dpi = 300)
ggsave("plot_unique_coding_snp_classes_by_model.png", p_unique_class_model, width = 8, height = 5, dpi = 300)
ggsave("plot_unique_three_category_by_model.png", p_three_cat, width = 8, height = 5, dpi = 300)

cat("\nSaved updated CSVs and plots successfully.\n")
