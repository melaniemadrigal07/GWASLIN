# GWAS SNP annotation summary
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)

#1) Clean FungiDB SNP context so the SNP base is marked as [A/C/G/T] 
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

#2) Strand helpers 
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

# 3) Strand-aware codon extraction from ctx + cds_pos + alt base
# ctx is genomic-forward and contains exactly one [X] for the SNP base
get_codons_from_ctx_strand <- function(ctx, cds_pos, alt_base, strand = "+") {
  if (is.na(ctx) || ctx == "N/A" || is.na(cds_pos) || is.na(alt_base) || is.na(strand)) return(NULL)
  
  alt_base <- toupper(alt_base)
  if (!alt_base %in% c("A","C","G","T")) return(NULL)
  
  # position within codon along the coding sequence
  pos_in_codon_coding <- ((cds_pos - 1) %% 3) + 1
  
  # if reverse strand, the SNP position within the genomic-forward triplet is mirrored
  pos_in_codon_genomic <- if (strand == "-") 4 - pos_in_codon_coding else pos_in_codon_coding
  
  codon_raw <- dplyr::case_when(
    pos_in_codon_genomic == 1 ~ stringr::str_extract(ctx, "\\[[ACGT]\\][ACGT]{2}"),
    pos_in_codon_genomic == 2 ~ stringr::str_extract(ctx, "[ACGT]\\[[ACGT]\\][ACGT]"),
    pos_in_codon_genomic == 3 ~ stringr::str_extract(ctx, "[ACGT]{2}\\[[ACGT]\\]"),
    TRUE ~ NA_character_
  )
  if (is.na(codon_raw)) return(NULL)
  
  refCodon_gen <- stringr::str_replace_all(codon_raw, "\\[|\\]", "")
  if (nchar(refCodon_gen) != 3) return(NULL)
  
  # substitute alt base in the genomic-forward triplet
  v <- strsplit(refCodon_gen, "")[[1]]
  v[pos_in_codon_genomic] <- alt_base
  altCodon_gen <- paste(v, collapse = "")
  
  # convert to coding-strand codons if reverse strand
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

# ---------- 4) Codon -> amino acid (standard code) ----------
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

# ---------- 5) Load your GWAS SNP table ----------
snps <- read_csv("NgsSnpBySourceId_Summary.csv", show_col_types = FALSE)

# ---------- 6) Coding vs non-coding (simple) ----------
coding_breakdown <- snps %>%
  mutate(Coding = tolower(Coding)) %>%
  count(Coding, name = "n") %>%
  mutate(prop = n / sum(n))

cat("\n=== Coding vs non-coding breakdown ===\n")
print(coding_breakdown)

# ---------- 7) Build SNP-level codon + AA + class (coding SNPs with cds_pos + strand) ----------
# NOTE: Adjust these column names if yours differ:
#   "Position in CDS", "Gene strand", "SNP context", "Minor Allele"
snps_codons <- snps %>%
  mutate(
    Coding = tolower(Coding),
    cds_pos = parse_integer(as.character(na_if(`Position in CDS`, "N/A"))),
    strand  = normalize_strand(`Gene strand`),
    ctx     = clean_context_bracket(`SNP context`),
    alt     = str_extract(toupper(as.character(`Minor Allele`)), "[ACGT]")
  ) %>%
  filter(Coding == "coding", !is.na(cds_pos), !is.na(strand), !is.na(alt)) %>%
  rowwise() %>%
  mutate(tmp = list(get_codons_from_ctx_strand(ctx, cds_pos, alt, strand))) %>%
  ungroup() %>%
  unnest_wider(tmp) %>%
  mutate(
    RefAA = codon_to_aa(RefCodon),
    AltAA = codon_to_aa(AltCodon),
    
    # SNP-level coding class
    SNP_Class = case_when(
      is.na(RefAA) | is.na(AltAA) ~ "unknown",
      RefAA == AltAA              ~ "synonymous",
      RefAA != "*" & AltAA != "*" ~ "missense",
      RefAA != "*" & AltAA == "*" ~ "nonsense",
      RefAA == "*" & AltAA != "*" ~ "stop_loss",
      TRUE                        ~ "other"
    )
  )

cat("\n=== Coding SNPs with codons extracted ===\n")
cat("n =", nrow(snps_codons), "\n")

# ---------- 8) SNP-level class counts (coding only) ----------
coding_class_summary <- snps_codons %>%
  count(SNP_Class, name = "n") %>%
  mutate(prop = n / sum(n)) %>%
  arrange(desc(n))

cat("\n=== SNP-level classes (coding only) ===\n")
print(coding_class_summary)

# ---------- 9) Non-synonymous subclass breakdown (missense/nonsense/stop_loss only) ----------
nonsyn_summary <- snps_codons %>%
  filter(SNP_Class %in% c("missense", "nonsense", "stop_loss")) %>%
  count(SNP_Class, name = "n") %>%
  mutate(prop = n / sum(n)) %>%
  arrange(desc(n))

cat("\n=== Non-synonymous subclasses (coding) ===\n")
print(nonsyn_summary)

# ---------- 10) One combined table including non-coding ----------
full_summary <- bind_rows(
  snps %>%
    mutate(
      Coding = tolower(Coding),
      Broad_Class = if_else(Coding == "coding", NA_character_, "non_coding")
    ) %>%
    filter(Broad_Class == "non_coding") %>%
    count(Broad_Class, name = "n"),
  snps_codons %>%
    mutate(Broad_Class = paste0("coding_", SNP_Class)) %>%
    count(Broad_Class, name = "n")
) %>%
  group_by(Broad_Class) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(prop = n / sum(n)) %>%
  arrange(desc(n))

cat("\n=== Full summary (non-coding + coding subclasses) ===\n")
print(full_summary)

# plotting non-synonymous subclasses 
# Comment out if you don't want figures.
cb_subclass <- c(
  missense = "#0072B2",   # blue
  nonsense = "#D55E00",   # vermillion
  stop_loss = "#E69F00"   # orange
)

if (nrow(nonsyn_summary) > 0) {
  p1 <- ggplot(nonsyn_summary, aes(x = SNP_Class, y = n, fill = SNP_Class)) +
    geom_col(width = 0.65) +
    scale_fill_manual(values = cb_subclass) +
    theme_classic(base_size = 13) +
    labs(
      x = "Non-synonymous SNP subclass",
      y = "Count",
      title = "Non-synonymous GWAS SNP subclasses (coding)"
    ) +
    theme(legend.position = "none")
  
  print(p1)
}

# ---------- 12) Save outputs (optional) ----------
write_csv(coding_breakdown, "gwas_snp_coding_vs_noncoding_summary.csv")
write_csv(coding_class_summary, "gwas_snp_coding_snp_class_summary.csv")
write_csv(nonsyn_summary, "gwas_snp_nonsyn_subclass_summary.csv")
write_csv(full_summary, "gwas_snp_full_class_summary.csv")
write_csv(snps_codons, "gwas_snps_with_codons_and_snp_class.csv")

cat("\nSaved CSVs:\n")
cat("- gwas_snp_coding_vs_noncoding_summary.csv\n")
cat("- gwas_snp_coding_snp_class_summary.csv\n")
cat("- gwas_snp_nonsyn_subclass_summary.csv\n")


cat("- gwas_snp_full_class_summary.csv\n")
cat("- gwas_snps_with_codons_and_snp_class.csv\n")



# ============================
# ONE unified GWAS SNP summary table
# ============================

# ---- 1) Broad coding vs non-coding ----
broad_summary <- snps %>%
  mutate(
    Coding = tolower(Coding),
    SNP_Class = if_else(Coding == "coding", "coding", "non_coding")
  ) %>%
  count(SNP_Class, name = "n") %>%
  mutate(
    Category_Level = "broad",
    prop = n / sum(n)
  )

# ---- 2) Coding SNP classes (synonymous / missense / etc.) ----
coding_summary <- snps_codons %>%
  count(SNP_Class, name = "n") %>%
  mutate(
    Category_Level = "coding",
    prop = n / sum(n)
  )

# ---- 3) Combine into ONE table ----
gwas_snp_summary <- bind_rows(
  broad_summary %>% select(Category_Level, SNP_Class, n, prop),
  coding_summary %>% select(Category_Level, SNP_Class, n, prop)
) %>%
  arrange(Category_Level, desc(n))

# ---- 4) Inspect ----
print(gwas_snp_summary, n = Inf)

# ---- 5) Save ONE CSV ----
write_csv(gwas_snp_summary, "gwas_snp_summary_all_classes.csv")

library(ggplot2)
library(dplyr)
library(scales)

# Okabe–Ito color-blind safe palette
cb <- c(
  non_coding  = "#999999",
  synonymous  = "#0072B2",
  missense    = "#D55E00",
  nonsense    = "#E69F00",
  stop_loss   = "#009E73",
  other       = "#56B4E9",
  unknown     = "#CC79A7"
)

p_coding <- gwas_snp_summary %>%
  filter(Category_Level == "coding") %>%
  mutate(
    SNP_Class = factor(SNP_Class, levels = c("synonymous","missense","nonsense","stop_loss","other","unknown"))
  ) %>%
  ggplot(aes(x = SNP_Class, y = n, fill = SNP_Class)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = cb, drop = FALSE) +
  theme_classic(base_size = 13) +
  labs(
    title = "GWAS coding SNP classes",
    x = "Coding SNP class",
    y = "Count"
  ) +
  theme(legend.position = "none")

p_coding

p_all <- gwas_snp_summary %>%
  mutate(
    SNP_Class = case_when(
      Category_Level == "broad" & SNP_Class == "coding" ~ "coding_total",
      TRUE ~ SNP_Class
    ),
    SNP_Class = factor(
      SNP_Class,
      levels = c("non_coding","coding_total","synonymous","missense","nonsense","stop_loss","other","unknown")
    )
  ) %>%
  ggplot(aes(x = SNP_Class, y = n, fill = SNP_Class)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Category_Level, scales = "free_x") +
  scale_fill_manual(values = c(cb, coding_total = "#000000"), drop = FALSE) +
  theme_classic(base_size = 13) +
  labs(
    title = "Functional classes of GWAS SNPs",
    x = NULL,
    y = "Count"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

p_all

#if missense SNP present

# 1) Amino acid property lookup (for "type" of missense change)
aa_props <- tibble::tribble(
  ~AA, ~AA_Group,
  "A", "nonpolar",
  "V", "nonpolar",
  "L", "nonpolar",
  "I", "nonpolar",
  "M", "nonpolar",
  "P", "nonpolar",
  "G", "special",
  "F", "aromatic",
  "W", "aromatic",
  "Y", "aromatic",
  "S", "polar",
  "T", "polar",
  "N", "polar",
  "Q", "polar",
  "C", "polar",
  "K", "basic",
  "R", "basic",
  "H", "basic",
  "D", "acidic",
  "E", "acidic",
  "*", "stop"
)

# 2) Extract reference base from bracketed context (e.g., A[C]T -> C)
#    Your 'ctx' already contains the cleaned bracket format in snps_codons
missense_snps <- snps_codons %>%
  filter(SNP_Class == "missense") %>%
  mutate(
    RefBase = stringr::str_match(ctx, "\\[([ACGT])\\]")[,2],
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
      Ref_AA_Group == Alt_AA_Group              ~ "conservative_like",
      TRUE                                      ~ "nonconservative_like"
    )
  )

# 3) View detailed missense table
missense_snps %>%
  select(
    `Gene ID`, Chromosome, Location, `Gene strand`,
    RefBase, AltBase, NT_Change, NT_Change_Type,
    pos_in_codon,
    RefCodon, AltCodon, Codon_Change,
    RefAA, AltAA, AA_Change,
    Ref_AA_Group, Alt_AA_Group, Missense_Type
  ) %>%
  arrange(`Gene ID`, Chromosome, Location) %>%
  print(n = Inf, width = Inf)

# 4) Save detailed missense table
write_csv(
  missense_snps %>%
    select(
      `Gene ID`, Chromosome, Location, `Gene strand`,
      RefBase, AltBase, NT_Change, NT_Change_Type,
      pos_in_codon,
      RefCodon, AltCodon, Codon_Change,
      RefAA, AltAA, AA_Change,
      Ref_AA_Group, Alt_AA_Group, Missense_Type
    ),
  "gwas_missense_snp_details.csv"
)

# ============================
# SUMMARIES: which genes + what type
# ============================

# A) Which genes have missense SNPs (and how many)
missense_genes_summary <- missense_snps %>%
  count(`Gene ID`, sort = TRUE, name = "n_missense_snps")

print(missense_genes_summary, n = Inf)

write_csv(missense_genes_summary, "gwas_missense_genes_summary.csv")

# B) Missense type summary (conservative-like vs nonconservative-like)
missense_type_summary <- missense_snps %>%
  count(Missense_Type, sort = TRUE)

print(missense_type_summary)

write_csv(missense_type_summary, "gwas_missense_type_summary.csv")

# C) Nucleotide change type summary (transition/transversion)
missense_nt_type_summary <- missense_snps %>%
  count(NT_Change_Type, sort = TRUE)

print(missense_nt_type_summary)

write_csv(missense_nt_type_summary, "gwas_missense_nt_change_type_summary.csv")

# D) Amino acid substitutions observed
aa_change_summary <- missense_snps %>%
  count(AA_Change, sort = TRUE)

print(aa_change_summary, n = Inf)

write_csv(aa_change_summary, "gwas_missense_aa_change_summary.csv")


# AA property annotations: charge + size


aa_properties <- tibble::tribble(
  ~AA, ~Charge,   ~SizeClass, ~PolarityClass,
  "A", "neutral", "small",    "nonpolar",
  "V", "neutral", "medium",   "nonpolar",
  "L", "neutral", "large",    "nonpolar",
  "I", "neutral", "large",    "nonpolar",
  "M", "neutral", "large",    "nonpolar",
  "P", "neutral", "medium",   "special",
  "G", "neutral", "small",    "special",
  "F", "neutral", "large",    "aromatic",
  "W", "neutral", "large",    "aromatic",
  "Y", "neutral", "large",    "aromatic",
  "S", "neutral", "small",    "polar",
  "T", "neutral", "medium",   "polar",
  "N", "neutral", "medium",   "polar",
  "Q", "neutral", "large",    "polar",
  "C", "neutral", "small",    "polar",
  "K", "positive","large",    "basic",
  "R", "positive","large",    "basic",
  "H", "positive","medium",   "basic",   # simplified; histidine can be context-dependent
  "D", "negative","medium",   "acidic",
  "E", "negative","large",    "acidic",
  "*", "stop",    "na",       "stop"
)

missense_props <- missense_snps %>%
  left_join(
    aa_properties %>%
      rename(
        RefAA = AA,
        Ref_Charge = Charge,
        Ref_Size = SizeClass,
        Ref_Polarity = PolarityClass
      ),
    by = "RefAA"
  ) %>%
  left_join(
    aa_properties %>%
      rename(
        AltAA = AA,
        Alt_Charge = Charge,
        Alt_Size = SizeClass,
        Alt_Polarity = PolarityClass
      ),
    by = "AltAA"
  ) %>%
  mutate(
    Charge_Change = paste0(Ref_Charge, " -> ", Alt_Charge),
    Charge_Shift = case_when(
      is.na(Ref_Charge) | is.na(Alt_Charge) ~ "unknown",
      Ref_Charge == Alt_Charge ~ "no_charge_change",
      TRUE ~ "charge_change"
    ),
    
    Size_Change = paste0(Ref_Size, " -> ", Alt_Size),
    Size_Shift = case_when(
      is.na(Ref_Size) | is.na(Alt_Size) ~ "unknown",
      Ref_Size == Alt_Size ~ "no_size_change",
      TRUE ~ "size_change"
    ),
    
    Polarity_Change = paste0(Ref_Polarity, " -> ", Alt_Polarity),
    Polarity_Shift = case_when(
      is.na(Ref_Polarity) | is.na(Alt_Polarity) ~ "unknown",
      Ref_Polarity == Alt_Polarity ~ "no_polarity_change",
      TRUE ~ "polarity_change"
    )
  )

# View detailed table
missense_props %>%
  select(
    `Gene ID`, Chromosome, Location,
    RefAA, AltAA, AA_Change,
    Ref_Charge, Alt_Charge, Charge_Change, Charge_Shift,
    Ref_Size, Alt_Size, Size_Change, Size_Shift,
    Ref_Polarity, Alt_Polarity, Polarity_Change, Polarity_Shift
  ) %>%
  arrange(`Gene ID`, Chromosome, Location) %>%
  print(n = Inf, width = Inf)

# Summaries
missense_props %>% count(Charge_Shift, sort = TRUE)
missense_props %>% count(Size_Shift, sort = TRUE)
missense_props %>% count(Polarity_Shift, sort = TRUE)

# Optional: exact transition summaries (e.g., positive->neutral)
missense_props %>% count(Charge_Change, sort = TRUE)
missense_props %>% count(Size_Change, sort = TRUE)

# Save
write_csv(missense_props, "gwas_missense_snp_with_charge_size_annotations.csv")

# ==========================================
# 3-category summary for UNIQUE SNPs only
# noncoding / coding-synonymous / coding-missense
# ==========================================

library(dplyr)
library(readr)
# ==========================================
# 3-category summary for UNIQUE SNPs only
# using your FungiDB column names
# ==========================================

snps_unique_base <- snps %>%
  mutate(
    SNP_ID = paste(Chromosome, Location, `Minor Allele`, sep = "_"),
    Coding = tolower(Coding)
  ) %>%
  distinct(SNP_ID, .keep_all = TRUE)

cat("\n=== Total unique SNPs ===\n")
cat(nrow(snps_unique_base), "\n")

# noncoding unique SNPs
noncoding_unique <- snps_unique_base %>%
  filter(Coding != "coding") %>%
  transmute(
    SNP_ID,
    Category3 = "noncoding"
  )

# coding unique SNPs from your codon-classified table
coding_unique <- snps_codons %>%
  mutate(
    SNP_ID = paste(Chromosome, Location, alt, sep = "_")
  ) %>%
  distinct(SNP_ID, .keep_all = TRUE)

# make 3 categories
coding_3cat <- coding_unique %>%
  transmute(
    SNP_ID,
    Category3 = case_when(
      SNP_Class == "synonymous" ~ "coding-synonymous",
      SNP_Class %in% c("missense", "nonsense", "stop_loss", "other", "unknown") ~ "coding-missense"
    )
  )

three_cat_summary <- bind_rows(
  noncoding_unique,
  coding_3cat
) %>%
  count(Category3, name = "n") %>%
  mutate(
    prop = n / sum(n),
    percent = round(prop * 100, 1)
  )

print(three_cat_summary)
cat("\nSum = ", sum(three_cat_summary$n), "\n")

chr2_snps <- snps %>%
  filter(Chromosome == "bcin_chr_2")

chr2_unique <- chr2_snps %>%
  mutate(SNP_ID = paste(Chromosome, Location, `Minor Allele`, sep = "_")) %>%
  distinct(SNP_ID, .keep_all = TRUE)

nrow(chr2_unique)

chr2_unique %>%
  select(`Gene ID`, Chromosome, Location, Coding, `Minor Allele`) %>%
  arrange(`Gene ID`, Location) %>%
  print(n = Inf)
chr2_unique %>%
  count(`Gene ID`, sort = TRUE)
chr2_unique %>%
  mutate(Coding = tolower(Coding)) %>%
  count(Coding, sort = TRUE)


all_chr_summary <- all_chr_unique %>%
  mutate(Coding = tolower(Coding)) %>%
  group_by(Chromosome) %>%
  summarise(
    n_unique_snps = n(),
    n_genes = n_distinct(`Gene ID`, na.rm = TRUE),
    n_coding = sum(Coding == "coding", na.rm = TRUE),
    n_noncoding = sum(Coding != "coding", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_unique_snps))

all_chr_summary



