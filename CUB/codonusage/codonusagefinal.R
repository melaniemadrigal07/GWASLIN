library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)

# cleaning fungidb csv output into bracket format
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

#functions to read in strands and snps based on data, i.e: gene on reverse strand, coding sequence is reverse complement
rev_comp <- function(seq) {
  if (is.na(seq)) return(NA_character_)
  s <- chartr("ACGT", "TGCA", toupper(seq))
  paste(rev(strsplit(s, "")[[1]]), collapse = "")
}
#normalize strand into a clean one
normalize_strand <- function(x) {
  x2 <- tolower(trimws(as.character(x)))
  case_when(
    x2 %in% c("+", "forward", "fwd", "plus", "sense") ~ "+",
    x2 %in% c("-", "reverse", "rev", "minus", "antisense") ~ "-",
    TRUE ~ NA_character_
  )
}

# extract codon sequences based on the if strand on fungidb is reverse or forward
get_codons_from_ctx_strand <- function(ctx, cds_pos, alt_base, strand = "+") {
  if (is.na(ctx) || ctx == "N/A" || is.na(cds_pos) || is.na(alt_base) || is.na(strand)) return(NULL)
  
  alt_base <- toupper(alt_base)
  if (!alt_base %in% c("A","C","G","T")) return(NULL)
#finding codon position  
  pos_in_codon_coding <- ((cds_pos - 1) %% 3) + 1
#mirrors codon position for reverse strand
  pos_in_codon_genomic <- if (strand == "-") 4 - pos_in_codon_coding else pos_in_codon_coding
#extract the codon for the selected SNP  
  codon_raw <- dplyr::case_when(
    pos_in_codon_genomic == 1 ~ stringr::str_extract(ctx, "\\[[ACGT]\\][ACGT]{2}"),
    pos_in_codon_genomic == 2 ~ stringr::str_extract(ctx, "[ACGT]\\[[ACGT]\\][ACGT]"),
    pos_in_codon_genomic == 3 ~ stringr::str_extract(ctx, "[ACGT]{2}\\[[ACGT]\\]"),
    TRUE ~ NA_character_
  )
  if (is.na(codon_raw)) return(NULL)
#strip brackets to get the reference codon  
  refCodon_gen <- stringr::str_replace_all(codon_raw, "\\[|\\]", "")
  if (nchar(refCodon_gen) != 3) return(NULL)
#construct alternate codon by subsitiuting with the alternate  
  v <- strsplit(refCodon_gen, "")[[1]]
  v[pos_in_codon_genomic] <- alt_base
  altCodon_gen <- paste(v, collapse = "")
 #converting coding-strand codons if it is a reverse strand 
  if (strand == "-") {
    refCodon <- rev_comp(refCodon_gen)
    altCodon <- rev_comp(altCodon_gen)
  } else {
    refCodon <- refCodon_gen
    altCodon <- altCodon_gen
  }
  
  tibble(RefCodon = refCodon, AltCodon = altCodon, pos_in_codon = pos_in_codon_coding)
}

#loading in a codon table to pull from later
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
#codes codon into the single amino acid letter
codon_to_aa <- function(codon) unname(codon_table[toupper(codon)])

#load SNP data from fungidb
snps <- read_csv("NgsSnpBySourceId_Summary.csv", show_col_types = FALSE)

# FungiDB effect parsing (SAFE: exact labels) 
snps2 <- snps %>%
  mutate(
    Coding = tolower(Coding),
    fungidb_raw = tolower(trimws(as.character(`Has Non-synonymous`))),
    FungiDB_Effect = case_when(
      Coding != "coding" ~ NA_character_,
      fungidb_raw == "synonymous" ~ "synonymous",
      fungidb_raw == "has non-synonymous" ~ "nonsynonymous",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-fungidb_raw)

# Add codons/AA using Minor Allele
snps_codons <- snps2 %>%
  mutate(
    cds_pos = parse_integer(as.character(na_if(`Position in CDS`, "N/A"))),
    strand  = normalize_strand(`Gene strand`),     # adjust name if needed
    ctx     = clean_context_bracket(`SNP context`),
    alt     = str_extract(toupper(as.character(`Minor Allele`)), "[ACGT]")
  ) %>%
  filter(Coding == "coding", !is.na(FungiDB_Effect), !is.na(cds_pos), !is.na(strand)) %>%
  rowwise() %>%
  mutate(tmp = list(get_codons_from_ctx_strand(ctx, cds_pos, alt, strand))) %>%
  ungroup() %>%
  unnest_wider(tmp) %>%
  mutate(
    RefAA = codon_to_aa(RefCodon),
    AltAA = codon_to_aa(AltCodon)
  )

# keep only synonymous SNPS to look at codon usage
fungidb_syn_codons <- snps_codons %>%
  filter(FungiDB_Effect == "synonymous") %>%
  select(`Gene ID`, Chromosome, Location, `Gene strand`, cds_pos,
         FungiDB_Effect, RefCodon, AltCodon, RefAA, AltAA, ctx)

# check that amino acids are preserved
fungidb_syn_codons %>%
  summarise(
    n = n(),
    all_AA_same = all(RefAA == AltAA),
    n_mismatch = sum(RefAA != AltAA)
  ) %>% print()

fungidb_syn_codons

#checking if bias changes
codon_usage <- read_csv("b.cincodonusage.csv")


#joining codon usage with snp data
syn_with_usage <- fungidb_syn_codons %>%
  mutate(
    RefCodon = toupper(RefCodon),
    AltCodon = toupper(AltCodon)
  ) %>%
  left_join(
    codon_usage %>%
      select(Codon, Usage_B_cinerea) %>%
      rename(RefCodon = Codon, RefUsage = Usage_B_cinerea),
    by = "RefCodon"
  ) %>%
  left_join(
    codon_usage %>%
      select(Codon, Usage_B_cinerea) %>%
      rename(AltCodon = Codon, AltUsage = Usage_B_cinerea),
    by = "AltCodon"
  ) %>%
  mutate(
    DeltaUsage = AltUsage - RefUsage,
    TowardPreferred = case_when(
      DeltaUsage > 0 ~ "toward_preferred",
      DeltaUsage < 0 ~ "toward_rare",
      TRUE ~ "no_change"
    )
  )
#checking that join was successful
syn_with_usage %>%
  summarise(
    n = n(),
    missing_ref = sum(is.na(RefUsage)),
    missing_alt = sum(is.na(AltUsage))
  )
#view ranked table
syn_with_usage %>%
  select(`Gene ID`, Chromosome, Location, `Gene strand`,
         RefCodon, AltCodon, RefUsage, AltUsage, DeltaUsage, TowardPreferred) %>%
  arrange(desc(DeltaUsage)) %>%
  print(n = Inf, width = Inf)
#count how many SNPs go each way
syn_with_usage %>%
  count(TowardPreferred, sort = TRUE)
#calculate descriptive statistics
syn_with_usage %>%
  summarise(
    n = n(),
    mean_delta = mean(DeltaUsage),
    median_delta = median(DeltaUsage),
    prop_toward_preferred = mean(DeltaUsage > 0),
    prop_toward_rare = mean(DeltaUsage < 0)
  )


#making dot+density plot 
ggplot(syn_with_usage, aes(x = DeltaUsage)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_density(fill = "grey80", alpha = 0.6) +
  geom_jitter(y = 0, height = 0.02, size = 2, alpha = 0.9) +
  theme_classic(base_size = 14) +
  labs(
    x = expression(Delta~Codon~Usage~"(" * Alt - Ref * ", % )"),
    y = NULL,
    title = "Direction of codon-usage change for synonymous GWAS SNPs"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
#signed bar plot
cb_palette <- c(
  toward_preferred = "#0072B2",  # blue
  toward_rare      = "#D55E00",  # vermillion
  no_change        = "#999999"   # neutral grey
)

ggplot(syn_with_usage, aes(
  y = reorder(`Gene ID`, DeltaUsage),
  x = DeltaUsage,
  fill = TowardPreferred
)) +
  geom_col(width = 0.55) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      toward_preferred = "#0072B2",
      toward_rare      = "#D55E00",
      no_change        = "#B0B0B0"
    )
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Synonymous GWAS SNPs shift codon usage bias",
    x = expression(Delta~Codon~Usage~"(" * Alt - Ref * ", % )"),
    y = "Gene",
    fill = "Direction of change"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t = 8)),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )




#output results
write_csv(syn_with_usage, "synonymous_snps_with_codon_usage_bias.csv")


#what codons does it shift 2
rare_shifts <- syn_with_usage %>%
  filter(DeltaUsage < 0) %>%
  select(
    `Gene ID`,
    Chromosome,
    Location,
    `Gene strand`,
    RefCodon,
    AltCodon,
    RefUsage,
    AltUsage,
    DeltaUsage
  ) %>%
  arrange(DeltaUsage)
print(rare_shifts, n = Inf, width = Inf)
#summarize what appears
rare_codon_counts <- rare_shifts %>%
  count(AltCodon, sort = TRUE) %>%
  rename(n_snps = n)

rare_codon_counts

#output csv files
write_csv(rare_shifts, "synonymous_snps_shifting_to_rare_codons.csv")
write_csv(rare_codon_counts, "rare_codon_targets_summary.csv")

dir_2cat <- syn_with_usage %>%
  filter(TowardPreferred %in% c("toward_preferred", "toward_rare")) %>%
  count(TowardPreferred)

dir_2cat

# H0: equal frequency (50:50)
chisq.test(dir_2cat$n, p = c(0.5, 0.5))
