library(dplyr)
library(tidyverse)
library(data.table)
library(topr)
library(ggplot2)
library(ComplexUpset)
library(qqman)

prefs <- read.table(
  "Scripts/### Preferences ###",
  header = FALSE,
  sep = "=",
  skip = 1
)

SNPset <- as.character(prefs[2, 2])
pheno.name <- as.character(prefs[1, 2])
multcomp <- 3870
ploidy <- 1
suggthresh <- 0.001

envs <- as.character(read.table("environments_to_run.txt")[, 1])
traits <- as.character(read.table("traits_to_run.txt")[, 1])

Genome_annotation <- fread("Software/parsed_geneannotation_output.tsv")

header <- readLines(paste0("data/", pheno.name), n = 1)
n_cols <- length(strsplit(header, ",")[[1]])
col_classes <- c("character", rep("", n_cols - 1))

pheno.data <- fread(
  paste0("data/", pheno.name),
  colClasses = col_classes
)

models <- data.frame(
  Model = c("K", "PK"),
  AssocDir = c("Tables/Assoc_files_K", "Tables/Assoc_files_PK")
)

AllGenomewideSNPS <- data.frame()
AllSuggestiveSNPS <- data.frame()
LambdaTable <- data.frame()

for (m in 1:nrow(models)) {
  
  model_name <- models$Model[m]
  assoc_dir <- models$AssocDir[m]
  
  for (q in 1:length(envs)) {
    for (i in 1:length(traits)) {
      
      trait <- traits[i]
      env <- envs[q]
      
      message("Running: ", model_name, " | ", trait, "_", env)
      
      assoc_file <- paste0(
        assoc_dir,
        "/",
        trait,
        "_",
        env,
        ".assoc.txt"
      )
      
      if (!file.exists(assoc_file)) {
        message("Skipping missing file: ", assoc_file)
        next
      }
      
      phenotyperange <- range(pheno.data[[trait]], na.rm = TRUE)[2] -
        range(pheno.data[[trait]], na.rm = TRUE)[1]
      
      snips <- fread(assoc_file, header = TRUE)
      snips <- snips[!is.na(snips$beta), ]
      
      snips4manhattan <- data.frame(
        CHROM = snips$chr,
        POS = snips$ps,
        ID = snips$rs,
        REF = NA,
        ALT = NA,
        P = snips$p_wald,
        BETA = snips$beta,
        AF = snips$af,
        SE = snips$se
      ) %>%
        filter(!is.na(P), P > 0, P <= 1)
      
      snips4manhattan$EFsize <- (ploidy * snips4manhattan$BETA) / phenotyperange
      snips4manhattan$SEstandard <- (ploidy * snips4manhattan$SE) / phenotyperange
      
      #-----------------------------
      # Manhattan plot
      #-----------------------------
      dir.create(
        paste0("Plots/Manhattans/", model_name),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      mantoprint <- manhattanExtra(
        snips4manhattan,
        title = paste(model_name, trait, env, sep = "_"),
        genome_wide_thresh = 0.05 / multcomp,
        annotate = 0.05 / multcomp,
        suggestive_thresh = quantile(snips4manhattan$P, suggthresh, na.rm = TRUE),
        flank_size = 1000,
        region_size = 1,
        sign_thresh_label_size = 0.001,
        sign_thresh_color = c("#005cb9", "#ab2328"),
        color = c("grey50", "#005cb9", "#ab2328"),
        build = Genome_annotation
      )
      
      pdf(
        paste0(
          "Plots/Manhattans/",
          model_name,
          "/manhattanplot_",
          model_name,
          "_",
          trait,
          "_",
          env,
          ".pdf"
        ),
        onefile = TRUE
      )
      print(mantoprint)
      dev.off()
      
      #-----------------------------
      # QQ plot + lambda
      #-----------------------------
      dir.create(
        paste0("Plots/QQplots/", model_name),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      qq_pvalues <- snips4manhattan$P
      
      chisq <- qchisq(1 - qq_pvalues, 1)
      lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
      
      LambdaTable <- rbind(
        LambdaTable,
        data.frame(
          Model = model_name,
          Trait = trait,
          Environment = env,
          Lambda = lambda
        )
      )
      
      pdf(
        paste0(
          "Plots/QQplots/",
          model_name,
          "/qqplot_",
          model_name,
          "_",
          trait,
          "_",
          env,
          ".pdf"
        )
      )
      
      qqman::qq(
        qq_pvalues,
        main = paste0(
          "QQ Plot: ",
          model_name,
          "_",
          trait,
          "_",
          env,
          "\nLambda = ",
          round(lambda, 3)
        )
      )
      
      dev.off()
      
      #-----------------------------
      # Significant/suggestive SNPs
      #-----------------------------
      SNPSforregion_table <- get_sign_and_sugg_loci(
        snips4manhattan,
        suggestive_thresh = quantile(snips4manhattan$P, suggthresh, na.rm = TRUE),
        genome_wide_thresh = 0.05 / multcomp,
        flank_size = 1,
        region_size = 1
      )
      
      if (dim(SNPSforregion_table$genome_wide_snps)[1] > 0) {
        tmp <- SNPSforregion_table$genome_wide_snps
        tmp$Trait <- trait
        tmp$Environment <- env
        tmp$Model <- model_name
        
        AllGenomewideSNPS <- rbind(AllGenomewideSNPS, tmp)
      }
      
      if (dim(SNPSforregion_table$suggestive_snps)[1] > 0) {
        tmp <- SNPSforregion_table$suggestive_snps
        tmp$Trait <- trait
        tmp$Environment <- env
        tmp$Model <- model_name
        
        AllSuggestiveSNPS <- rbind(AllSuggestiveSNPS, tmp)
      }
    }
  }
}

# Save summary tables
dir.create("Tables", recursive = TRUE, showWarnings = FALSE)

write.csv(
  AllGenomewideSNPS,
  "Tables/GenomewideSNPS_K_vs_PK.csv",
  row.names = FALSE
)

write.csv(
  AllSuggestiveSNPS,
  "Tables/SuggestiveSNPS_K_vs_PK.csv",
  row.names = FALSE
)

write.csv(
  LambdaTable,
  "Tables/LambdaSummary_K_vs_PK.csv",
  row.names = FALSE
)

#-----------------------------
# UpSet plot: SNPs shared by K vs PK
#-----------------------------
dir.create("Plots/Colocalization", recursive = TRUE, showWarnings = FALSE)

if (nrow(AllGenomewideSNPS) > 0) {
  
  wide_model_df <- AllGenomewideSNPS %>%
    select(ID, Model) %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(
      names_from = Model,
      values_from = value,
      values_fill = 0
    )
  
  model_sets <- colnames(wide_model_df)[-1]
  
  pdf("Plots/Colocalization/Upset_Genomewide_K_vs_PK.pdf")
  
  ComplexUpset::upset(
    wide_model_df,
    intersect = model_sets,
    name = "Model",
    width_ratio = 0.25,
    base_annotations = list(
      "Intersection size" = intersection_size(
        text = list(size = 3)
      )
    ),
    set_sizes = upset_set_size()
  )
  
  dev.off()
}

#-----------------------------
# UpSet plot: SNPs shared by Trait × Model
#-----------------------------
if (nrow(AllGenomewideSNPS) > 0) {
  
  wide_trait_model_df <- AllGenomewideSNPS %>%
    mutate(ModelTrait = paste(Model, Trait, sep = "_")) %>%
    select(ID, ModelTrait) %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(
      names_from = ModelTrait,
      values_from = value,
      values_fill = 0
    )
  
  trait_model_sets <- colnames(wide_trait_model_df)[-1]
  
  pdf(
    "Plots/Colocalization/Upset_Genomewide_TraitModel_K_vs_PK.pdf",
    width = 12,
    height = 8
  )
  
  ComplexUpset::upset(
    wide_trait_model_df,
    intersect = trait_model_sets,
    name = "Trait_Model",
    width_ratio = 0.2,
    base_annotations = list(
      "Intersection size" = intersection_size(
        text = list(size = 3)
      )
    ),
    set_sizes = upset_set_size()
  )
  
  dev.off()
}

#-----------------------------
# Combined QQ plot: K vs PK

dir.create(
  "Plots/QQplots_Combined",
  recursive = TRUE,
  showWarnings = FALSE
)

# Read K model
assoc_file_K <- paste0(
  "Tables/Assoc_files_K/",
  trait,
  "_",
  env,
  ".assoc.txt"
)

# Read PK model
assoc_file_PK <- paste0(
  "Tables/Assoc_files_PK/",
  trait,
  "_",
  env,
  ".assoc.txt"
)

if (file.exists(assoc_file_K) & file.exists(assoc_file_PK)) {
  
  snips_K <- fread(assoc_file_K)
  snips_PK <- fread(assoc_file_PK)
  
  pvals_K <- snips_K$p_wald[
    !is.na(snips_K$p_wald) &
      snips_K$p_wald > 0 &
      snips_K$p_wald <= 1
  ]
  
  pvals_PK <- snips_PK$p_wald[
    !is.na(snips_PK$p_wald) &
      snips_PK$p_wald > 0 &
      snips_PK$p_wald <= 1
  ]
  
  # expected values
  expected_K <- -log10(ppoints(length(pvals_K)))
  observed_K <- -log10(sort(pvals_K))
  
  expected_PK <- -log10(ppoints(length(pvals_PK)))
  observed_PK <- -log10(sort(pvals_PK))
  
  qq_df <- rbind(
    data.frame(
      Expected = expected_K,
      Observed = observed_K,
      Model = "K"
    ),
    data.frame(
      Expected = expected_PK,
      Observed = observed_PK,
      Model = "P+K"
    )
  )
  
  # lambda K
  chisq_K <- qchisq(1 - pvals_K, 1)
  lambda_K <- median(chisq_K, na.rm = TRUE) / qchisq(0.5, 1)
  
  # lambda PK
  chisq_PK <- qchisq(1 - pvals_PK, 1)
  lambda_PK <- median(chisq_PK, na.rm = TRUE) / qchisq(0.5, 1)
  
  p_qq <- ggplot(qq_df, aes(Expected, Observed, color = Model)) +
    geom_point(size = 0.8, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(
        "QQ Plot: ",
        trait,
        "_",
        env,
        "\nK λ = ",
        round(lambda_K, 3),
        " | P+K λ = ",
        round(lambda_PK, 3)
      ),
      x = expression(Expected~~-log[10](p)),
      y = expression(Observed~~-log[10](p))
    ) +
    scale_color_manual(values = c("K" = "#d95f02", "P+K" = "#1b9e77"))
  
  ggsave(
    paste0(
      "Plots/QQplots_Combined/QQ_",
      trait,
      "_",
      env,
      "_K_vs_PK.pdf"
    ),
    p_qq,
    width = 6,
    height = 5
  )
}






