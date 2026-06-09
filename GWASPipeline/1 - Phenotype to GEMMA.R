##### Botrytis Kliebenstein population GWAS using GEMMA -
#####   a re imagined pipeline borrowing some ideas from the version in Temme et al. 2021 & Masalia et al 2018 (PLOS).
##### Genome: HA412HO
##### SNPs: Bo510 SNP set built by Susanna Atwell
##### Algorithm: GEMMA
library(data.table)


setwd("GWAS_pipeline-master")
#### read in preferences
prefs<-read.table("Scripts/### Preferences ###",header=F,sep="=",skip=1)
SNPset<-as.character(prefs[2,2])
pheno.name<-as.character(prefs[1,2])
multcomp<-as.numeric(as.character(prefs[3,2]))

#translate Isolate Names using the key (ignore this section if already run)
#PhenotypeDataset<-read.csv(paste0("data/",pheno.name),header = T)  
#IsolateKey<-read.csv("Software/IsolateKey.csv",header = T)
#PhenotypeDataset$Isolate <- IsolateKey$Isolate[match(PhenotypeDataset$Isolate, IsolateKey$Sample)]
#write.csv(PhenotypeDataset,paste0("data/",pheno.name),row.names = F)

## Read in traits and environments to run
## Read in traits and environments to run
envs <- as.character(read.table("environments_to_run.txt")[, 1])
traits <- as.character(read.table("traits_to_run.txt")[, 1])

pheno.data <- fread(paste("data/", pheno.name, sep = ""))

library(data.table)
library(ggplot2)

# Calculate/load population structure PCA using LD-pruned SNPs

pca_prefix <- paste0("Software/", SNPset, "_LDpruned_PCA")
pca_file <- paste0(pca_prefix, ".eigenvec")

prune_prefix <- paste0("Software/", SNPset, "_LDpruned")
prune_file <- paste0(prune_prefix, ".prune.in")

covar_file_full <- paste0("Software/", SNPset, "_LDpruned_PCA_covariates.txt")
covar_file_gemma <- paste0(SNPset, "_LDpruned_PCA_covariates.txt")

# Step 1: LD pruning
if (!file.exists(prune_file)) {
  
  system(paste(
    "./Software/plink",
    "--bfile", paste0("Software/", SNPset),
    "--indep-pairwise 50'kb' 10 0.1",
    "--out", prune_prefix
  ))
  
  if (!file.exists(prune_file)) {
    stop("LD pruning failed. Could not find: ", prune_file)
  }
}

# Step 2: Run PCA on LD-pruned SNPs
if (!file.exists(pca_file)) {
  
  system(paste(
    "./Software/plink",
    "--bfile", paste0("Software/", SNPset),
    "--extract", prune_file,
    "--pca 10",
    "--out", pca_prefix
  ))
  
  if (!file.exists(pca_file)) {
    stop("PLINK PCA failed. Could not find: ", pca_file)
  }
}

# Step 3: Load PCA file
pca <- fread(pca_file, header = FALSE)
# Step 3b: Load eigenvalues and make elbow/scree plot

eigenval_file <- paste0(pca_prefix, ".eigenval")
eigenvals <- fread(eigenval_file, header = FALSE)

scree_df <- data.frame(
  PC = seq_len(nrow(eigenvals)),
  Eigenvalue = eigenvals$V1
)

scree_df$Variance_Explained <- scree_df$Eigenvalue / sum(scree_df$Eigenvalue)
scree_df$Cumulative_Variance <- cumsum(scree_df$Variance_Explained)

p_elbow <- ggplot(scree_df, aes(x = PC, y = Variance_Explained)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks = scree_df$PC) +
  labs(
    title = "PCA elbow plot",
    x = "Principal component",
    y = "Proportion variance explained"
  )

print(p_elbow)

ggsave(
  "Plots/PCA/PCA_LDpruned_elbow_plot.pdf",
  p_elbow,
  width = 6,
  height = 5
)







# Step 4: Make GEMMA covariate file using PC1-PC3
if (!file.exists(covar_file_full)) {
  
  covar <- pca[, .(V3, V4, V5)]
  
  write.table(
    covar,
    file = covar_file_full,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

# Step 5: Plot PCA
pca_plot_df <- data.frame(
  Sample = pca$V2,
  PC1 = pca$V3,
  PC2 = pca$V4,
  PC3 = pca$V5
)

dir.create("Plots/PCA", recursive = TRUE, showWarnings = FALSE)

p1 <- ggplot(pca_plot_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = "Population structure PCA using LD-pruned SNPs",
    x = "PC1",
    y = "PC2"
  )

print(p1)

ggsave(
  "Plots/PCA/PCA_LDpruned_PC1_PC2.pdf",
  p1,
  width = 6,
  height = 5
)
#-----------------------------
# Move into Software folder to run GEMMA

setwd("Software")



# Move into Software folder to run GEMMA
setwd("Software")

# Create output folders
dir.create("../Tables/Assoc_files_K", recursive = TRUE, showWarnings = FALSE)
dir.create("../Tables/Assoc_files_PK", recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(envs)) {
  
  env <- envs[i]
  
  for (q in 1:length(traits)) {
    
    trait <- traits[q]
    print(paste(trait, env, sep = "_"))
    
    select_cols <- c("Isolate", trait)
    
    if (!select_cols[2] %in% names(pheno.data)) {
      print("phenotype missing")
      next
    }
    
    trait.data <- pheno.data[, ..select_cols]
    
    fam.file <- fread(paste(SNPset, ".fam", sep = ""))
    
    # erase what was in the trait value column for .fam
    fam.file$V6 <- NULL
    
    # merge based on isolate names
    fam.file <- merge(
      fam.file,
      trait.data,
      by.x = "V1",
      by.y = "Isolate",
      all.x = TRUE
    )
    
    # write new fam file for GEMMA
    write.table(
      file = paste(SNPset, ".fam", sep = ""),
      fam.file,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
    
    #-----------------------------
    # Model 1: K-only
    # Kinship correction only
    #-----------------------------
    system(paste(
      "./gemma -bfile ", SNPset,
      " -k ", SNPset, ".cXX.txt",
      " -lmm 1",
      " -outdir ../Tables/Assoc_files_K/",
      " -o ", paste(trait, env, sep = "_"),
      sep = ""
    ))
    
    #-----------------------------
    # Model 2: P + K
    # PCA  + kinship
    #-----------------------------
    system(paste(
      "./gemma -bfile ", SNPset,
      " -k ", SNPset, ".cXX.txt",
      " -c ", covar_file_gemma,
      " -lmm 1",
      " -outdir ../Tables/Assoc_files_PK/",
      " -o ", paste(trait, env, sep = "_"),
      sep = ""
    ))
  }
}

setwd("..")
