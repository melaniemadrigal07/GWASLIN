# Data Types

This repository includes scripts that integrate multiple phenotypic and genomic data types to characterize variation in *Botrytis cinerea* responses to Linalool.
## Metabolic Traits (MetabolicTraits/)
Metabolic traits are derived from time-series imaging using LAB color space and reflect physiological and metabolic responses.

Variables include:
- max_DeltaE_5, auc_5, slope_5
- max_b, tp_max_b, auc_b, slope_b_5, r2_b_5, range_b

Notes:
- Background correction performed using M+R controls
- Time-series data summarized into features

## Structural Traits (Structural/)
Structural traits quantify fungal morphology from SkelPy and further skeltonization cleaning from SkelPyR.

Variables include:
- total_length, mean_edge_length, median_edge_length, max_edge_length
- n_edges, n_nodes, n_tips, n_branches
- tip_density, tip_fraction, branch_fraction
- slope_fractalDimension, auc_fractalDimension

Notes:
- Derived from image skeletonization SkelPy
- Capture growth architecture and network organization


## Genomic Data (SNPs)
Genetic variation data used for GWAS.

Includes:
- Chromosome position
- Reference and alternate alleles
- SNP classification (synonymous, non-synonymous, non-coding)

Derived annotations:
- Gene association
- Codon position and amino acid change
- Functional annotation (PFam categories)

## Codon Usage Bias (CUB/)
Analysis of synonymous SNP effects on codon usage.

Includes:
- RefCodon vs AltCodon
- Codon frequency in reference genome
- Codon preference classification

Notes:
- Used to assess effects of synonymous variation on translation efficiency

## Descriptive Statistics (Descriptive stats/)
Summary and exploratory statistics across traits.

Includes:
- Trait distributions
- Summary statistics
- Outlier detection
- Normalization checks

## Integration
These data types are integrated for:
- GWAS (trait–SNP associations)
- PCA and clustering analyses
- Functional interpretation of candidate genes
