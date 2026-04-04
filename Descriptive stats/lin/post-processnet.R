install.packages("devtools")
devtools::install_github("melaniemadrigal07/SkelPyR")
library(SkelPyR)
#set wd
setwd("/Volumes/Zayn/linkskel/Calculations8/skelpy_output")

results <- run_fyskel_pipeline()

#building csv
library(dplyr)
library(readr)
library(stringr)

parent_dir <- "/Volumes/Zayn/linkskel"

calc_dirs <- list.dirs(parent_dir, recursive = FALSE, full.names = TRUE)
calc_dirs <- calc_dirs[grepl("^Calculation", basename(calc_dirs))]

network_all <- bind_rows(lapply(calc_dirs, function(dir) {
  f <- file.path(dir, "skelpy_output", "hyphae_network_summary.csv")
  if (!file.exists(f)) return(NULL)
  
  plate_label <- basename(dir)                 # "Calculations2"
  plate_id    <- str_extract(plate_label, "\\d+")  # "2"
  
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      PlateID     = plate_id,
      PlateLabel  = plate_label
    )
}))
# What plates did we extract?
unique(network_all$PlateID)

# Counts per plate
table(network_all$PlateID)

# Quick peek
network_all %>% select(PlateID, PlateLabel) %>% distinct()

library(tidyverse)

# --- A) Read plate-map templates (LayoutID + Well -> Sample) ---
read_plate_map <- function(file, layout_id) {
  read_csv(file, skip = 1, col_names = FALSE, show_col_types = FALSE) %>%
    filter(str_detect(X1, "^[A-H]$")) %>%
    setNames(c("Row", as.character(1:12))) %>%
    mutate(across(`1`:`12`, as.character)) %>%
    pivot_longer(cols = -Row, names_to = "Column", values_to = "Sample") %>%
    mutate(
      Well = paste0(Row, Column),
      LayoutID = as.character(layout_id)
    ) %>%
    select(LayoutID, Well, Sample)
}

plate_map_all <- bind_rows(
  read_plate_map("Plate_1_with_Controls.csv", layout_id = "1"),
  read_plate_map("Plate_2_with_Controls.csv", layout_id = "2")
)

# --- B) Plate -> layout mapping (edit if yours differs) ---
plate_layout_key <- tibble(
  PlateID  = as.character(1:8),
  LayoutID = if_else(as.integer(PlateID) <= 4, "2", "1")
)

# Expand to per-plate map: PlateID + Well -> Sample
plate_map_by_plate <- plate_layout_key %>%
  left_join(plate_map_all, by = "LayoutID", relationship = "many-to-many") %>%
  select(PlateID, Well, Sample)

# Extract Well from SkelPy output (you may need to adjust column name) ---
# Try common possibilities in order:
candidate_cols <- c("Filename", "filename", "image", "Image", "file", "File")

col_in_data <- candidate_cols[candidate_cols %in% names(network_all)][1]

if (is.na(col_in_data)) {
  stop("Couldn't find a filename-like column in network_all. Run colnames(network_all) and tell me what looks like the file/well ID.")
}

network_all <- network_all %>%
  mutate(
    Well = str_extract(.data[[col_in_data]], "(?<![A-Z])[A-H](?:[1-9]|1[0-2])(?!\\d)")
  )

# --- D) Join sample labels onto SkelPy metrics ---
network_annotated <- network_all %>%
  left_join(plate_map_by_plate, by = c("PlateID", "Well"))

# --- E) Sanity checks ---
table(is.na(network_annotated$Well))     # should be mostly FALSE
table(is.na(network_annotated$Sample))   # should be mostly FALSE (except blanks/unused wells)

# --- F)averages by PlateID + Sample ---
means_by_plate_sample <- network_annotated %>%
  group_by(PlateID, Sample) %>%
  summarise(
    n = n(),
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Write outputs if you want
write_csv(network_annotated, file.path(parent_dir, "ALL_network_metrics_with_plate_and_sample.csv"))
write_csv(means_by_plate_sample, file.path(parent_dir, "MEANS_network_metrics_by_plate_sample.csv"))
#convert samples to names
library(dplyr)
library(readr)

name_key <- read_csv("NameKey.csv", show_col_types = FALSE) %>%
  rename(
    IsolateNumber = `Isolate Number`,
    IsolateName   = Isolate
  ) %>%
  mutate(IsolateNumber = as.character(IsolateNumber))

# sanity check
head(name_key)
colnames(name_key)

network_annotated <- network_annotated %>%
  mutate(Sample = as.character(Sample)) %>%
  left_join(
    name_key,
    by = c("Sample" = "IsolateNumber")
  ) %>%
  mutate(
    Isolate = IsolateName
  ) %>%
  select(-IsolateName)


means_by_plate_sample <- means_by_plate_sample %>%
  mutate(Sample = as.character(Sample)) %>%
  left_join(
    name_key,
    by = c("Sample" = "IsolateNumber")
  ) %>%
  mutate(Isolate = IsolateName) %>%
  select(-IsolateName)

write_csv(
  network_annotated,
  file.path(parent_dir, "ALL_network_metrics_with_plate_well_isolate.csv")
)


write_csv(
  means_by_plate_sample,
  file.path(parent_dir, "MEANS_network_metrics_by_plate_isolate.csv")
)

means_by_isolate <- means_by_plate_sample %>%
  group_by(Isolate) %>%
  summarise(
    n_plates = n(),
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

write_csv(
  means_by_isolate,
  file.path(parent_dir, "MEANS_network_metrics_by_isolate.csv")
)

traits_gwas <- means_by_isolate_w %>%
  transmute(
    Isolate,
    total_length_mean,          # growth phenotype
    tip_fraction_resid,         # size-independent allocation phenotype
    tip_fraction_weighted,      # optional: keep for reference (don’t GWAS if confounded)
    tip_fraction_unweighted     # optional
  )
write_csv(
  traits_gwas,
  file.path(parent_dir, "GWAS_traits_isolate_level.csv")
)
