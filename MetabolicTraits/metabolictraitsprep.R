# ============================================================
# FINAL SCRIPT — Lab color metrics (ΔE + b*) traits (TP1–5)
#
# This script:
# - Reads LAB_Means*.csv across plates
# - Reads plate maps and assigns isolate IDs
# - Corrects known plate issue (plate 3, F9)
# - Applies background correction (M+R + Solvent controls)
# - Converts discrete timepoints into real experimental time
#   using plate-specific start times (accounts for ~26 min delay)
# - Computes ΔE (color change magnitude)
# - Extracts traits (early response, TP1–5)
# - Outputs timecourse + trait tables
# ============================================================

library(tidyverse)
library(lubridate)
library(stringr)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)

# ---------------------------
# SETTINGS
# only using early timepoints (first 5)
TP_MAX <- 5

# define control wells used for background correction
controls <- c("M+R", "Solvent_1", "Solvent_2", "Solvent_3")

# regex patterns to identify background wells
background_patterns <- c("M\\+R", "^Solvent_")

# ---------------------------
# FILES
# measurement files (plate outputs)
measurement_files <- list(
  list(file = "LAB_Means2.csv", plate_id = "2", layout_id = "2"),
  list(file = "LAB_Means1.csv", plate_id = "1", layout_id = "2"),
  list(file = "LAB_Means3.csv", plate_id = "3", layout_id = "2"),
  list(file = "LAB_Means4.csv", plate_id = "4", layout_id = "2"),
  list(file = "LAB_Means5.csv", plate_id = "5", layout_id = "1"),
  list(file = "LAB_Means6.csv", plate_id = "6", layout_id = "1"),
  list(file = "LAB_Means7.csv", plate_id = "7", layout_id = "1"),
  list(file = "LAB_Means8.csv", plate_id = "8", layout_id = "1")
)

# plate layouts (maps wells → isolate IDs)
plate_map_files <- list(
  list(file = "Plate_1_with_Controls.csv", layout_id = "1"),
  list(file = "Plate_2_with_Controls.csv", layout_id = "2")
)

# optional name key (maps isolate number → name)
name_key_file <- "NameKey.csv"

# ---------------------------
# HELPERS

# read Cytation output and tag plate + layout
read_lab_measurements <- function(file, plate_id, layout_id) {
  read_csv(file, show_col_types = FALSE) %>%
    mutate(
      plateid  = as.character(plate_id),
      layoutid = as.character(layout_id)
    )
}

# convert plate map format into long format (Well → Sample)
read_plate_map <- function(file, layout_id) {
  read_csv(file, skip = 1, col_names = FALSE, show_col_types = FALSE) %>%
    filter(str_detect(X1, "^[A-H]$")) %>%
    setNames(c("Row", as.character(1:12))) %>%
    mutate(across(`1`:`12`, as.character)) %>%
    pivot_longer(cols = -Row, names_to = "Column", values_to = "Sample") %>%
    mutate(
      Well     = paste0(Row, Column),
      layoutid = as.character(layout_id)
    ) %>%
    select(layoutid, Well, Sample)
}

# trapezoid rule for AUC calculation
trapz_vec <- function(x, y) {
  o <- order(x); x <- as.numeric(x[o]); y <- as.numeric(y[o])
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# ---------------------------
# 1) READ DATA

# combine all measurement files
combined_data <- map_dfr(measurement_files, ~ read_lab_measurements(.x$file, .x$plate_id, .x$layout_id))

# combine plate maps
plate_maps <- bind_rows(
  purrr::map_dfr(plate_map_files, ~ read_plate_map(.x$file, .x$layout_id))
)

# ---------------------------
# 2) JOIN SAMPLE IDS + FIX KNOWN ISSUE

combined_data <- combined_data %>%
  left_join(plate_maps, by = c("layoutid", "Well")) %>%
  mutate(
    # fix incorrect/missing sample in plate 3 F9
    Sample = if_else(plateid == "3" & Well == "F9", "23", Sample)
  )

# ---------------------------
# 3) ADD REAL TIME (IMPORTANT)

# define actual start times for each plate read
# this accounts for the stagger (~26 min delay between plates)
read1_start_times <- c(
  "1" = ymd_hms("2025-05-14 19:36:00"),
  "2" = ymd_hms("2025-05-14 20:03:00"),
  "3" = ymd_hms("2025-05-14 20:31:00"),
  "4" = ymd_hms("2025-05-14 20:58:00"),
  "5" = ymd_hms("2025-05-14 21:25:00"),
  "6" = ymd_hms("2025-05-14 21:52:00"),
  "7" = ymd_hms("2025-05-14 22:16:00"),
  "8" = ymd_hms("2025-05-14 22:47:00")
)

combined_data <- combined_data %>%
  mutate(
    Read1_Start = read1_start_times[plateid],
    
    # convert discrete timepoints → real time (4 hr spacing)
    ReadTime = Read1_Start + minutes((Timepoint - 1) * 240)
  )

# ---------------------------
# 4) BACKGROUND CORRECTION

combined_norm <- combined_data %>%
  group_by(plateid, Timepoint) %>%
  mutate(
    # identify background wells
    is_bg_sample = (
      str_detect(Sample, str_c(background_patterns, collapse = "|")) |
        (plateid == "3" & Well == "F9")
    ) & !(plateid == "3" & Well == "F10"),
    
    # compute background per plate × timepoint
    BG_L = mean(Mean_L[is_bg_sample], na.rm = TRUE),
    BG_a = mean(Mean_a[is_bg_sample], na.rm = TRUE),
    BG_b = mean(Mean_b[is_bg_sample], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    # subtract background
    Norm_L = Mean_L - BG_L,
    Norm_a = Mean_a - BG_a,
    Norm_b = Mean_b - BG_b
  )

# ---------------------------
# 5) CREATE SHARED TIME AXIS

combined_norm <- combined_norm %>%
  mutate(
    # align all plates onto one timeline
    HoursSinceStart = as.numeric(
      difftime(ReadTime, min(ReadTime, na.rm = TRUE), units = "hours")
    )
  )

# ---------------------------
# 6) SAMPLE × TIME MEANS + ΔE

lab_means <- combined_norm %>%
  group_by(Sample, Timepoint, HoursSinceStart) %>%
  summarise(
    Avg_Lc = mean(Norm_L, na.rm = TRUE),
    Avg_ac = mean(Norm_a, na.rm = TRUE),
    Avg_bc = mean(Norm_b, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # ΔE = magnitude of color change
    DeltaE_c = sqrt(Avg_Lc^2 + Avg_ac^2 + Avg_bc^2)
  ) %>%
  arrange(Sample, HoursSinceStart)

# save full timecourse
write_csv(lab_means, "deltaE_timecourse_real_time.csv")

# ---------------------------
# 7) TRAIT EXTRACTION (TP1–5 ONLY)

traits_early <- lab_means %>%
  filter(Timepoint <= TP_MAX) %>%
  filter(!Sample %in% controls) %>%
  group_by(Sample) %>%
  summarise(
    
    # ΔE traits
    max_DeltaE_5 = max(DeltaE_c, na.rm = TRUE),
    
    # area under curve (real time)
    auc_5 = trapz_vec(HoursSinceStart, DeltaE_c),
    
    # slope over time
    slope_5 = {
      ok <- is.finite(HoursSinceStart) & is.finite(DeltaE_c)
      x <- HoursSinceStart[ok]; y <- DeltaE_c[ok]
      if (length(x) >= 2) coef(lm(y ~ x))[2] else NA_real_
    },
    
    # b* traits
    max_b    = max(Avg_bc, na.rm = TRUE),
    
    # time to peak (now in HOURS, not timepoint index)
    tp_max_b = HoursSinceStart[which.max(Avg_bc)][1],
    
    auc_b = trapz_vec(HoursSinceStart, Avg_bc),
    
    slope_b_5 = {
      ok <- is.finite(HoursSinceStart) & is.finite(Avg_bc)
      x <- HoursSinceStart[ok]; y <- Avg_bc[ok]
      if (length(x) >= 2) coef(lm(y ~ x))[2] else NA_real_
    },
    
    r2_b_5 = {
      ok <- is.finite(HoursSinceStart) & is.finite(Avg_bc)
      x <- HoursSinceStart[ok]; y <- Avg_bc[ok]
      if (length(x) >= 2) summary(lm(y ~ x))$r.squared else NA_real_
    },
    
    range_b = max(Avg_bc, na.rm = TRUE) - min(Avg_bc, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(Sample)

# save traits
write_csv(traits_early, "deltaE_bstar_traits_TP1-5_real_time.csv")

# ---------------------------
# 8)NAME KEY JOIN

if (file.exists(name_key_file)) {
  name_key <- read_csv(name_key_file, show_col_types = FALSE) %>%
    rename(
      SampleName = Isolate,
      Sample     = `Isolate Number`
    ) %>%
    mutate(Sample = as.character(Sample))
  
  traits_named <- traits_early %>%
    mutate(Sample = as.character(Sample)) %>%
    left_join(name_key, by = "Sample") %>%
    relocate(SampleName, .after = Sample) %>%
    arrange(Sample)
  
  write_csv(traits_named, "deltaE_traits_with_namekey_real_time.csv")
}