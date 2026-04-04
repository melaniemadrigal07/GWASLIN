library(SkelPyR)

parent_dir <- "/Volumes/Zayn/linkskel"

# Find candidate folders
calc_dirs <- list.dirs(parent_dir, recursive = FALSE, full.names = TRUE)
calc_dirs <- calc_dirs[grepl("^Calculations", basename(calc_dirs))]

# Skip folders already completed
needs_run <- !file.exists(
  file.path(calc_dirs, "skelpy_output", "hyphae_network_summary.csv")
)

calc_dirs_to_run <- calc_dirs[needs_run]

message("Folders to process:")
print(basename(calc_dirs_to_run))

# Run pipeline
results <- lapply(calc_dirs_to_run, function(dir) {
  out_dir <- file.path(dir, "skelpy_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Processing ", basename(dir))
  
  tryCatch(
    run_fyskel_pipeline(
      base_dir   = dir,
      output_dir = out_dir
    ),
    error = function(e) {
      message("❌ Failed ", basename(dir), ": ", conditionMessage(e))
      NULL
    }
  )
})
