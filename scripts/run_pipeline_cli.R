invocation_wd <- getwd()
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) return(default)
  args[[hit + 1]]
}
normalize_cli_path <- function(path) {
  if (is.null(path) || grepl("^s3://", path)) return(path)
  if (grepl("^/", path)) return(path)
  normalizePath(file.path(invocation_wd, path), mustWork = FALSE)
}

raw_path <- normalize_cli_path(get_arg("--raw", "data/simulated/raw_plate_reader.csv"))
plate_map_path <- normalize_cli_path(get_arg("--plate-map", "data/simulated/plate_map.csv"))
metadata_path <- normalize_cli_path(get_arg("--metadata", "data/simulated/run_metadata.csv"))
out_dir <- normalize_cli_path(get_arg("--out", file.path("results", paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))))

repo_root <- Sys.getenv("HEKBLUER_HOME", unset = NA_character_)
if (is.na(repo_root) || repo_root == "") {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg)) {
    repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."))
  } else {
    repo_root <- getwd()
  }
}
old_wd <- setwd(repo_root)
on.exit(setwd(old_wd), add = TRUE)

suppressPackageStartupMessages({
  source("R/qc_metrics.R")
  source("R/export.R")
  source("R/analysis.R")
})

raw_data <- read.csv(raw_path, stringsAsFactors = FALSE)
plate_map <- read.csv(plate_map_path, stringsAsFactors = FALSE)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

results <- run_hekblue_analysis(raw_data, plate_map, metadata, output_dir = out_dir)
message("Analysis complete: ", out_dir)
print(results$final_qc_table)
