source("functions/prepare_env.R")
PrepareEnv()
sc <- reticulate::import("scanpy")

log_message("Start loading data...")
data_dir <- "../../data/BrainOmicsData/raw/"
res_dir <- "../../data/BrainOmicsData/processed/"

data_record <- read.csv("data/dataset_used.csv")
accessions <- unique(data_record$Accession)

success_list <- c()
fail_list <- c()

for (acc in accessions) {
  log_message("Processing {.val {acc}}")

  raw_path <- file.path(data_dir, acc)
  processed_path <- file.path(res_dir, acc)
  output_file <- file.path(processed_path, paste0(acc, ".rds"))

  if (file.exists(output_file)) {
    log_message("Skipping {.val {acc}} - already processed")
    success_list <- c(success_list, acc)
    next
  }

  if (!dir.exists(raw_path)) {
    log_message("No data directory found for {.val {acc}}")
    fail_list <- c(fail_list, acc)
    next
  }

  h5ad_files <- list.files(
    raw_path,
    pattern = "\\.h5ad$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(h5ad_files) == 0) {
    log_message("No h5ad file found for {.val {acc}}")
    fail_list <- c(fail_list, acc)
    next
  }

  h5ad_file <- h5ad_files[1]

  tryCatch(
    {
      adata <- sc$read_h5ad(h5ad_file)
      srt <- adata_to_srt(adata)
      dir.create(processed_path, recursive = TRUE, showWarnings = FALSE)
      saveRDS(srt, output_file)
      log_message("Successfully processed {.val {acc}}")
      success_list <- c(success_list, acc)
    },
    error = function(e) {
      log_message("Failed to process {.val {acc}} : {.error {e$message}}")
      fail_list <<- c(fail_list, acc)
    }
  )
  rm(adata)
  rm(srt)
  gc()
}

log_message(
  "Completed. Success: {.val {length(success_list)}} Failed: {.val {length(fail_list)}}"
)
write.csv(
  success_list,
  file.path(res_dir, "BCAtlas/success_accessions.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  fail_list,
  file.path(res_dir, "BCAtlas/failed_accessions.csv"),
  row.names = FALSE,
  quote = FALSE
)

if (length(fail_list) > 0) {
  log_message(
    "Failed accessions: {.val {fail_list}}"
  )
}
