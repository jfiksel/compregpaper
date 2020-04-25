### Simulations for comparison of estimates of E[Y|X]
files_dir <- "cluster_output_compare_preds"
files <- list.files(files_dir, full.names = TRUE)
results <- do.call(rbind, lapply(files, readRDS))
output_file <- "comparePredsResults.rds"
saveRDS(results, output_file)

### Simulations for estimation of just M
files_dir <- "cluster_output_sims_no_covar"
files <- list.files(files_dir, full.names = TRUE)
results <- lapply(files, readRDS)
output_file <- "noCovarResults.rds"
saveRDS(results, output_file)
