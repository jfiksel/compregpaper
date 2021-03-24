### Simulations for comparison of estimates of E[Y|X]
files_dir <- "cluster_output_compare_preds"
files <- list.files(files_dir, full.names = TRUE)
results <- do.call(rbind, lapply(files, readRDS))
output_file <- "comparePredsResults.rds"
saveRDS(results, output_file)

### Simulations for comparison of estimates of E[Y|X] for revision
files_dir <- "cluster_output_compare_preds_revision"
files <- list.files(files_dir, full.names = TRUE)
results <- do.call(rbind, lapply(files, readRDS))
output_file <- "comparePredsRevisionResults.rds"
saveRDS(results, output_file)

### Also for two additional situations for supplement
files_dir <- "cluster_output_compare_preds_pairing_me"
files <- list.files(files_dir, full.names = TRUE)
results <- do.call(rbind, lapply(files, readRDS))
output_file <- "comparePredsPairingMEResults.rds"
saveRDS(results, output_file)

### Simulations for estimation of just M
files_dir <- "cluster_output_sims_no_covar"
files <- list.files(files_dir, full.names = TRUE)
results <- lapply(files, readRDS)
output_file <- "noCovarResults.rds"
saveRDS(results, output_file)

### Equal rows test
files_dir <- "cluster_output_equal_rows_B_test"
files <- list.files(files_dir, full.names = TRUE)
results <- do.call(rbind, lapply(files, readRDS))
output_file <- "equalRowBTestResults.rds"
saveRDS(results, output_file)
