output_dir <- "cluster_output_independence_test"
files <- list.files(output_dir, full.names = TRUE)
pval_df <- do.call(rbind, lapply(files, readRDS))
saveRDS(pval_df, "pvalsIndependence.rds")

output_dir <- "cluster_output_independence_test_loss_function"
files <- list.files(output_dir, full.names = TRUE)
pval_df <- do.call(rbind, lapply(files, readRDS))
saveRDS(pval_df, "pvalsIndependenceLLR.rds")


output_dir <- "cluster_output_independence_test_multinomial"
files <- list.files(output_dir, full.names = TRUE)
pval_df <- do.call(rbind, lapply(files, readRDS))
saveRDS(pval_df, "pvalsIndependenceMultinomial.rds")