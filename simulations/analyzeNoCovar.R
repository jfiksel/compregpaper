library(here)
library(tidyverse)
library(ggpubr)
results_list <- readRDS(here("simulations", "noCovarResults.rds"))
### three different transition matrices
C <- 4
M1 <- diag(C)
M2 <- matrix(c(1, 0, 0, 0,
               .65, .35, 0, 0,
               0, 0, 0.5, 0.5,
               0, 0, 0, 1),
             nrow = C, byrow = TRUE)
M3 <- .6*diag(C) + .1
M_list <- list(M1, M2, M3)
### Function for KL distance
kld <- function(true, est) {
    d <- 0
    for(i in 1:length(true)) {
        if(true[i] == 0 | est[i] == 0) {
            d <- d 
        } else {
            d <- d + true[i] * log(true[i] / est[i])
        }
    }
    return(d)
}

####### KLD Distance results and plots
if(!file.exists(here("simulations", "noCovarKLD.rds"))) {
    results_df <- do.call(rbind, lapply(1:length(results_list), function(i) {
        sims <- results_list[[i]]
        setting <- sims$setting
        Mtrue <- M_list[[setting$Mind]]
        m_est_list <- sims$m_est_list
        kld_df_all <- do.call(rbind, lapply(1:length(m_est_list), function(j) {
            Mest <- m_est_list[[j]]
            kld_vec <- sapply(1:C, function(c) { 
                kld(Mtrue[c,], Mest[c,])
            })
            kld_df <- data.frame(kld = kld_vec, row = 1:C, sim = j)
            return(kld_df)
        }))
        kld_df_all <- cbind(kld_df_all, setting)
        return(kld_df_all)
    }))
    saveRDS(results_df, here("simulations", "noCovarKLD.rds"))
} else {
    results_df <- readRDS(here("simulations", "noCovarKLD.rds"))
}

kld_avg <-
    results_df %>%
    group_by(row, dist, phi, nmax, N, Mind) %>%
    summarize(kld_avg = mean(kld)) %>%
    mutate(Mind = paste0("M", Mind))

### Now plot results for N=250 and for different data generating scenarios
dirichlet_plot_N_250 <-
    kld_avg %>%
    filter(dist == "dirichlet", N == 250) %>%
    ggplot(aes(x = phi, y = kld_avg, color = factor(row))) +
    geom_line() +
    geom_point() +
    facet_wrap(~Mind, nrow = 1) +
    ggtitle("Dirichlet distributed outcome") +
    labs(color = element_text("Row of M"))

multinomial_plot_N_250 <-
    kld_avg %>%
    filter(dist == "multinomial", N == 250) %>%
    ggplot(aes(x = nmax, y = kld_avg, color = factor(row))) +
    geom_line() +
    geom_point() +
    facet_wrap(~Mind, nrow = 1) +
    ggtitle("Multinomial distributed outcome") +
    labs(color = element_text("Row of M"))

multinomial_dirichlet_plot_fixed_phi_N_250 <-
    kld_avg %>%
    filter(dist == "multinomial-dirichlet", N == 250, phi == 10) %>%
    ggplot(aes(x = nmax, y = kld_avg, color = factor(row))) +
    geom_line() +
    geom_point() +
    facet_wrap(~Mind, nrow = 1) +
    ggtitle("Dirichlet-Multinomial distributed outcome, phi = 10") +
    labs(color = element_text("Row of M"))

multinomial_dirichlet_plot_fixed_nmax_N_250 <-
    kld_avg %>%
    filter(dist == "multinomial-dirichlet", N == 250, nmax == 11) %>%
    ggplot(aes(x = phi, y = kld_avg, color = factor(row))) +
    geom_line() +
    geom_point() +
    facet_wrap(~Mind, nrow = 1) +
    ggtitle("Dirichlet-Multinomial distributed outcome, nmax = 11") +
    labs(color = element_text("Row of M"))

plots_combined_N_250 <-
    ggarrange(dirichlet_plot_N_250,
              multinomial_plot_N_250,
              multinomial_dirichlet_plot_fixed_phi_N_250,
              multinomial_dirichlet_plot_fixed_nmax_N_250,
              common.legend = TRUE,
              legend = "bottom")
figs_dir <- here('simulations', 'figs')
if(!dir.exists(figs_dir)) {
    dir.create(figs_dir, recursive = TRUE)
}
ggsave(file.path(figs_dir, "kld_N250_no_covar.pdf"),
       plots_combined_N_250, width = 10, height = 10)

### Now plot results over N=250 and for different data generating scenarios
dirichlet_plot_phi_10 <-
    kld_avg %>%
    filter(dist == "dirichlet", phi == 10) %>%
    ggplot(aes(x = N, y = kld_avg, color = factor(row))) +
    geom_line() +
    geom_point() +
    facet_wrap(~Mind, nrow = 1) +
    ggtitle("Dirichlet distributed outcome phi = 10") +
    labs(color = element_text("Row of M"))

multinomial_plot_nmax_11 <-
    kld_avg %>%
    filter(dist == "multinomial", nmax == 11) %>%
    ggplot(aes(x = N, y = kld_avg, color = factor(row))) +
    geom_line() +
    geom_point() +
    facet_wrap(~Mind, nrow = 1) +
    ggtitle("Multinomial distributed outcome nmax = 11") +
    labs(color = element_text("Row of M"))

multinomial_dirichlet_plot_fixed_phi_fixed_nmax <-
    kld_avg %>%
    filter(dist == "multinomial-dirichlet", nmax == 11, phi == 10) %>%
    ggplot(aes(x = N, y = kld_avg, color = factor(row))) +
    geom_line() +
    geom_point() +
    facet_wrap(~Mind, nrow = 1) +
    ggtitle("Dirichlet-Multinomial distributed outcome, phi = 10, nmax = 11") +
    labs(color = element_text("Row of M"))

plots_combined_all_N <-
    ggarrange(dirichlet_plot_phi_10,
              multinomial_plot_nmax_11,
              multinomial_dirichlet_plot_fixed_phi_fixed_nmax,
              common.legend = TRUE,
              legend = "bottom")
figs_dir <- here('simulations', 'figs')
if(!dir.exists(figs_dir)) {
    dir.create(figs_dir, recursive = TRUE)
}
ggsave(file.path(figs_dir, "kld_all_N_no_covar.pdf"),
       plots_combined_all_N, width = 10, height = 10)


##############################################
### Data frame with specific values of M versus truth
if(!file.exists(here("simulations", "noCovarMdf.rds"))) {
    results_df <- do.call(rbind, lapply(1:length(results_list), function(i) {
        sims <- results_list[[i]]
        setting <- sims$setting
        Mtrue <- M_list[[setting$Mind]]
        Mtruevec <- as.vector(Mtrue)
        Mnames <- paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), "]")
        m_est_list <- sims$m_est_list
        M_df_all <- do.call(rbind, lapply(1:length(m_est_list), function(j) {
            Mestvec <- as.vector(m_est_list[[j]])
            M_df <- data.frame(mtrue = Mtruevec,
                               mest = Mestvec,
                               mname = Mnames)
            return(M_df)
        }))
        M_df_all <- cbind(M_df_all, setting)
        return(M_df_all)
    }))
    saveRDS(results_df, here("simulations", "noCovarMdf.rds"))
} else {
    results_df <- readRDS(here("simulations", "noCovarMdf.rds"))
}

#### Look at estimates of M for dirichlet-multinomial phi = 10
M1_dir_mult_plot <-
    results_df %>%
    filter(Mind == 1, phi == 10, dist == "multinomial-dirichlet") %>%
    ggplot(aes(x = mest, color = factor(nmax))) +
    geom_line(stat = 'density') +
    geom_vline(aes(xintercept = mtrue)) +
    facet_wrap(~mname, scales = "free_y") +
    labs(color = "nmax") +
    ggtitle("Estimates of M1, Dirichlet-multinomial phi = 10, N = 250")
M2_dir_mult_plot <-
    results_df %>%
    filter(Mind == 2, phi == 10, dist == "multinomial-dirichlet") %>%
    ggplot(aes(x = mest, color = factor(nmax))) +
    geom_line(stat = 'density') +
    geom_vline(aes(xintercept = mtrue)) +
    facet_wrap(~mname, scales = "free_y") +
    labs(color = "nmax") +
    ggtitle("Estimates of M2, Dirichlet-multinomial phi = 10, N = 250")
M3_dir_mult_plot <-
    results_df %>%
    filter(Mind == 3, phi == 10, dist == "multinomial-dirichlet") %>%
    ggplot(aes(x = mest, color = factor(nmax))) +
    geom_line(stat = 'density') +
    geom_vline(aes(xintercept = mtrue)) +
    facet_wrap(~mname) +
    labs(color = "nmax") +
    ggtitle("Estimates of M3, Dirichlet-multinomial phi = 10, N = 250")
pdf(file.path(figs_dir, "M_estimates_dirichlet_multinomial.pdf"), width = 8, height = 8)
print(M1_dir_mult_plot)
print(M2_dir_mult_plot)
print(M3_dir_mult_plot)
dev.off()
