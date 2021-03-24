library(tidyverse)
library(here)
library(gtable)
library(grid)
library(cowplot)
results_df <- readRDS(here("simulations", "comparePredsALRResults.rds"))
avg_kld_df <-
    results_df %>%
    group_by(model, N, param_index, true_model, datagen) %>%
    summarize(avg_kld = mean(kld), max_kld = max(kld))
avg_kld_df$param_index <- factor(avg_kld_df$param_index,
                                 levels = c(1,2),
                                 labels = c("Education", "White Cells"))
avg_kld_df$model <- factor(avg_kld_df$model,
                           levels = c("compreg", "alr"),
                           labels = c("Direct Regression", "ALR"))
avg_kld_df$true_model <- factor(avg_kld_df$true_model,
                           levels = c("compreg", "ALR"),
                           labels = c("Direct Regression", "ALR"))
avg_kld_df$datagen <- factor(avg_kld_df$datagen,
                             levels = c('dirichlet', "multinomial", 'dirichlet-multinomial', 'alr-normal'),
                             labels = c('Dirichlet', "Multinomial", 'Dirichlet-Multinomial', 'Logistic-Normal'))


compare_true_models_plot <-
    avg_kld_df %>%
    filter(true_model == "Direct Regression") %>%
    ggplot(aes(x = N, y = log(avg_kld), color = model, linetype = factor(param_index))) +
    geom_point(alpha = .8) +
    geom_line() +
    facet_grid(~datagen, scales = "free_y") +
    xlab("N") +
    ylab("Log KLD") +
    scale_color_discrete(name = "Fitted model") +
    scale_linetype_discrete(name = "Dataset") +
    theme(legend.position = "none")
compare_true_models_plot_alr <-
    avg_kld_df %>%
    filter(true_model == "ALR", datagen == "Logistic-Normal") %>%
    ggplot(aes(x = N, y = log(avg_kld), color = model, linetype = factor(param_index))) +
    geom_point(alpha = .8) +
    geom_line() +
    xlab("N") +
    ylab("Log KLD") +
    scale_color_discrete(name = "Fitted model") +
    scale_linetype_discrete(name = "Dataset")
kld_plot_alr <- cowplot::plot_grid(compare_true_models_plot, compare_true_models_plot_alr,
                   labels = "AUTO", ncol = 1, rel_widths = c(2,1))

ggsave(here("simulations", "figs", "modelComparisonALR.pdf"), kld_plot_alr ,
       width = 8, height = 6)

ggsave(here("simulations", "figs", "modelComparisonALR1.pdf"), compare_true_models_plot,
       width = 8, height = 4)

z <- ggplotGrob(p)

#  New strip at the top
z <- gtable_add_rows(z, z$height[7], pos = 6)  # New row added below row 6

# Check the layout
gtable_show_layout(z)   # New strip goes into row 7 
# New strip spans columns 5 to 9

z <- gtable_add_grob(z, 
                     list(rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
                          textGrob("True model", gp = gpar(cex = .75, fontface = 'bold', col = "black"))), 
                     t=7, l=5, b=7, r=9, name = c("a", "b"))

# Add small gap between strips - below row 6
z <- gtable_add_rows(z, unit(2/10, "line"), 7)

# Draw it
grid.newpage()
grid.draw(z)


index_1_supp_plot <-
    avg_kld_df %>%
    filter(param_index == 1) %>%
    ggplot(aes(x = N, y = log(avg_kld), color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(true_model ~ datagen, scales = "free_y")

index_2_plot <-
    plot_df %>%
    filter(param_index == 2) %>%
    ggplot(aes(x = N, y = log(avg_kld), color = model)) +
    geom_point() +
    geom_line() +
    facet_wrap(~true_model, scales = "free_y")

index_2_supp_plot <-
    avg_kld_df %>%
    filter(param_index == 2) %>%
    ggplot(aes(x = N, y = log(avg_kld), color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(true_model ~ datagen, scales = "free_y")


dir_compreg_plot <-
    avg_kld_df %>%
    filter(dist == "dirichlet", true_model == "compreg") %>%
    ggplot(aes(x = N, y = avg_kld, color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(Mind ~ phi) +
    ggtitle("True model dirichlet, compreg")

dir_mflr_plot <-
    avg_kld_df %>%
    filter(dist == "dirichlet", true_model == "mflr") %>%
    ggplot(aes(x = N, y = avg_kld, color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(Mind ~ phi) +
    ggtitle("True model dirichlet, mflr")

mult_compreg_plot <-
    avg_kld_df %>%
    filter(dist == "multinomial", true_model == "compreg") %>%
    ggplot(aes(x = N, y = avg_kld, color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(Mind ~ nmax) +
    ggtitle("True model multinomial, compreg")

mult_mflr_plot <-
    avg_kld_df %>%
    filter(dist == "multinomial", true_model == "mflr") %>%
    ggplot(aes(x = N, y = avg_kld, color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(Mind ~ nmax) +
    ggtitle("True model multinomial, mflr")

mult_dir_compreg_plot <-
    avg_kld_df %>%
    filter(dist == "multinomial-dirichlet", true_model == "compreg") %>%
    mutate(nmax_phi = paste0("nmax ", nmax, "phi ", phi)) %>%
    ggplot(aes(x = N, y = avg_kld, color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(Mind ~ nmax_phi) +
    ggtitle("True model multinomial-dirichlet, compreg")

mult_dir_mflr_plot <-
    avg_kld_df %>%
    filter(dist == "multinomial-dirichlet", true_model == "mflr") %>%
    mutate(nmax_phi = paste0("nmax ", nmax, "phi ", phi)) %>%
    ggplot(aes(x = N, y = avg_kld, color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(Mind ~ nmax_phi) +
    ggtitle("True model multinomial-dirichlet, mflr")

ilr_plot <-
    avg_kld_df %>%
    filter(dist == "ilr_regression") %>%
    ggplot(aes(x = N, y = avg_kld, color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(Mind ~ sigma2) +
    ggtitle("True model ILR")
pdf(here("simulations", "figs", "comparePredsFirstPlots.pdf"), width = 10, height = 10)
print(dir_compreg_plot)
print(dir_mflr_plot)
print(mult_compreg_plot)
print(mult_mflr_plot)
print(mult_dir_compreg_plot)
print(mult_dir_mflr_plot)
print(ilr_plot)
dev.off()