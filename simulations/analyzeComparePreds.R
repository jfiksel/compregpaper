library(tidyverse)
library(here)
library(gtable)
library(grid)
results_df <- readRDS(here("simulations", "comparePredsResults.rds"))
avg_kld_df <-
    results_df %>%
    group_by(model, N, param_index, true_model, datagen) %>%
    summarize(avg_kld = mean(kld), max_kld = max(kld))
avg_kld_df$param_index <- factor(avg_kld_df$param_index,
                                 levels = c(1,2),
                                 labels = c("Education", "White Cells"))
avg_kld_df$model <- factor(avg_kld_df$model,
                           levels = c("compreg", "mflr", "ilr"),
                           labels = c("Direct Regression", "Alenazi (2019)", "Chen et al. (2017)"))
avg_kld_df$datagen <- factor(avg_kld_df$datagen,
                             levels = c('dirichlet', "multinomial", 'dirichlet-multinomial', 'ilr-normal'),
                             labels = c('Dirichlet', "Multinomial", 'Dirichlet-Multinomial', 'ilr-normal'))

### First plot showing when true mean is CompReg
fig_true_compreg <-
    avg_kld_df %>%
    filter(true_model == "compreg", datagen != "ilr-normal", model == "Direct Regression") %>%
    ggplot(aes(x = N, y = log(avg_kld), linetype = datagen)) +
    geom_point(alpha = .8) +
    geom_line() +
    facet_grid(~param_index, scales = "free_y") +
   # theme(legend.position="bottom") +
    xlab("N") +
    ylab("Log KLD") +
    scale_linetype_manual(name = "Data generating mechanism",
                          values = c("Dirichlet" = "solid",
                                     "Multinomial" = "dotted",
                                     "Dirichlet-Multinomial" = "dashed"))
ggsave(here("simulations", "figs", "fig2.pdf"), fig_true_compreg,
       width = 8, height = 4)

### compare to other models
fig_true_compreg_supp <-
    avg_kld_df %>%
    filter(true_model == "compreg", datagen != "ilr-normal") %>%
    ggplot(aes(x = N, y = log(avg_kld), linetype = datagen, color = model)) +
    geom_point(alpha = .8) +
    geom_line() +
    facet_grid(~param_index, scales = "free_y") +
    # theme(legend.position="bottom") +
    xlab("N") +
    ylab("Log KLD") +
    scale_color_discrete(name = "Fitted Model") +
    scale_linetype_manual(name = "Data generating mechanism",
                          values = c("Dirichlet" = "solid",
                                     "Multinomial" = "dotted",
                                     "Dirichlet-Multinomial" = "dashed"))
ggsave(here("simulations", "figs", "suppfig1.pdf"), fig_true_compreg_supp,
       width = 8, height = 4)

### Now plot showing across different true models

### Data frame where true mean is CompReg, datagen is dirichlet

compreg_df <- filter(avg_kld_df, true_model == "compreg", datagen == "Dirichlet")

### Data frame where true mean is MFLR, datagen is dirichlet

mflr_df <- filter(avg_kld_df, true_model == "mflr", datagen == "Dirichlet")

### Data frame where true mean is ILR, datagen is ILR-Normal

ilr_df <- filter(avg_kld_df, true_model == "ILR", datagen == "ilr-normal")

plot_df <- bind_rows(compreg_df, mflr_df) %>% bind_rows(ilr_df)
plot_df$true_model <- factor(plot_df$true_model,
                           levels = c("compreg", "mflr", "ILR"),
                           labels = c("True model \n Direct Regression",
                                      "True model \n Alenazi (2019)",
                                      "True model \n Chen et al. (2017)"))

compare_true_models_plot <-
    plot_df %>%
    ggplot(aes(x = N, y = log(avg_kld), color = model, linetype = factor(param_index))) +
    geom_point(alpha = .8) +
    geom_line() +
    facet_grid(~true_model, scales = "free_y") +
    xlab("N") +
    ylab("Log KLD") +
    scale_color_discrete(name = "Fitted model") +
    scale_linetype_discrete(name = "Dataset")
ggsave(here("simulations", "figs", "fig1.pdf"), compare_true_models_plot,
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