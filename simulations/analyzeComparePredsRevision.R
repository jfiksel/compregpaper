library(tidyverse)
library(here)
library(gtable)
library(grid)
library(ggpubr)
results_df <- readRDS(here("simulations", "comparePredsRevisionResults.rds"))
avg_kld_df <-
    results_df %>%
    group_by(model, N, param_index, true_model, datagen) %>%
    summarize(avg_kld = mean(kld), max_kld = max(kld))
avg_kld_df$param_index <- factor(avg_kld_df$param_index,
                                 levels = c(1,2),
                                 labels = c("Education", "White Cells"))
avg_kld_df$model <- factor(avg_kld_df$model,
                           levels = c("compreg", "ilr", "alr", "mflr"),
                           labels = c("Direct Regression",
                                      "Log-ratio",
                                      "Log-ratio (ALR)",
                                      "Pseudo-ML"))
avg_kld_df$datagen <- factor(avg_kld_df$datagen,
                             levels = c('dirichlet', "multinomial", 'dirichlet-multinomial', 'ilr-normal', 'alr-normal'),
                             labels = c('Dirichlet', "Multinomial", 'Dirichlet-Multinomial', 'ilr-normal', 'Logistic-Normal'))
avg_kld_df$true_model <- factor(avg_kld_df$true_model,
                             levels = c("compreg", "ILR", "ALR", "mflr"),
                             labels = c("True model \n Direct Regression",
                                        "True model \n Log-ratio",
                                        "True model \n Log-ratio (ALR)",
                                        "True model \n Pseudo-ML"))

### First plot: Different true models
fig1a <-
    avg_kld_df %>%
    filter((true_model == "True model \n Direct Regression" & datagen == "Dirichlet") | true_model != "True model \n Direct Regression") %>%
    mutate(true_model = as.character(true_model)) %>%
    mutate(true_model = ifelse(true_model == "True model \n Direct Regression",
                               "True model \n Direct Regression \n (Dirichlet)",
                               true_model)) %>%
    filter(!grepl('ALR', model), !grepl('ALR', true_model)) %>%
    ggplot(aes(x = N, y = log(avg_kld), color = model, shape = model, linetype = factor(param_index))) +
    geom_point(alpha = .8, size = 1.75) +
    geom_line() +
    facet_grid(~true_model, scales = "free_y") +
    xlab("N") +
    ylab("Log KLD") +
    scale_color_discrete(name = "Fitted model") +
    scale_shape_discrete(name = "Fitted model") +
    scale_linetype_discrete(name = "Dataset")
    

fig1b <-
    avg_kld_df %>%
    filter(true_model == "True model \n Direct Regression") %>%
    filter(!grepl('ALR', model), datagen != 'ilr-normal', datagen != 'Dirichlet') %>%
    mutate(datagen = paste0("True model \n Direct Regression \n (", datagen, ")")) %>%
    mutate(datagen = factor(datagen,
                            levels = paste0("True model \n Direct Regression \n (",
                                            c("Multinomial",
                                            "Dirichlet-Multinomial",
                                            "Logistic-Normal"),
                                            ")"))) %>%
    ggplot(aes(x = N, y = log(avg_kld), color = model, shape = model, linetype = factor(param_index))) +
    geom_point(alpha = .8, size = 1.75) +
    geom_line() +
    facet_grid(~datagen, scales = "free_y") +
    xlab("N") +
    ylab("Log KLD") +
    scale_color_discrete(name = "Fitted model") +
    scale_shape_discrete(name = "Fitted model") +
    scale_linetype_discrete(name = "Dataset")

fig1 <- ggarrange(fig1a, fig1b,
          labels=c('A', 'B'),
          common.legend = T, nrow = 2,
          legend = 'right')
ggsave(here('simulations', 'figs','fig1resub.pdf'), fig1,
       width = 8, height = 6)


### Now for additional scenarios of paired measurements and measurement error
results_df <- readRDS(here("simulations", "comparePredsPairingMEResults.rds"))
avg_kld_df <-
    results_df %>%
    group_by(model, N, param_index, datagen) %>%
    summarize(avg_kld = mean(kld), max_kld = max(kld)) %>%
    filter(model != "alr")
avg_kld_df$model <- factor(avg_kld_df$model,
                           levels = c("compreg", "mflr", "ilr"),
                           labels = c("Direct Regression",
                                      "Pseudo-ML",
                                      "Log-ratio"))
avg_kld_df$datagen <- factor(avg_kld_df$datagen,
                             levels = c('paired_aggregation', 'measurement_error'),
                             labels = c('Paired Aggregation', 'Measurement Error'))
avg_kld_df$param_index <- factor(avg_kld_df$param_index,
                                 levels = 1:2,
                                 labels = c("B1", "B2"))

pairing_me_plot <-
    ggplot(avg_kld_df, aes(x = N, y = log(avg_kld), color = model,shape = model, linetype = param_index)) +
    geom_point(alpha = .8, size = 1.75) +
    geom_line() +
    facet_grid(~datagen, scales = "free_y") +
    xlab("N") +
    ylab("Log KLD") +
    scale_color_discrete(name = "Fitted model") +
    scale_shape_discrete(name = "Fitted model") +
    scale_linetype_discrete(name = "Value of B")
ggsave(here('simulations', 'figs','suppfig1resub.pdf'), pairing_me_plot,
       width = 8, height = 4)
