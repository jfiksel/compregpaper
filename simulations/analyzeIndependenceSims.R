library(tidyverse)
library(here)
library(xtable)
pval_df <- readRDS(here("simulations", "pvalsIndependenceLLR.rds"))
pval_summary <-
    pval_df %>%
    group_by(datagen, Mindex, n, method, true_model) %>%
    summarise(rejectProb = mean(p <= .05),
              sumReject = sum(p <= .05))

### Tables for paper
### Table 1
### Independence settings
pval_summary %>%
    filter(method == "compreg", Mindex == 4, true_model == 'compreg')
pval_summary %>%
    filter(method == "compreg", Mindex == 4, true_model == 'ILR')
pval_summary %>%
    filter(method == "ILR", Mindex == 4)


pval_summary  %>%
    filter(Mindex == 1, method == 'compreg') %>%
    mutate(typeIIerror = 1 - rejectProb)

pval_summary  %>%
    filter(Mindex == 1, method == 'ILR') %>%
    mutate(typeIIerror = 1 - rejectProb)

pval_summary  %>%
    filter(Mindex == 2, method == 'compreg') %>%
    mutate(typeIIerror = 1 - rejectProb)

pval_summary  %>%
    filter(Mindex == 2, method == 'ILR') %>%
    mutate(typeIIerror = 1 - rejectProb)

pval_summary  %>%
    filter(Mindex == 3, method == 'compreg') %>%
    mutate(typeIIerror = 1 - rejectProb)

pval_summary  %>%
    filter(Mindex == 3, method == 'ILR') %>%
    mutate(typeIIerror = 1 - rejectProb)


###################################
### Equal rows b test
equal_rows_df <- readRDS(here("simulations", "equalRowBTestResults.rds"))
rejection_rates <-
    equal_rows_df %>%
    group_by(datagen, pairindex) %>%
    summarise(reject = mean(p < .05))

# pval_chisq_df <- readRDS(here("simulations", "pvalsIndependence.rds"))
# 
# pval_df <- bind_rows(mutate(pval_llr_df, method = "LLR"),
#                      mutate(pval_chisq_df, method = "ChiSq"))
# 

# ggplot(pval_summary, aes(x = datagen, y = rejectProb, color = method)) +
#     geom_point() +
#     facet_grid(n ~ Mindex)
# indep_setting <- pval_summary %>% filter(Mindex == 4)

### Tables for paper


    


