here::i_am("code/reproduce_score_plots_2.R")
library(here)

library(dplyr)
library(ggplot2)
library(survival)
library(pch)
library(flexsurv)
library(nphRCT)


source(here("utils", "parametric_milestone.R"))

## get data
source(here("code", "recensoring.R"))


## get scores
df_lr2 <- find_scores(formula=Surv(time, event) ~ arm, data=df2, method = "lr")
df_surv_pseudo_p2 <- get_scores_surv_p(df2, tau = 18, breaks = c(0,100))
df_rmst2 <- find_scores(formula=Surv(time, event) ~ arm, tau=24, data=df2, method = "rmst")
df_mwlrt2 <- find_scores(formula=Surv(time, event) ~ arm,  data=df2, s_star=0.5,method = "mw")
df_milestone_pseudo2 <- find_scores(formula=Surv(time, event) ~ arm, tau=18,data=df2,method = "ms")
df_flex2 <- get_scores_flexsurv_p(df2, 18, k = 2)
df_wmst2 <- find_pseudovalues_wmst(formula=Surv(time, event) ~ arm,  data=df2, tau_1 = 12, tau_2 = 24)
df_ahsw2 <- find_pseudovalues_ahsw(formula=Surv(time, event) ~ arm,  data=df2, tau = 18)

# tidy up
df_lr2 <- df_lr2$df %>% mutate(label = "Log-rank", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_surv_pseudo_p2 <- df_surv_pseudo_p2$df %>% mutate(label = "Param. milestone (tau = 18)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_rmst2 <- df_rmst2$df %>% mutate(label = "KM RMST(tau = 24)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_mwlrt2 <- df_mwlrt2$df %>% mutate(label = "Modestly-weighted log-rank (s* = 0.5)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_surv_pseudo2 <- df_milestone_pseudo2$df %>% mutate(label = "K-M milestone (tau = 18)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_flex2 <- df_flex2$df %>% mutate(label = "Flexible param. milestone (tau = 18)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_wmst2 <- df_wmst2$df %>% mutate(label = "WMST(12,24)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_ahsw2 <- df_ahsw2$df %>% mutate(label = "AHSW(tau = 18)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))

# stack
cens_comp = rbind(df_lr2,
                  df_mwlrt2,
                  df_rmst2, 
                  df_surv_pseudo_p2,
                  df_surv_pseudo2,
                  df_flex2,
                  df_wmst2,
                  df_ahsw2) %>%
  mutate(type = ifelse(event == 1, "Event", "Censored"))

# re-order
cens_comp$label = factor(cens_comp$label,
                         levels = unique(cens_comp$label))

# plot
p_cens_comp = ggplot() +
  geom_point(data = cens_comp,  aes(x = time, y = standardized_score, color = group, alpha = type)) + 
  scale_alpha_discrete(range = c(0.2, 1)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  ylab("Score") +
  labs(color = "Arm", alpha = "Observation type") +
  geom_hline(data = cens_comp %>% group_by(group, label) %>% dplyr::summarize(mean_score = mean(standardized_score)), 
             aes(yintercept = mean_score, colour = group), linetype = 2) +
  facet_wrap(~ label)



ggsave(here("output", "p_cens_comp.pdf"), p_cens_comp, width = 8, height = 7, dpi = 300) 




