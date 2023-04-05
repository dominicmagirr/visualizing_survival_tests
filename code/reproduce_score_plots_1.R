here::i_am("code/reproduce_score_plots_1.R")
library(here)

library(dplyr)
library(ggplot2)
library(survival)
library(pch)
library(flexsurv)
library(nphRCT)


source(here("utils", "parametric_milestone.R"))
source(here("utils", "window_mean_survival.R"))

## get data
df <- readxl::read_excel(here("data", "41591_2018_134_MOESM3_ESM.xlsx"),
                         sheet = 2) %>% 
  select(PtID, ECOGGR, OS, OS.CNSR, TRT01P) %>%
  mutate(event = -1 * (OS.CNSR - 1),
         time = OS,
         arm = ifelse(TRT01P == "Docetaxel", "0", "1")) %>% 
  select(time, event, arm)

## ------------------------------------------------------------------------------

## Figure 2

## get scores
df_lr <- find_scores(formula=Surv(time, event) ~ arm, data=df, method = "lr")
df_milestone_pseudo <- find_scores(formula=Surv(time, event) ~ arm, tau=18, data=df, method = "ms")
df_mwlrt <- find_scores(formula=Surv(time, event) ~ arm,  data=df, s_star=0.5, method = "mw")

df_surv_pseudo_p <- get_scores_surv_p(df, tau = 18, breaks = c(0,100))
df_piece <- get_scores_surv_p(df, tau = 18, breaks = c(0,2,4,6,8,100))
df_flex <- get_scores_flexsurv_p(df, 18, k = 2)


# tidy up 
df_lr <- df_lr$df %>% mutate(label = "Log-rank", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_mwlrt <- df_mwlrt$df %>% mutate(label = "Modestly-weighted log-rank (s* = 0.5)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_surv_pseudo <- df_milestone_pseudo$df %>% mutate(label = "KM milestone (tau = 18)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))

df_surv_pseudo_p <- df_surv_pseudo_p$df %>% mutate(label = "Param. milestone (tau = 18)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_piece <- df_piece$df %>% mutate(label = "Piecewise param. milestone (tau = 18)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_flex <- df_flex$df %>% mutate(label = "Flexible param. milestone (tau = 18)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))

# stack data frames
modest_comp = rbind(df_lr,
                    df_mwlrt, 
                    df_surv_pseudo,
                    df_surv_pseudo_p,
                    df_piece,
                    df_flex) %>%
  mutate(type = ifelse(event == 1, "Event", "Censored"))

# re-order 
modest_comp$label = factor(modest_comp$label,
                           levels = unique(modest_comp$label))

# plot
p_modest_comp = ggplot() +
  geom_point(data = modest_comp,  aes(x = time, y = standardized_score, color = group, alpha = type)) + 
  scale_alpha_discrete(range = c(0.2, 1)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  ylab("Score") +
  labs(color = "Arm", alpha = "Observation type") +
  geom_hline(data = modest_comp %>% group_by(group, label) %>% dplyr::summarize(mean_score = mean(standardized_score)), 
             aes(yintercept = mean_score, colour = group), linetype = 2) +
  facet_wrap(~ label)


#save
ggsave(here("output", "p_modest_comp.pdf"),
       p_modest_comp, width = 8, height = 5, dpi = 300) 


## ------------------------------------------------------------------------------

## Figure 3

# get scores
df_rmst <- find_scores(formula=Surv(time, event) ~ arm, tau=24, data=df, method = "rmst")
df_fh_1_0 <- find_scores(formula=Surv(time, event) ~ arm,  data=df, rho=1, gamma = 0, method = "fh")

# tidy up
df_rmst <- df_rmst$df %>% mutate(label = "KM RMST(tau = 24)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_fh_1_0 <- df_fh_1_0$df %>% mutate(label = "FH(1,0)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))

# stack data frames
rmst_comp = rbind(df_lr,
                  df_rmst, 
                  df_fh_1_0) %>%
  mutate(type = ifelse(event == 1, "Event", "Censored"))


# re-order
rmst_comp$label = factor(rmst_comp$label,
                         levels = unique(rmst_comp$label))

#plot
p_rmst_comp = ggplot() +
  geom_point(data = rmst_comp,  aes(x = time, y = standardized_score, color = group, alpha = type)) + 
  scale_alpha_discrete(range = c(0.2, 1)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  ylab("Score") +
  labs(color = "Arm", alpha = "Observation type") +
  geom_hline(data = rmst_comp %>% group_by(group, label) %>% dplyr::summarize(mean_score = mean(standardized_score)), 
             aes(yintercept = mean_score, colour = group), linetype = 2) +
  facet_wrap(~ label)



ggsave(here("output", "p_rmst_comp.pdf"), p_rmst_comp, width = 8, height = 3, dpi = 300) 



## ------------------------------------------------------------------------------

## Figure 4

# get scores
df_wmst <- find_pseudovalues_wmst(formula=Surv(time, event) ~ arm,  data=df, tau_1 = 12, tau_2 = 24)
df_ahsw <- find_pseudovalues_ahsw(formula=Surv(time, event) ~ arm,  data=df, tau = 18)


# tidy up
df_wmst <- df_wmst$df %>% mutate(label = "WMST(12,24)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_ahsw <- df_ahsw$df %>% mutate(label = "AHSW(tau = 18)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))

# stack
novel_comp = rbind(df_wmst,
                   df_ahsw) %>%
  mutate(type = ifelse(event == 1, "Event", "Censored"))

# re-order
novel_comp$label = factor(novel_comp$label,
                          levels = unique(novel_comp$label))

# plot
p_novel_comp = ggplot() +
  geom_point(data = novel_comp,  aes(x = time, y = standardized_score, color = group, alpha = type)) + 
  scale_alpha_discrete(range = c(0.2, 1)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  ylab("Score") +
  labs(color = "Arm", alpha = "Observation type") +
  geom_hline(data = novel_comp %>% group_by(group, label) %>% dplyr::summarize(mean_score = mean(standardized_score)), 
             aes(yintercept = mean_score, colour = group), linetype = 2) +
  facet_wrap(~ label)



ggsave(here("output", "p_novel_comp.pdf"), p_novel_comp, width = 6, height = 3, dpi = 300) 


## ------------------------------------------------------------------------------

## Figure 5

# get scores
df_fh_0_1 <- find_scores(formula=Surv(time, event) ~ arm,  data=df, rho=0, gamma = 1, method = "fh")
df_fh_1_1 <- find_scores(formula=Surv(time, event) ~ arm,  data=df, rho=1, gamma = 1, method = "fh")
df_fh_0_half <- find_scores(formula=Surv(time, event) ~ arm,  data=df, rho=0, gamma = 0.5, method = "fh")

# tidy up
df_fh_0_1 <- df_fh_0_1$df %>% mutate(label = "FH(0,1)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_fh_1_1 <- df_fh_1_1$df %>% mutate(label = "FH(1,1)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))
df_fh_0_half <- df_fh_0_half$df %>% mutate(label = "FH(0,0.5)", time = t_j) %>% dplyr::select(c(time, event, group, standardized_score, label))

# stack
fh_comp = rbind(df_fh_0_1,
                df_fh_1_1, 
                df_fh_0_half) %>%
  mutate(type = ifelse(event == 1, "Event", "Censored"))

# re-order
fh_comp$label = factor(fh_comp$label,
                       levels = unique(fh_comp$label))

# plot
p_fh_comp = ggplot() +
  geom_point(data = fh_comp,  aes(x = time, y = standardized_score, color = group, alpha = type)) + 
  scale_alpha_discrete(range = c(0.2, 1)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  ylab("Score") +
  labs(color = "Arm", alpha = "Observation type") +
  geom_hline(data = fh_comp %>% group_by(group, label) %>% dplyr::summarize(mean_score = mean(standardized_score)), 
             aes(yintercept = mean_score, colour = group), linetype = 2) +
  facet_wrap(~ label)



ggsave(here("output", "p_fh_comp.pdf"), p_fh_comp, width = 8, height = 3, dpi = 300) 


