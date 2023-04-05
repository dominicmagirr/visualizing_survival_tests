here::i_am("code/recensoring.R")
library(here)

library(dplyr)
library(ggplot2)
library(survival)

df <- readxl::read_excel(here("data", "41591_2018_134_MOESM3_ESM.xlsx"),
                         sheet = 2) %>% 
  select(PtID, ECOGGR, OS, OS.CNSR, TRT01P) %>%
  mutate(event = -1 * (OS.CNSR - 1),
         time = OS,
         arm = ifelse(TRT01P == "Docetaxel", "0", "1")) %>% 
  select(time, event, arm)


df2 =df
head(df2)
length(df2$time)
set.seed(314)

df2$add_cens <- runif(287) * 26
df2$event <- ifelse(df2$add_cens < df2$time, 0, df2$event)
df2$time <- ifelse(df2$add_cens < df2$time, df2$add_cens, df2$time)

km2 <- survfit(Surv(time, event) ~ arm,
               data = df2)

p_km2 <- survminer::ggsurvplot(km2, 
                               data = df2, 
                               conf.int = TRUE,
                               risk.table = TRUE, 
                               break.x.by = 2,
                               legend.title = "",
                               legend.labs = c("0","1"),
                               xlab = "Time (months)",
                               ylab = "Overall survival",
                               risk.table.fontsize = 4,
                               legend = c(0.8,0.8))



pdf(here("output", "p_km2.pdf"))
print(p_km2, newpage = FALSE)
dev.off()
