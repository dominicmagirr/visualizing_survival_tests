here::i_am("code/POPLAR_trial.R")
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

km <- survfit(Surv(time, event) ~ arm,
              data = df)

p_km <- survminer::ggsurvplot(km, 
                              data = df, 
                              conf.int = TRUE,
                              risk.table = TRUE, 
                              break.x.by = 2,
                              legend.title = "",
                              legend.labs = c("Docetaxel","Atezolizumab"),
                              xlab = "Time (months)",
                              ylab = "Overall survival",
                              risk.table.fontsize = 4,
                              legend = c(0.8,0.8))

pdf(here("output", "p_km.pdf"))
print(p_km, newpage = FALSE)
dev.off()
