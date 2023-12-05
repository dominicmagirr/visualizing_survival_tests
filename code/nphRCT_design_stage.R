library(dplyr)
## Suppose we anticipate on the control arm that survival 
## is exponentially distributed with median 12 months.
## On the experimental arm, we expect survival to follow the 
## same distribution as the control arm for the first six months.
## Thereafter, we expect the HR to be 0.6
## Also, we expect uniform recruitment over 12 months, with a further
## 24 months of follow-up...

N_1 <- 200
N_0 <- 200

G_inverse <- function(u, min_follow_up, trial_end){
  min_follow_up + u * (trial_end - min_follow_up)
}


S_inverse <- function(u, lambda_c, hr, t_star){
  
 u_star <- exp(-lambda_c * t_star)
 ifelse(u > u_star, 
        -log(u) / lambda_c,
        -log(u / u_star) / (lambda_c * hr) + t_star)
  
}

times_1_unc <- S_inverse((1:N_1)/N_1, log(2) / 12, 0.6, 6)
times_0_unc <- S_inverse((1:N_0)/N_0, log(2) / 12, 1, 6)

C_1 <- G_inverse((1:N_1)/N_1, 24, 36)
C_0 <- G_inverse((1:N_0)/N_0, 24, 36)

times_1 <- pmin(times_1_unc, C_1)
times_0 <- pmin(times_0_unc, C_0)

plot(times_1, (1:N_1)/N_1)

event_1 <- times_1_unc < C_1
event_0 <- times_0_unc < C_0

dat <- data.frame(time = c(times_1, times_0),
                  event = c(event_1, event_0),
                  arm = rep(c("1", "0"), c(N_1, N_0)))

p_1 <- nphRCT::find_scores(Surv(time, event) ~ arm, 
                    data = dat,
                    method = "mw",
                    s_star = 1) %>% plot() 

p_2 <- nphRCT::find_scores(Surv(time, event) ~ arm, 
                    tau=24, data=dat,
                    method = "rmst") %>% plot()

cowplot::plot_grid(p_1, p_2)

