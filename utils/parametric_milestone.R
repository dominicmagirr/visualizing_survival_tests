##################################
## parametric milestone
##################################
get_scores_surv_p <- function(df, tau, breaks = c(0,2,4,6,8,100)){
  
  n <- length(df$time)
  surv_both <- pchreg(Surv(time,event)~1, data = df, breaks = breaks)
  s_full <- predict(surv_both, newdata = data.frame(time = tau, event = 1))$Surv
  p_v_surv <- numeric(length(df$time))
  
  for (i in seq_along(p_v_surv)){
    
    df_minus_i <- df[-i,]
    
    ### parametric estimate of survival probability at tau
    surv_minus_i <- pchreg(Surv(time, event)~1, data = df_minus_i, breaks = breaks)
    s_full_minus_i <- predict(surv_minus_i, newdata = data.frame(time = tau, event = 1))$Surv
    
    p_v_surv[i] <- n * s_full - (n - 1) * s_full_minus_i
    
  }
  max_p_s <- max(-p_v_surv)
  min_p_s <- min(-p_v_surv)
  
  A = 2 / (max_p_s - min_p_s)
  B = 1 - A * max_p_s
  
  df_surv_pseudo <- data.frame(t_j = df$time,
                               event = df$event,
                               group = factor(df$arm),
                               score = -p_v_surv,
                               standardized_score = -p_v_surv * A + B)
  
  out = list(df = df_surv_pseudo)
  class(out) <- "df_score"
  out
}

##################################
## flexible parametric milestone
##################################
get_scores_flexsurv_p <- function(df, tau, k = 3){
  
  n <- length(df$time)
  
  ### flex parametric
  surv_flex <- flexsurvspline(Surv(time, event) ~ 1, data = df, k = k, scale = "hazard")
  
  s_full <- predict(surv_flex,
                    newdata = data.frame(time = tau, event = 1), 
                    type = "survival",
                    times = tau)$.pred_survival
  
  p_v_surv <- numeric(length(df$time))
  
  for (i in seq_along(p_v_surv)){
    
    df_minus_i <- df[-i,]
    
    ### parametric estimate of survival probability at tau
    surv_minus_i <- flexsurvspline(Surv(time, event) ~ 1, 
                                   data = df_minus_i, 
                                   scale = "hazard",
                                   knots = surv_flex$knots[-c(1,k+2)],
                                   bknots = surv_flex$knots[c(1,k+2)])
    
    s_full_minus_i <- predict(surv_minus_i,
                              newdata = data.frame(time = tau, event = 1), 
                              type = "survival",
                              times = tau)$.pred_survival
    
    p_v_surv[i] <- n * s_full - (n - 1) * s_full_minus_i
    
  }
  max_p_s <- max(-p_v_surv)
  min_p_s <- min(-p_v_surv)
  
  A = 2 / (max_p_s - min_p_s)
  B = 1 - A * max_p_s
  
  df_surv_pseudo <- data.frame(t_j = df$time,
                               event = df$event,
                               group = factor(df$arm),
                               score = -p_v_surv,
                               standardized_score = -p_v_surv * A + B)
  
  out = list(df = df_surv_pseudo)
  class(out) <- "df_score"
  out
}
