### rmst1 code from survRM2 package
rmst1<-function(time, status, tau, alpha=0.05){
  #-- time
  #-- statuts
  #-- tau -- truncation time
  #-- alpha -- gives (1-alpha) confidence interval
  Surv<-survival::Surv
  
  ft= survival::survfit(Surv(time, status)~1)
  idx=ft$time<=tau
  
  wk.time=sort(c(ft$time[idx],tau))
  wk.surv=ft$surv[idx]
  wk.n.risk =ft$n.risk[idx]
  wk.n.event=ft$n.event[idx]
  
  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst = sum(areas)
  rmst
  
  wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                   wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
  wk.var =c(wk.var,0)
  rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se  = sqrt(rmst.var)
  
  #--- check ---
  # print(ft, rmean=tau)
  
  #--- output ---
  out=matrix(0,2,4)
  out[1,]=c(rmst, rmst.se, rmst-stats::qnorm(1-alpha/2)*rmst.se, rmst+stats::qnorm(1-alpha/2)*rmst.se)
  out[2,]=c(tau-out[1,1], rmst.se, tau-out[1,4], tau-out[1,3])
  rownames(out)=c("RMST","RMTL")
  colnames(out)=c("Est.", "se", paste("lower .",round((1-alpha)*100, digits=0), sep=""), 
                  paste("upper .",round((1-alpha)*100, digits=0), sep=""))
  
  Z=list()
  Z$result=out
  Z$rmst = out[1,]
  Z$rmtl = out[2,]
  Z$tau=tau
  Z$rmst.var = rmst.var
  Z$fit=ft
  class(Z)="rmst1"
  
  return(Z)
  
}


#################################################################################


find_pseudovalues_wmst <- function(formula,
                                   data,
                                   method,
                                   tau_1 = 1,
                                   tau_2 = 2){

  ### extract terms from formula
  nphRCT:::check_formula(formula=formula,data=data)
  
  formula_vars <- all.vars(formula)
  time_col <- formula_vars[1]
  status_col <- formula_vars[2]
  terms_vars<-formula_vars[-(1:2)]
  Terms <- stats::terms(formula,"strata")
  strata_index <- survival::untangle.specials(Terms,"strata")$terms
  
  if (length(strata_index)>0){
    strata_col<-terms_vars[strata_index]
    group_col<-terms_vars[-strata_index]  
  }else{
    group_col<-terms_vars
  }
  
  if (length(group_col)!=1) {
    stop("Formula must contain only one treatment arm indicator. All other terms must be specified as strata.")
  }
  
  data[[group_col]]<-as.factor(data[[group_col]])
  Surv<-survival::Surv
  
  rmst_full_tau_1 <- rmst1(time = data[[time_col]],
                           status = data[[status_col]],
                           tau = tau_1)$rmst[1]
  
  rmst_full_tau_2 <- rmst1(time = data[[time_col]],
                           status = data[[status_col]],
                           tau = tau_2)$rmst[1]
  
  wmst_full <- rmst_full_tau_2 - rmst_full_tau_1
  
  n <- nrow(data)
  
  ### also need the size of each arm...
  
  p_v_wmst <- numeric(nrow(data))
  
  for (i in seq_along(p_v_wmst)){
    data_minus_i <- data[-i,]
    
    rmst_minus_i_tau_1 <- rmst1(time = data_minus_i[[time_col]],
                                status = data_minus_i[[status_col]],
                                tau = tau_1)$rmst[1]
    
    rmst_minus_i_tau_2 <- rmst1(time = data_minus_i[[time_col]],
                                status = data_minus_i[[status_col]],
                                tau = tau_2)$rmst[1]
    
    wmst_minus_i <- rmst_minus_i_tau_2 - rmst_minus_i_tau_1
    
    p_v_wmst[i] <- n * wmst_full - (n - 1) * wmst_minus_i
    
  }
  max_p_s <- max(-p_v_wmst)
  min_p_s <- min(-p_v_wmst)
  
  A = 2 / (max_p_s - min_p_s)
  B = 1 - A * max_p_s
  
  df_wmst_pseudo <- data.frame(t_j = data[[time_col]],
                               event = data[[status_col]],
                               group = data[[group_col]],
                               score = -p_v_wmst,
                               standardized_score = -p_v_wmst * A + B)
  
  df_wmst_pseudo <- df_wmst_pseudo[order(df_wmst_pseudo[["t_j"]]),]
  
    
  out = list(df = df_wmst_pseudo)
  class(out) <- "df_score"
  out
 
}


#####################################################
find_pseudovalues_ahsw <- function(formula,
                                   data,
                                   method,
                                   tau){
  
  ### extract terms from formula
  nphRCT:::check_formula(formula=formula,data=data)
  
  formula_vars <- all.vars(formula)
  time_col <- formula_vars[1]
  status_col <- formula_vars[2]
  terms_vars<-formula_vars[-(1:2)]
  Terms <- stats::terms(formula,"strata")
  strata_index <- survival::untangle.specials(Terms,"strata")$terms
  
  if (length(strata_index)>0){
    strata_col<-terms_vars[strata_index]
    group_col<-terms_vars[-strata_index]  
  }else{
    group_col<-terms_vars
  }
  
  if (length(group_col)!=1) {
    stop("Formula must contain only one treatment arm indicator. All other terms must be specified as strata.")
  }
  
  data[[group_col]]<-as.factor(data[[group_col]])
  Surv<-survival::Surv
  
  rmst_full <- rmst1(time = data[[time_col]],
                           status = data[[status_col]],
                           tau = tau)$rmst[1]
  
 

  
  n <- length(data[[time_col]])
  
  ### non-parametric estimate of survival
  formula_km<-stats::as.formula(paste0("Surv(",time_col,",",status_col,")~1"))
  surv <- survival::survfit(formula_km, data = data)
  
  s_full <- summary(surv, time = tau)$surv
  
  ahsw_full <- log(rmst_full / (1 - s_full))
  
  
  ### also need the size of each arm...
  
  p_v_ahsw <- numeric(nrow(data))
  
  for (i in seq_along(p_v_ahsw)){
    data_minus_i <- data[-i,]
    
    rmst_minus_i <- rmst1(time = data_minus_i[[time_col]],
                                status = data_minus_i[[status_col]],
                                tau = tau)$rmst[1]
    
    
    ### non-parametric estimate of survival probability at tau
    formula_km<-stats::as.formula(paste0("Surv(",time_col,",",status_col,")~1"))
    surv_minus_i <- survival::survfit(formula_km, data = data_minus_i)
    s_full_minus_i <- summary(surv_minus_i, time = tau)$surv
    
    
    ahsw_minus_i <- log(rmst_minus_i / (1 - s_full_minus_i))
    
    p_v_ahsw[i] <- n * ahsw_full - (n - 1) * ahsw_minus_i
    
  }
  max_p_s <- max(-p_v_ahsw)
  min_p_s <- min(-p_v_ahsw)
  
  A = 2 / (max_p_s - min_p_s)
  B = 1 - A * max_p_s
  
  df_ahsw_pseudo <- data.frame(t_j = data[[time_col]],
                               event = data[[status_col]],
                               group = data[[group_col]],
                               score = -p_v_ahsw,
                               standardized_score = -p_v_ahsw * A + B)
  
  df_ahsw_pseudo <- df_ahsw_pseudo[order(df_ahsw_pseudo[["t_j"]]),]
  
  
  out = list(df = df_ahsw_pseudo)
  class(out) <- "df_score"
  out
  
}









