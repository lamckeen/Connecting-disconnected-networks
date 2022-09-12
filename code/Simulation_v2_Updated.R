# Simulation Procedure: Variance of muhat estimate
# Update: Including direct trial information
# Two large disconnected networks from BRD

library(tidyverse)
library(readxl)
library(gemtc)
library(igraph)
library(magic)
library(xtable)
library(tictoc)
source("Functions.R")

# Data --------------------------------------------------------------------

nw_1_dat = read.csv("Data/Network_1.csv")[,1:28]
nw_2_dat = read.csv("Data/Network_2.csv")[,1:28]
nw_combined_dat = rbind(nw_1_dat,nw_2_dat)

nw_combined_risks = rbind(read.csv("Data/Network_1_risks.csv"),
                          read.csv("Data/Network_2_risks.csv"))

network_1_trt = c("NAC","A-G","BOVA4PRE","A","H","I",
                  "L","M","P","T","Y","Z","Q","R","O",
                  "AUTONONAME3","PRESBOV")
network_2_trt = c("E", "E-O","E-P","E-Q","G","PMHP4")


# True Effects ------------------------------------------------------------

odds_l = nw_combined_risks[nw_combined_risks$trt=="L",2]/
  (1-nw_combined_risks[nw_combined_risks$trt=="L",2])
odds_g = nw_combined_risks[nw_combined_risks$trt=="G",2]/
  (1-nw_combined_risks[nw_combined_risks$trt=="G",2])

true_u_l_g = log(odds_g/odds_l)

odds_nac = nw_combined_risks[nw_combined_risks$trt=="NAC",2]/
  (1-nw_combined_risks[nw_combined_risks$trt=="NAC",2])
odds_e = nw_combined_risks[nw_combined_risks$trt=="E",2]/
  (1-nw_combined_risks[nw_combined_risks$trt=="E",2])

true_u_nac_e = log(odds_e/odds_nac)

# NMA Matrices ------------------------------------------------------------

get_mats = nma_mats(nw_combined_dat)
Y_mat = get_mats$Y
design_mat = get_mats$X
s_mat = get_mats$S
baseline = get_mats$baseline
main_trts = get_mats$main_trts
non_baseline_trts = get_mats$non_baseline_trts


# Trials ------------------------------------------------------------------

# trial_sim = data.frame(n = c(seq(3,10),25,50,100,250,500),
#                        trt_1 = c("PRESBOV", rep("H",6),rep("NAC",6)),
#                        trt_2 = rep("E",13))%>%
#     rbind(data.frame(n = c(seq(3,10),25,50,100,250,500),
#                      trt_1 = rep("NAC",13),
#                      trt_2 = rep("E", 13))) %>%
#     mutate(trial = paste0(trt_1, "-", trt_2))

# trial_sim = data.frame(n = c(seq(3,10),25,50,100,250,500),
#                        trt_1 = rep("NAC",13),
#                        trt_2 = rep("E",13)) %>%
#   mutate(trial = paste0(trt_1, "-", trt_2))

trial_sim = data.frame(n = c(seq(3,10),25,50,100,250,500),
                       trt_1 = c("PRESBOV", rep("H",6),rep("NAC",1),rep("L",5)),
                       trt_2 = rep("G",13)) %>%
  rbind(data.frame(n = c(seq(3,10),25,50,100,250,500),
                   trt_1 = rep("L",13),
                   trt_2 = rep("G", 13))) %>%
  mutate(trial = paste0(trt_1, "-", trt_2))

# Simulation --------------------------------------------------------------

## Apply a function to each row that simulates data and gets sampling var

var_tab = sapply(1:nrow(trial_sim), function(x){
  
  n_connect = trial_sim[x,"n"]
  x_new = rep(0, length(non_baseline_trts))
  baseline_arm = as.character(trial_sim[x,"trt_1"])
  other_arm = as.character(trial_sim[x,"trt_2"])
  
  if(baseline_arm == baseline){
    x_new[which(non_baseline_trts == other_arm)] <- 1
  }
  
  if(baseline_arm != baseline){
    trt1 <- which(non_baseline_trts == baseline_arm)
    trt2 <- which(non_baseline_trts == other_arm)
    x_new[trt1] <- -1
    x_new[trt2] <- 1
  }
  
  x_all = rbind(design_mat, x_new)
  
  p_1 = nw_combined_risks[which(nw_combined_risks$trt == baseline_arm), "risk_mean"]
  p_2 = nw_combined_risks[which(nw_combined_risks$trt == other_arm), "risk_mean"]
  
  u_l_g = c()
  l_g_rmse = c()
  var_l_g = c()
  # u_nac_e = c()
  # nac_e_rmse = c()
  # var_nac_e = c()
  
  for (i in 1:1000){
    
    # Simulate new trial data (and correct for 0)
    
    r_1 = rbinom(1,n_connect,p_1)
    r_2 = rbinom(1,n_connect,p_2)
    # r_1[r_1 == 0] <- .5
    # r_2[r_2 == 0] <- .5
    # r_1[r_1 == n_connect] <- n_connect-.5
    # r_2[r_2 == n_connect] <- n_connect-.5
    # 
    # p_1_hat = r_1/n_connect
    # p_2_hat = r_2/n_connect
    
    p_1_hat = (r_1+1)/(n_connect+2)
    p_2_hat = (r_2+1)/(n_connect+2)
    
    odds_1 = p_1_hat/(1-p_1_hat)
    odds_2 = p_2_hat/(1-p_2_hat)
    
    logodds_12 = log(odds_2/odds_1) #new y value
    
    sigma2_new = 1/(n_connect*p_1_hat*(1-p_1_hat)) + 1/(n_connect*p_2_hat*(1-p_2_hat))
    
    # Calculate new mu hat vector
    
    ## New S
    s_all = rbind(s_mat,c(rep(0,nrow(s_mat)))) %>%
      cbind(c(rep(0,nrow(.))))
    
    s_all[nrow(s_all),nrow(s_all)] <- sigma2_new
    
    ## New Y
    
    y_all = c(Y_mat,logodds_12)
    
    ## Mu hat
    
    mu_hat = solve(t(x_all)%*%solve(s_all)%*%x_all) %*% t(x_all)%*%solve(s_all)%*% y_all
    names(mu_hat) = non_baseline_trts
    
    # Save difference for specific contrast
    
    # u_nac_e[i] = mu_hat[3]
    # nac_e_rmse[i] = (true_u_nac_e-u_nac_e[i])^2
    u_l_g[i] = mu_hat[15]-mu_hat[8]
    l_g_rmse[i] = (true_u_l_g-u_l_g[i])^2
    
    
    # Save variance for contrast
    var_mu = solve(t(x_all)%*%solve(s_all)%*%x_all)
    c_vec = rep(0,length(non_baseline_trts))
    c_vec[15] = 1
    c_vec[8] = -1
    var_l_g[i] = t(c_vec) %*% var_mu %*% (c_vec)
  }
  
  # true_trial_sigma2 = 1/(n_connect*p_1*(1-p_1)) + 1/(n_connect*p_2*(1-p_2))
  
  return(list(vv = var(u_l_g),
              mm = mean(u_l_g),
              rmse=mean(l_g_rmse),
              v_mean = mean(var_l_g)))
  
  
  
})


# Summary Table -----------------------------------------------------------

# calculated_var = c(
#   2.8260,
#   2.1221,
#   1.6988,
#   1.4166,
#   1.2151,
#   1.0639,
#   0.9463 ,
#   0.8519,
#   0.3408,
#   0.1704,
#   0.0852,
#   0.0341,
#   0.0170
# )

calculated_var=c(
  4.5914,
  3.7483,
  3.2414,
  2.9035,
  2.6621,
  2.4811,
  2.3403,
  2.2273,
  1.4054,
  0.7027,
  0.3513,
  0.1405,
  0.0703
)

tt = var_tab %>%
  as.data.frame() %>% 
  t() %>%
  data.frame() %>%
  mutate(n = 2*trial_sim$n,
         trial = trial_sim$trial,
         bias = as.numeric(mm)-true_u_l_g)

tt_1 = tt %>%
  head(13) %>%
  left_join(tt %>% tail(13), by="n") %>%
  mutate(calc.x = calculated_var,
         calc.y = c(
                    # 11.7115,
                    # 8.7837,
                    # 7.0269,
                    # 5.8558,
                    # 5.0192,
                    # 4.3918,
                    # 3.9038,
                    # 3.5135,
                    rep(NA,13)),
         blank = rep(NA,13)) %>%
  dplyr::select(c(n, trial.x,vv.x,v_mean.x,calc.x,bias.x,rmse.x,blank,
                  vv.y,v_mean.y,calc.y,bias.y,rmse.y))

xtable(tt_1, digits=c(0,0,0,4,4,4,4,4,4,4,4,4,4,4))
