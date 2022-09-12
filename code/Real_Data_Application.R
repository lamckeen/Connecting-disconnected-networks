# Real Data Application: Finding the best connecting trial
# Two large disconnected networks from BRD

library(tidyverse)
library(readxl)
library(gemtc)
library(igraph)
library(magic)
library(xtable)
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

trials_connect = expand.grid(network_1_trt,network_2_trt)

n_connect = 3 #Sample size in connecting trial

# Get Variance Matrix -----------------------------------------------------

get_mats = nma_mats(nw_combined_dat)
design_mat = get_mats$X
s_mat = get_mats$S
baseline = get_mats$baseline
main_trts = get_mats$main_trts
non_baseline_trts = get_mats$non_baseline_trts

connect_var_all = sapply(1:nrow(trials_connect), function(x){
  
  ## Get new row in design and new design
  
  x_new = rep(0, length(non_baseline_trts))
  baseline_arm = as.character(trials_connect[x,"Var1"])
  other_arm = as.character(trials_connect[x,"Var2"])
  
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
  
  ## Get new element in s and new s
  
  p_1 = nw_combined_risks[which(nw_combined_risks$trt == baseline_arm), "risk_mean"]
  p_2 = nw_combined_risks[which(nw_combined_risks$trt == other_arm), "risk_mean"]
  sigma2_new = 1/(n_connect*p_1*(1-p_1)) + 1/(n_connect*p_2*(1-p_2))
  
  s_all = rbind(s_mat,c(rep(0,nrow(s_mat)))) %>%
    cbind(c(rep(0,nrow(.))))
  
  s_all[nrow(s_all),nrow(s_all)] <- sigma2_new
                
  
  ## Get variance of mu 
  
  var_mu = solve(t(x_all)%*%solve(s_all)%*%x_all)
  rownames(var_mu) = non_baseline_trts
  colnames(var_mu) = non_baseline_trts
  
  return(var_mu)
  
}, simplify=FALSE)


# Get Contrast Variance ---------------------------------------------------

## Get variance of contrast of interest - start with NAC to E

coi_NAC_E = lapply(connect_var_all, function(x){
  c_vec = rep(0,length(non_baseline_trts))
  c_vec[3] = 1 #starting with this specific cvec
  coi_var = t(c_vec) %*% x %*% (c_vec)
  return(as.numeric(coi_var))
})

coi_NAC_E_v = do.call(rbind,coi_NAC_E)
trials_connect[which.min(coi_NAC_E_v),] # minimum variance

## All contrasts...use my trials connect to get all c vectors for all contrasts!
## Then for each 102 elements of my list, get all 102 variances of contrasts

## Finding a trial that is completely indirect? None of those for NAC-E
trials_connect[which((coi_NAC_E_v)<coi_NAC_E_v[[1]]),]

## For comparison of interest L-G, there do exist those that are better than L-G directly...

coi_test = lapply(connect_var_all, function(x){
  c_vec = rep(0,length(non_baseline_trts))
  c_vec[15] = -1 #starting with this specific cvec
  c_vec[8] = 1 #starting with this specific cvec
  coi_var = t(c_vec) %*% x %*% (c_vec)
  return(as.numeric(coi_var))})

coi_test_v = do.call(rbind,coi_test) %>%
  as.data.frame() %>%
  cbind(as.character(trials_connect$Var1),as.character(trials_connect$Var2))
colnames(coi_test_v) = c("Var", "T1", "T2")
coi_L_G = (coi_test_v %>% filter(T1 == "L" & T2 == "G"))[1,1]
tt = coi_test_v %>%
  filter(Var < coi_L_G) %>%
  filter(T1 != "L") %>%
  filter(T2 != "G")
tt[which.min(tt$Var),]

## Now varying sample size to illustrate this indirect
  
n_connect_all = c(seq(3,10),25,50,100,250,500)
best_trial = list()

for (i in 1:length(n_connect_all)){
  
  n_connect = n_connect_all[i]
  connect_var_all = sapply(1:nrow(trials_connect), function(x){
    
    ## Get new row in design and new design
    
    x_new = rep(0, length(non_baseline_trts))
    baseline_arm = as.character(trials_connect[x,"Var1"])
    other_arm = as.character(trials_connect[x,"Var2"])
    
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
    
    ## Get new element in s and new s
    
    p_1 = nw_combined_risks[which(nw_combined_risks$trt == baseline_arm), "risk_mean"]
    p_2 = nw_combined_risks[which(nw_combined_risks$trt == other_arm), "risk_mean"]
    sigma2_new = 1/(n_connect*p_1*(1-p_1)) + 1/(n_connect*p_2*(1-p_2))
    
    s_all = rbind(s_mat,c(rep(0,nrow(s_mat)))) %>%
      cbind(c(rep(0,nrow(.))))
    
    s_all[nrow(s_all),nrow(s_all)] <- sigma2_new
    
    
    ## Get variance of mu 
    
    var_mu = solve(t(x_all)%*%solve(s_all)%*%x_all)
    rownames(var_mu) = non_baseline_trts
    colnames(var_mu) = non_baseline_trts
    
    return(var_mu)
    
  }, simplify=FALSE)
  
  coi_L_G = lapply(connect_var_all, function(x){
    c_vec = rep(0,length(non_baseline_trts))
    c_vec[15] = -1 #starting with this specific cvec
    c_vec[8] = 1 #starting with this specific cvec
    coi_var = t(c_vec) %*% x %*% (c_vec)
    return(as.numeric(coi_var))
  })
  
  coi_test_v = do.call(rbind,coi_L_G) %>%
    as.data.frame() %>%
    cbind(as.character(trials_connect$Var1),as.character(trials_connect$Var2))
  colnames(coi_test_v) = c("Var", "T1", "T2")
  coi_L_G_num = (coi_test_v %>% filter(T1 == "L" & T2 == "G"))[1,1]
  tt = coi_test_v %>%
    filter(Var < coi_L_G_num) %>%
    filter(T1 != "L") %>%
    filter(T2 != "G")
  if (nrow(tt) == 0){
    min_best = as.data.frame(t(c(Var = coi_L_G_num,T1 = "L",T2="G")))
    min_best$Var = as.numeric(min_best$Var)
  }else{
    min_best = tt[which.min(tt$Var),]
  }
  
  best_trial[[i]] = c(2*n_connect,
                      paste0(min_best$T1,"-",
                             min_best$T2),
                      min_best$Var,
                      coi_L_G_num)
}

aa = do.call(rbind,best_trial) %>% as.data.frame()
aa[,3] = round(as.numeric(aa[,3]),4)
aa[,4] = round(as.numeric(aa[,4]),4)
xtable(aa,digits=4)
  

# Varying Sample Size -----------------------------------------------------

n_connect_all = c(seq(3,10),25,50,100,250,500)
best_trial = list()

for (i in 1:length(n_connect_all)){
  
  n_connect = n_connect_all[i]
  connect_var_all = sapply(1:nrow(trials_connect), function(x){
    
    ## Get new row in design and new design
    
    x_new = rep(0, length(non_baseline_trts))
    baseline_arm = as.character(trials_connect[x,"Var1"])
    other_arm = as.character(trials_connect[x,"Var2"])
    
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
    
    ## Get new element in s and new s
    
    p_1 = nw_combined_risks[which(nw_combined_risks$trt == baseline_arm), "risk_mean"]
    p_2 = nw_combined_risks[which(nw_combined_risks$trt == other_arm), "risk_mean"]
    sigma2_new = 1/(n_connect*p_1*(1-p_1)) + 1/(n_connect*p_2*(1-p_2))
    
    s_all = rbind(s_mat,c(rep(0,nrow(s_mat)))) %>%
      cbind(c(rep(0,nrow(.))))
    
    s_all[nrow(s_all),nrow(s_all)] <- sigma2_new
    
    
    ## Get variance of mu 
    
    var_mu = solve(t(x_all)%*%solve(s_all)%*%x_all)
    rownames(var_mu) = non_baseline_trts
    colnames(var_mu) = non_baseline_trts
    
    return(var_mu)
    
  }, simplify=FALSE)
  coi_NAC_E = lapply(connect_var_all, function(x){
    c_vec = rep(0,length(non_baseline_trts))
    c_vec[3] = 1 #starting with this specific cvec
    coi_var = t(c_vec) %*% x %*% (c_vec)
    return(as.numeric(coi_var))
  })
  
  coi_NAC_E_v = do.call(rbind,coi_NAC_E)
  p_NAC = nw_combined_risks[9,2]
  p_E = nw_combined_risks[18,2]
  best_trial[[i]] = c(2*n_connect,
                     paste0(trials_connect[which.min(coi_NAC_E_v),][,1],"-",
                           trials_connect[which.min(coi_NAC_E_v),][,2]),
                      round(min(coi_NAC_E_v),4),
                     round(1/(n_connect*p_NAC*(1-p_NAC))+1/(n_connect*p_E*(1-p_E)),4)
                     )
}

do.call(rbind,best_trial)
xtable(do.call(rbind,best_trial))

# Varying Sample Size -----------------------------------------------------
# Quick check of L-G

n_connect_all = c(seq(3,10),25,50,100,250,500)
best_trial = list()

for (i in 1:length(n_connect_all)){
  
  n_connect = n_connect_all[i]
  connect_var_all = sapply(1:nrow(trials_connect), function(x){
    
    ## Get new row in design and new design
    
    x_new = rep(0, length(non_baseline_trts))
    baseline_arm = as.character(trials_connect[x,"Var1"])
    other_arm = as.character(trials_connect[x,"Var2"])
    
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
    
    ## Get new element in s and new s
    
    p_1 = nw_combined_risks[which(nw_combined_risks$trt == baseline_arm), "risk_mean"]
    p_2 = nw_combined_risks[which(nw_combined_risks$trt == other_arm), "risk_mean"]
    sigma2_new = 1/(n_connect*p_1*(1-p_1)) + 1/(n_connect*p_2*(1-p_2))
    
    s_all = rbind(s_mat,c(rep(0,nrow(s_mat)))) %>%
      cbind(c(rep(0,nrow(.))))
    
    s_all[nrow(s_all),nrow(s_all)] <- sigma2_new
    
    
    ## Get variance of mu 
    
    var_mu = solve(t(x_all)%*%solve(s_all)%*%x_all)
    rownames(var_mu) = non_baseline_trts
    colnames(var_mu) = non_baseline_trts
    
    return(var_mu)
    
  }, simplify=FALSE)
  coi_L_G = lapply(connect_var_all, function(x){
    c_vec = rep(0,length(non_baseline_trts))
    c_vec[15] = 1 #starting with this specific cvec
    c_vec[8] = -1 #starting with this specific cvec
    coi_var = t(c_vec) %*% x %*% (c_vec)
    return(as.numeric(coi_var))
  })
  
  p_L = nw_combined_risks[7,2]
  p_G = nw_combined_risks[22,2]
  coi_L_G_v = do.call(rbind,coi_L_G)
  best_trial[[i]] = c(2*n_connect,
                      paste0(trials_connect[which.min(coi_L_G_v),][,1],"-",
                             trials_connect[which.min(coi_L_G_v),][,2]),
                      round(min(coi_L_G_v),4),
                      round(1/(n_connect*p_L*(1-p_L))+1/(n_connect*p_G*(1-p_G)),4))
}

aa = do.call(rbind,best_trial) %>% as.data.frame()
aa[,3] = round(as.numeric(aa[,3]),4)
xtable(aa,digits=4)


# Example: COI XX-YY --------------------------------------------------------

#EO-A, T-G, R-G(always best direct...),NAC-G,NAC-EP(bigger diff) < *NAC-EQ,H-EQ*
# I-EQ is found as biggest difference in all combinations script

n_connect_all = c(seq(3,10),25,50,100,250,500)
best_trial = list()

for (i in 1:length(n_connect_all)){
  
  n_connect = n_connect_all[i]
  connect_var_all = sapply(1:nrow(trials_connect), function(x){
    
    ## Get new row in design and new design
    
    x_new = rep(0, length(non_baseline_trts))
    baseline_arm = as.character(trials_connect[x,"Var1"])
    other_arm = as.character(trials_connect[x,"Var2"])
    
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
    
    ## Get new element in s and new s
    
    p_1 = nw_combined_risks[which(nw_combined_risks$trt == baseline_arm), "risk_mean"]
    p_2 = nw_combined_risks[which(nw_combined_risks$trt == other_arm), "risk_mean"]
    sigma2_new = 1/(n_connect*p_1*(1-p_1)) + 1/(n_connect*p_2*(1-p_2))
    
    s_all = rbind(s_mat,c(rep(0,nrow(s_mat)))) %>%
      cbind(c(rep(0,nrow(.))))
    
    s_all[nrow(s_all),nrow(s_all)] <- sigma2_new
    
    
    ## Get variance of mu 
    
    var_mu = solve(t(x_all)%*%solve(s_all)%*%x_all)
    rownames(var_mu) = non_baseline_trts
    colnames(var_mu) = non_baseline_trts
    
    return(var_mu)
    
  }, simplify=FALSE)
  coi_I_EQ = lapply(connect_var_all, function(x){
    c_vec = rep(0,length(non_baseline_trts))
    c_vec[22] = 1 #starting with this specific cvec
    c_vec[7] = -1
    coi_var = t(c_vec) %*% x %*% (c_vec)
    return(as.numeric(coi_var))
  })
  
  coi_I_EQ_v = do.call(rbind,coi_I_EQ)
  p_I = nw_combined_risks[6,2]
  p_EQ = nw_combined_risks[21,2]
  best_trial[[i]] = c(2*n_connect,
                      paste0(trials_connect[which.min(coi_I_EQ_v),][,1],"-",
                             trials_connect[which.min(coi_I_EQ_v),][,2]),
                      round(as.numeric(min(coi_I_EQ_v)),4),
                      round(1/(n_connect*p_I*(1-p_I))+1/(n_connect*p_EQ*(1-p_EQ)),4)
  )
}

do.call(rbind,best_trial)
xtable(do.call(rbind,best_trial))



# Formula Verification ----------------------------------------------------

## Connecting networks one and two through NAC-E with n=500each (simulated data version)

n_connect=500
x_new = rep(0, length(non_baseline_trts))
baseline_arm = as.character(trials_connect[1,"Var1"])
other_arm = as.character(trials_connect[1,"Var2"])
x_new[which(non_baseline_trts == other_arm)] <- 1

x_all = rbind(design_mat, x_new)

p_1 = nw_combined_risks[which(nw_combined_risks$trt == baseline_arm), "risk_mean"]
p_2 = nw_combined_risks[which(nw_combined_risks$trt == other_arm), "risk_mean"]
r_1 = rbinom(1,n_connect,p_1)
r_2 = rbinom(1,n_connect,p_2)
if (r_1 == 0){r_1 = .5}
if (r_2 == 0){r_2 = .5}
if (r_1 == n_connect){r_1 = n_connect-.5}
if (r_2 == n_connect){r_2 = n_connect-.5}

p_1_hat = r_1/n_connect
p_2_hat = r_2/n_connect
sigma2_new = 1/(n_connect*p_1_hat*(1-p_1_hat)) + 1/(n_connect*p_2_hat*(1-p_2_hat))

s_all = rbind(s_mat,c(rep(0,nrow(s_mat)))) %>%
  cbind(c(rep(0,nrow(.))))

s_all[nrow(s_all),nrow(s_all)] <- sigma2_new

var_mu = solve(t(x_all)%*%solve(s_all)%*%x_all)
rownames(var_mu) = non_baseline_trts
colnames(var_mu) = non_baseline_trts

c_vec = rep(0,length(non_baseline_trts))
c_vec[22] = 1
c_vec[4] = -1
t(c_vec) %*% var_mu %*% (c_vec)
