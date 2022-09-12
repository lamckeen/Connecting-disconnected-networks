# All functions for NMA

##  Wide to long function to get contrast data into form we need for netmeta
## Works for 4 arms max

wide_to_long_contrast <- function(MTCdata){
  N <- max(MTCdata$Number.of.arms)
  study <- rep(MTCdata$study,N-1)
  trt2 <- NULL
  se <- NULL
  lor <- NULL
  v <- NULL
  num_arm <- NULL
  for(i in c(2,3,4)){
    lor <- c(lor, eval(parse(text = paste0("MTCdata$lor.",i, sep = ""))))
    se <- c(se, eval(parse(text = paste0("MTCdata$se.",i, sep = ""))))
    trt2 <- c(trt2, eval(parse(text = paste0("MTCdata$Arm.",i, sep = ""))))
    v <- c(v, MTCdata$V)
    num_arm <- c(num_arm, MTCdata$Number.of.arms)
  }
  res <- data.frame(id = study,trt2=trt2,lor=lor,se=se,v=v,num_arm=num_arm)
  res <- res %>% dplyr::filter(!is.na(trt2)) %>% arrange(id) %>%
    mutate(base = sapply(1:nrow(.),function(x){
      MTCdata[which(MTCdata$study == .[x,1]),2]}))
  
  # Getting comparisons for multiarm
  n_unique = length(unique(res$id))
  add_rows = list()
  
  for (j in 1:n_unique){
    stud = unique(res$id)[j]
    df = res %>% filter(id == stud)
    
    if (df$num_arm[1] == 2){
      add_rows[j] = NULL
    }else if (df$num_arm[1]==3){
      # Add 1 row
      t_1 = df$trt2[1]
      t_2 = df$trt2[2]
      
      add_rows[[j]] = data.frame(id = unique(df$id), 
                             trt2 = t_2, 
                             lor = df[2,"lor"]-df[1,"lor"], 
                             se = sqrt(df[2,"se"]^2+df[1,"se"]^2-2*df[2,"v"]), 
                             v = NA, 
                             base =t_1)
      
    }else if (df$num_arm[1]==4){
      # Add 3 rows
      t_1 = df$trt2[1]
      t_2 = df$trt2[2]
      t_3 = df$trt2[3]
      
      add_row_1 = data.frame(id = unique(df$id), 
                             trt2 = t_2, 
                             lor = df[2,"lor"]-df[1,"lor"], 
                             se = sqrt(df[2,"se"]^2+df[1,"se"]^2-2*df[2,"v"]), 
                             v = NA, 
                             base =t_1)
      
      add_row_2 = data.frame(id = unique(df$id), 
                             trt2 = t_3, 
                             lor = df[3,"lor"]-df[1,"lor"], 
                             se = sqrt(df[3,"se"]^2+df[1,"se"]^2-2*df[3,"v"]), 
                             v = NA, 
                             base =t_1)
      
      add_row_3 = data.frame(id = unique(df$id), 
                             trt2 = t_3, 
                             lor = df[3,"lor"]-df[2,"lor"], 
                             se = sqrt(df[3,"se"]^2+df[2,"se"]^2-2*df[2,"v"]), 
                             v = NA, 
                             base =t_2)
      add_rows[[j]] = rbind(add_row_1,add_row_2,add_row_3)
      
    }
  }
  
  add_rows_all = do.call(rbind,add_rows)
  res_final = rbind(res[,-6],add_rows_all)
  return(res_final)
}



## Function to get design/error matrices to use for frequentist NMA

nma_mats <- function(MTC_data){
  
  main_trts= setdiff(unique(c(MTC_data$Arm.1,MTC_data$Arm.2,
                                 MTC_data$Arm.3,MTC_data$Arm.4)),NA)
  baseline = names(which(table(MTC_data$Arm.1) == max(table(MTC_data$Arm.1))))
  non_baseline_trts = setdiff(main_trts, baseline)
  
  data <- MTC_data
  Y <- list()
  S <- list()
  X <- list()
  
  for(i in 1:length(data)){
    
    n_arms <- sum(!is.na(with(data[i, ], c(Arm.1, Arm.2, Arm.3, Arm.4))))
    
    if(n_arms == 2){
      y <- data[i, "lor.2"]
      s <- (data[i, "se.2"])^2
      
      x <- rep(0, length(non_baseline_trts))
      baseline_arm <- data[i, "Arm.1"]
      if(baseline_arm == baseline){
        x[which(non_baseline_trts == data[i, "Arm.2"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.1"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.2"])
        x[trt1] <- -1
        x[trt2] <- 1
      }
      
      Y[[i]] <- y
      S[[i]] <- s
      X[[i]] <- x
    }
    
    if(n_arms == 3){
      y1 <- data[i, "lor.2"]
      s1 <- (data[i, "se.2"])^2
      
      y2 <- data[i, "lor.3"]
      s2 <- (data[i, "se.3"])^2
      
      v <- data[i, "V"]
      
      x1 <- rep(0, length(non_baseline_trts))
      x2 <- rep(0, length(non_baseline_trts))
      
      baseline_arm <- data[i, "Arm.1"]
      if(baseline_arm == baseline){
        x1[which(non_baseline_trts == data[i, "Arm.2"])] <- 1
        x2[which(non_baseline_trts == data[i, "Arm.3"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.1"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.2"])
        trt3 <- which(non_baseline_trts == data[i, "Arm.3"])
        x1[trt1] <- -1
        x1[trt2] <- 1
        x2[trt1] <- -1
        x2[trt3] <- 1
      }
      
      Y[[i]] <- c(y1, y2)
      S[[i]] <- matrix(data = c(s1, v, v, s2), nrow = 2, byrow = TRUE)
      X[[i]] <- rbind(x1, x2)
    }
    
    if(n_arms == 4){
      y1 <- data[i, "lor.2"]
      s1 <- (data[i, "se.2"])^2
      
      y2 <- data[i, "lor.3"]
      s2 <- (data[i, "se.3"])^2
      
      y3 <- data[i, "lor.4"]
      s3 <- (data[i, "se.4"])^2
      
      v <- data[i, "V"]
      
      x1 <- rep(0, length(non_baseline_trts))
      x2 <- rep(0, length(non_baseline_trts))
      x3 <- rep(0, length(non_baseline_trts))
      
      baseline_arm <- data[i, "Arm.1"]
      if(baseline_arm == baseline){
        x1[which(non_baseline_trts == data[i, "Arm.2"])] <- 1
        x2[which(non_baseline_trts == data[i, "Arm.3"])] <- 1
        x3[which(non_baseline_trts == data[i, "Arm.4"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.1"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.2"])
        trt3 <- which(non_baseline_trts == data[i, "Arm.3"])
        trt4 <- which(non_baseline_trts == data[i, "Arm.4"])
        x1[trt1] <- -1
        x1[trt2] <- 1
        x2[trt1] <- -1
        x2[trt3] <- 1
        x3[trt1] <- -1
        x3[trt4] <- 1
      }
      
      Y[[i]] <- c(y1, y2, y3)
      S[[i]] <- matrix(data = c(s1, v, v, 
                                v, s2, v, 
                                v, v, s3), nrow = 3, byrow = TRUE)
      X[[i]] <- rbind(x1, x2, x3)
    }
  }  
  
  X_final <- do.call(rbind, X)
  
  Y_final <- do.call(c, Y)
  
  S_final <- do.call(adiag, S)
  
  output <- list(Y = Y_final, S = S_final, X = X_final, 
                 main_trts = main_trts, baseline = baseline,
                 non_baseline_trts = non_baseline_trts)
  return(output)
} 
