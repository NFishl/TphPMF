test_generation_test_thresh <- function(control_imp_mat, crc_imp_mat, imp_method, test_thresh){
  sum(abs(control_imp_mat - sim_tab_zi_control) < 0.00001)/(300 * 60)
  sum(abs(control_imp_mat - log10(1.01)) < 0.00001)/(300 * 60)
  
  sum(abs(crc_imp_mat - sim_tab_zi_CRC) < 0.00001)/(300 * 60)
  sum(abs(crc_imp_mat - log10(1.01)) < 0.00001)/(300 * 60)
  
  
  control_imp_mat[1:6,1:6]
  imputed_mat <- sim_tab_zi
  imputed_mat[which(condition == 1), ] = as.matrix(crc_imp_mat)
  imputed_mat[which(condition == 0), ] = as.matrix(control_imp_mat)
  
  # sum(((sim_tab - imputed_mat)^2))
  
  #explore the tsne plot, PCA and heatmap
  pdf(file = paste0("plots/", imp_method, ".pdf", collapse = ""))
  # sim_tsne_out <- tsne(imputed_mat)
  # plot(sim_tsne_out, col = as.factor(condition))
  
  sim_rtsne_out <- Rtsne(imputed_mat)
  plot(sim_rtsne_out$Y, col = as.factor(condition))
  
  sim_pca_out <- prcomp(imputed_mat, center = TRUE)
  plot(sim_pca_out$x[,1:2], col = as.factor(condition))
  
  heatmap(as.matrix(imputed_mat), scale="row")
  
  distmat <- as.matrix(dist(as.matrix(imputed_mat)))
  diag(distmat) = max(distmat)
  heatmap(distmat)
  dev.off()
  #now check the wilcoxon test results
  less_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
  greater_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
  # less_identified_imp <- which(p.adjust(less_pval_imp, method = "fdr") < test_thresh)
  # greater_identified_imp <- which(p.adjust(greater_pval_imp, method = "fdr") < test_thresh)
  # 
  less_identified_imp <- which(less_pval_imp < test_thresh)
  greater_identified_imp <- which(greater_pval_imp < test_thresh)
  
  
  sum(less_identified_imp %in% less_identified_full)
  sum(greater_identified_imp %in% greater_identified_full)
  
  length(less_identified_imp)
  length(greater_identified_imp)
  
  length(less_identified_full)
  length(greater_identified_full)
  
  #make a .csv file that contains all the comparisons between t chosen taxa
  result <- matrix(NA, nrow = 3, ncol = 4)
  result[1,1] = length(less_identified_full)
  result[1,2] = length(greater_identified_full)
  result[2,1] = length(less_identified_zi)
  result[2,2] = length(greater_identified_zi)
  result[3,1] = length(less_identified_imp)
  result[3,2] = length(greater_identified_imp)
  result[2,3] = sum(less_identified_zi %in% less_identified_full)
  result[2,4] = sum(greater_identified_zi %in% greater_identified_full)
  result[3,3] = sum(less_identified_imp %in% less_identified_full)
  result[3,4] = sum(greater_identified_imp %in% greater_identified_full)
  # result[2,5] = sum(less_identified_zi %in% less_identified_imp)
  # result[3,5] = sum(greater_identified_zi %in% greater_identified_imp)
  
  colnames(result) <- c("# of taxa in CRC < control", "# of taxa in CRC > control", "overlap with full data in <", "overlap with full data in >")
  rownames(result) <- c("Full", "ZI", "Imputed")
  write.csv(result, file = paste0("results/", imp_method, ".csv", collapse = ""))
  return(list("precision" = (result[3,3] + result[3,4])/(result[3,1] + result[3,2]), "recall" = (result[3,3] + result[3,4])/(result[1,1] + result[1,2]) ))
}


test_generation_test_thresh_raw <- function(control_imp_mat, crc_imp_mat, imp_method, test_thresh){
  sum(abs(control_imp_mat - sim_tab_zi_control) < 0.00001)/(300 * 60)
  sum(abs(control_imp_mat - log10(1.01)) < 0.00001)/(300 * 60)
  
  sum(abs(crc_imp_mat - sim_tab_zi_CRC) < 0.00001)/(300 * 60)
  sum(abs(crc_imp_mat - log10(1.01)) < 0.00001)/(300 * 60)
  
  imputed_mat <- sim_tab_zi
  imputed_mat[which(condition == 1), ] = as.matrix(crc_imp_mat)
  imputed_mat[which(condition == 0), ] = as.matrix(control_imp_mat)
  
  sum(((sim_tab - imputed_mat)^2))
  
  #explore the tsne plot, PCA and heatmap
  pdf(file = paste0("plots/", imp_method, "_dm", delta_mu,"_avgmu", avg_mu,"_noise", noise, "_test_thresh", test_thresh,".pdf", collapse = ""))
  # sim_tsne_out <- tsne(imputed_mat)
  # plot(sim_tsne_out, col = as.factor(condition))
  
  sim_rtsne_out <- Rtsne(imputed_mat)
  plot(sim_rtsne_out$Y, col = as.factor(condition))
  
  sim_pca_out <- prcomp(imputed_mat, center = TRUE)
  plot(sim_pca_out$x[,1:2], col = as.factor(condition))
  
  heatmap(as.matrix(imputed_mat), scale="row")
  
  distmat <- as.matrix(dist(as.matrix(imputed_mat)))
  diag(distmat) = max(distmat)
  heatmap(distmat)
  dev.off()
  
  #now check the wilcoxon test results
  less_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "less", p.adjust.method = "bonf")$p.value)
  greater_pval_imp <- apply(imputed_mat, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "greater", p.adjust.method = "bonf")$p.value)
  less_identified_imp <- which(less_pval_imp <= test_thresh)
  greater_identified_imp <- which(greater_pval_imp <= test_thresh)
  sum(less_identified_imp %in% less_identified_full)
  sum(greater_identified_imp %in% greater_identified_full)
  
  length(less_identified_imp)
  length(greater_identified_imp)
  
  length(less_identified_full)
  length(greater_identified_full)
  
  #make a .csv file that contains all the comparisons between t chosen taxa
  result <- matrix(NA, nrow = 3, ncol = 4)
  result[1,1] = length(less_identified_full)
  result[1,2] = length(greater_identified_full)
  result[2,1] = length(less_identified_zi)
  result[2,2] = length(greater_identified_zi)
  result[3,1] = length(less_identified_imp)
  result[3,2] = length(greater_identified_imp)
  result[2,3] = sum(less_identified_zi %in% less_identified_full)
  result[2,4] = sum(greater_identified_zi %in% greater_identified_full)
  result[3,3] = sum(less_identified_imp %in% less_identified_full)
  result[3,4] = sum(greater_identified_imp %in% greater_identified_full)
  # result[2,5] = sum(less_identified_zi %in% less_identified_imp)
  # result[3,5] = sum(greater_identified_zi %in% greater_identified_imp)
  
  colnames(result) <- c("# of taxa in CRC < control", "# of taxa in CRC > control", "overlap with full data in <", "overlap with full data in >")
  rownames(result) <- c("Full", "ZI", "Imputed")
  
  # write.csv(result, file = paste0("results/summary_dm", delta_mu,"_avgmu", avg_mu,"_noise", noise, "_test_thresh", test_thresh, ".csv", collapse = ""))
  return(list("precision" = (result[3,3] + result[3,4])/(result[3,1] + result[3,2]), "recall" = (result[3,3] + result[3,4])/(result[1,1] + result[1,2]) ))
}


gamma_norm_mix <- function(y, X){
  loglik <- function(p, alpha, beta, cov_par, var1, X, y){
    n = length(y)
    lkval <- 0
    fgam <- dgamma(y, shape = alpha, rate = beta)
    for(i in 1:n){
      if(!is.vector(X)){
        lkval <- lkval + log10( p*fgam[i] + (1-p)*dnorm(y[i], mean = X[i,] %*% cov_par, sd = sqrt(var1)))
      }else{
        lkval <- lkval + log10( p*fgam[i] + (1-p)*dnorm(y[i], mean = X[i] * cov_par, sd = sqrt(var1)))
      }
      
    }
    return(lkval)
  }
  n = length(y)
  alpha_init <- 1
  beta_init <- 10
  p_init <- 0.5
  cov_par_init <- solve(t(X) %*% X) %*% t(X) %*% y
  var_init <- t(y - X %*% cov_par_init) %*% (y - X %*% cov_par_init) / n
  
  alpha_t <- alpha_init
  beta_t <- beta_init
  cov_par_t <- cov_par_init
  var_t <- var_init
  p_t <- p_init
  
  #update gamam param
  #Wei's Method
  ### root-finding equation
  fn = function(alpha, target){
    log(alpha) - digamma(alpha) - target
  }
  update_gmm_pars = function(x, wt){
    if(max(wt) > 0.00001){
      tp_s = sum(wt)
      tp_t = sum(wt * x)
      tp_u = sum(wt * log(x))
      tp_v = -tp_u / tp_s - log(tp_s / tp_t)
      if (tp_v <= 0){
        alpha = 20
      }else{
        alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
        if (alpha0 >= 20){alpha = 20
        }else{
          alpha = uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v,
                          extendInt = "yes")$root
        }
      }
      ## need to solve log(x) - digamma(x) = tp_v
      ## We use this approximation to compute the initial value
      beta = tp_s / tp_t * alpha
    }else{
      alpha = 0.001
      beta = 1000
    }
    
    return(c(alpha, beta))
  }
  #convergence criteria
  flag = TRUE
  maxitr = 300
  itr = 0
  while(flag){
    #E_step
    mean_t <- X %*% cov_par_t
    n = length(y)
    a_hat_t <- rep(0, n)
    dg_t <- dgamma(y, shape = alpha_t, rate = beta_t)
    for(i in 1:n){
      if(dg_t[i] == 0){
        a_hat_t[i] = 0
      }else{
        a_hat_t[i] <- p_t * dg_t[i]/(p_t * dg_t[i]+ (1-p_t)*dnorm(y[i], mean = mean_t[i], sd = sqrt(var_t)))
      }
    }
    #maximization
    #fit p
    p_t1 <- sum(a_hat_t)/n
    X_tilta <- sqrt(1-a_hat_t) * X
    y_tilta <- sqrt(1-a_hat_t) * y
    #fit normal
    out <- tryCatch(
      {
        # Just to highlight: if you want to use more than one
        # R expression in the "try" part then you'll have to
        # use curly brackets.
        # 'tryCatch()' will return the last evaluated expression
        # in case the "try" part was completed successfully
        cov_par_t1 <- solve(t(X_tilta) %*% X_tilta) %*% t(X_tilta) %*% y_tilta
        TRUE
      },
      error=function(cond) {
        FALSE
      }
    )
    if(!out){
      return(list("d" = rep(1, n)))
    }
    cov_par_t1 <- solve(t(X_tilta) %*% X_tilta) %*% t(X_tilta) %*% y_tilta
    var_t1 <- sum((1 - a_hat_t) * (y - X %*% cov_par_t)^2) / sum(1-a_hat_t)
    #fit gamma
    par_gamm <- update_gmm_pars(x = y, wt = a_hat_t)
    alpha_t1 = par_gamm[1]
    beta_t1 <- par_gamm[2]
    loglik1 <- loglik(p = p_t, alpha = alpha_t, beta = beta_t, cov_par = cov_par_t, var1 = var_t, X = X, y = y)
    loglik2 <- loglik(p = p_t1, alpha = alpha_t1, beta = beta_t1, cov_par = cov_par_t1, var1 = var_t1, X = X, y = y)
    if((abs(loglik1 - loglik2)) < 0.05 || itr > maxitr){
      flag = FALSE
    }else{
      alpha_t <- alpha_t1
      beta_t <- beta_t1
      cov_par_t <- cov_par_t1
      var_t <- var_t1
      p_t <- p_t1
      itr = itr + 1
    }
  }
  #fit a normal curve
  eta_hat <- cov_par_init
  omega_hat <- var_init
  #calculate new likelihood
  norm_log_lik = 0
  for(i in 1:n){
    if(!is.vector(X)){
      norm_log_lik = norm_log_lik + log10(dnorm(y[i], mean = X[i,] %*% eta_hat, sd = sqrt(omega_hat)))
    }else{
      norm_log_lik = norm_log_lik + log10(dnorm(y[i], mean = X[i] * eta_hat, sd = sqrt(omega_hat)))
    }
  }
  Dev = -2 * norm_log_lik - (-2 * loglik2)
  judgement = pchisq(Dev, df = 3, lower.tail = FALSE) < 0.05
  
  if(!judgement || (alpha_t/beta_t > 1)){
    p_t = 0
    a_hat_t <- rep(0,n)
  }
  return(list("p" = p_t, "alpha" = alpha_t, "beta" = beta_t, "cov_par" = cov_par_t, "var" = var_t, "d" = a_hat_t, "eta" = eta_hat, "omega" = omega_hat, "Deviance" = Dev))
}


ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    lapply(colnames(data_comp), FUN = function(x){
      
    })
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula) 
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
}



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

set.seed(1)
library(sparseDOSSA)
metadata <- matrix(rbinom(n = 100, size = 1, prob = 0.5), nrow = 1, ncol = 100)
simulated_data <- sparseDOSSA(number_features = 150, number_samples = 100,
                              percent_spiked = 0.3, UserMetadata = metadata)
str(simulated_data)


matrix_data <- simulated_data[["OTU_norm"]][[1]]  
df_data1 <- as.data.frame(matrix_data, stringsAsFactors = FALSE)
colnames(df_data1) <- as.character(df_data1[1, ])
rownames(df_data1) <- as.character(df_data1[, 1])
df_data1 <- df_data1[-1, -1]
df_data2 <- rbind(df_data1[1, ], df_data1[302:451, ])
df_data2 <-t(df_data2)


condition <- df_data2[, "Metadata1"]
condition <- as.numeric(condition)

sim_tab_zi_CRC = df_data2[condition == 1,]
sim_tab_zi_CRC <- sim_tab_zi_CRC[, -1]
sim_tab_zi_control = df_data2[condition == 0,]
sim_tab_zi_control <- sim_tab_zi_control[, -1]
str(sim_tab_zi_control)
sim_tab_zi_control <- matrix(as.numeric(sim_tab_zi_control), nrow = 52, ncol = 150)
sim_tab_zi_CRC <- matrix(as.numeric(sim_tab_zi_CRC), nrow = 48, ncol = 150)
imputed_count_mat_list_control <- mbImpute(otu_tab = sim_tab_zi_control)
str(imputed_count_mat_list_control)
imputed_count_mat_list_CRC <- mbImpute(otu_tab = sim_tab_zi_CRC)
CON_imp_mat <- imputed_count_mat_list_control$imp_count_mat_lognorm
CRC_imp_mat <- imputed_count_mat_list_CRC$imp_count_mat_lognorm


imputed_mat <- df_data2[, -1]
imputed_mat <- matrix(as.numeric(imputed_mat), nrow = 100, ncol = 150)
imputed_mat[which(condition == 1), ] = as.matrix(CRC_imp_mat)
imputed_mat[which(condition == 0), ] = as.matrix(CON_imp_mat)

sim_tab <- df_data2[, -1]
sim_tab <- matrix(as.numeric(sim_tab), nrow = 100, ncol = 150)
distmat <- as.matrix(dist(as.matrix(sim_tab)))
diag(distmat) = mean(distmat)
heatmap.2(distmat, col = heatcol, dendrogram = 'none', scale="row", Colv = FALSE, Rowv = FALSE, trace = 'none')

#use t-test to identify truly DE genes in full data
less_pval_full <- apply(sim_tab, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_full <- apply(sim_tab, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
less_identified_full <- which( less_pval_full < 0.05)
greater_identified_full <- which( greater_pval_full < 0.05)

less_pval_zi <- apply(sim_tab, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_zi <- apply(sim_tab, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
# less_identified_zi <- which( p.adjust(less_pval_zi, method = "fdr") < 0.1)
# greater_identified_zi <- which(  p.adjust(greater_pval_zi, method = "fdr") < 0.1)
less_identified_zi <- which( less_pval_zi  < 0.1)
greater_identified_zi <- which(  greater_pval_zi < 0.1)

truth_DA <- c(56,23,19,9,119,127,1,26,72,78,77,20,113,60,114,149,83,59,44,35,133,128,121,42,132,13,142,43,5,76,36,6,46,11,58,123,147,110,40,54,116,39,21,129,57)
length(truth_DA)

Wilcox_precision <- (sum(less_identified_zi %in% truth_DA) + sum(greater_identified_zi %in% truth_DA)) / (length(less_identified_zi) + length(greater_identified_zi))
Wilcox_recall <- (sum(less_identified_zi %in% truth_DA) + sum(greater_identified_zi %in% truth_DA)) / length(truth_DA)
Wilcox_F1_score <- 2 * (Wilcox_precision * Wilcox_recall) / (Wilcox_precision + Wilcox_recall)

Wilcox_metrics <- c(Wilcox_precision, Wilcox_recall, Wilcox_F1_score)

less_pval_mbImpute <- apply(imputed_mat, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_mbImpute <- apply(imputed_mat, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
less_identified_mbImpute <- which(less_pval_mbImpute < 0.1)
greater_identified_mbImpute <- which(greater_pval_mbImpute < 0.1)

Wilcox_mbImpute_precision <- (sum(less_identified_mbImpute %in% truth_DA) +
                                sum(greater_identified_mbImpute %in% truth_DA)) / (length(less_identified_mbImpute) + length(greater_identified_mbImpute))
Wilcox_mbImpute_recall <- (sum(less_identified_mbImpute %in% truth_DA) +
                             sum(greater_identified_mbImpute %in% truth_DA)) / length(truth_DA)
Wilcox_mbImpute_F1_score <- 2 * (Wilcox_mbImpute_precision * Wilcox_mbImpute_recall) / (Wilcox_mbImpute_precision + Wilcox_mbImpute_recall)
Wilcox_mbImpute_metrics <- c(Wilcox_mbImpute_precision, Wilcox_mbImpute_recall, Wilcox_mbImpute_F1_score)


##TphPMF
sim_tab_BHPMF <- sim_tab
sim_tab_BHPMF[sim_tab_BHPMF == 0] <- NA
sim_tab_BHPMF_transposed <- t(sim_tab_BHPMF)


trait.info<- sim_tab_BHPMF_transposed
head(trait.info)

library(usethis)
library(devtools)
library(Matrix)

library(BHPMF)

D <- matrix(0, ncol = 150, nrow = 150)
close1 <- 1: 50
close2 <- 51:100
close3 <- 101:150

for(i in 1:3){
  for(j in 1:3){
    if(i == j){
      tastr = eval(parse(text = paste0("close", i, collapse = "")))
      D[tastr, tastr] = 2
    }else{
      tastr1 = eval(parse(text = paste0("close", i, collapse = "")))
      tastr2 = eval(parse(text = paste0("close", j, collapse = "")))
      D[tastr1, tastr2] = 16
    }
  }
}


diag(D) = 0
dist_obj <- as.dist(D)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 2)
genus_clusters <- cutree(hc, h = 10)
family_clusters <- cutree(hc, h = 15)
microbe_data <- data.frame(
  plant_id = 1:150,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)
microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
hierarchy.info<-microbe_data



head(hierarchy.info)




setwd("/Users/hanxinyu/Desktop/bbhpmf_sim4/")
tmp.dir<-dirname("/Users/hanxinyu/Desktop/bbhpmf_sim4/tmp/")

GapFilling(trait.info,hierarchy.info,used.num.hierarchy.levels=0,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)

data <- read.table("/Users/hanxinyu/Desktop/bbhpmf_sim4/mean_gap_filled.txt", header = TRUE, sep = "\t")
bhmpf_pre <- as.matrix(data)
bhmpf_pre_data<-t(bhmpf_pre)


less_pval_BHPMF <- apply(bhmpf_pre_data, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_BHPMF <- apply(bhmpf_pre_data, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
less_identified_BHPMF <- which(less_pval_BHPMF < 0.1)
greater_identified_BHPMF <- which(greater_pval_BHPMF < 0.1)


Wilcox_BHPMF_precision <- (sum(less_identified_BHPMF %in% truth_DA) +
                             sum(greater_identified_BHPMF %in% truth_DA)) / (length(less_identified_BHPMF) + length(greater_identified_BHPMF))
Wilcox_BHPMF_recall <- (sum(less_identified_BHPMF %in% truth_DA) +
                          sum(greater_identified_BHPMF %in% truth_DA)) / length(truth_DA)
Wilcox_BHPMF_F1_score <- 2 * (Wilcox_BHPMF_precision * Wilcox_BHPMF_recall) / (Wilcox_BHPMF_precision + Wilcox_BHPMF_recall)
Wilcox_BHPMF_metrics <- c(Wilcox_BHPMF_precision, Wilcox_BHPMF_recall, Wilcox_BHPMF_F1_score)

Wilcox_BHPMF_metrics


############ ANCOM ##############

Vardat <- cbind(rnorm(100,0,1), condition)
colnames(Vardat) <- c("rand_cov", "condition")
Vardat <- cbind(1:100, Vardat)
colnames(Vardat)[1] <- "Sample.ID"
OTUdat <- 10^sim_tab - 1.01
taxanames <- unlist( lapply(1:dim(sim_tab)[2], FUN = function(x){
  paste("taxa", x, sep = "")
}))
colnames(OTUdat) <- taxanames
OTUdat <- cbind(1:100, OTUdat)
colnames(OTUdat)[1] = "Sample.ID"

Vardat <- as.data.frame(Vardat)
OTUdat <- as.data.frame(OTUdat)

comparison_test=ANCOM.main(OTUdat,
                           Vardat,
                           adjusted=F,
                           repeated=F,
                           main.var="condition",
                           adj.formula=NULL,
                           repeat.var=NULL,
                           multcorr=2,
                           sig=0.9,
                           prev.cut=0.90,
                           longitudinal = F)

ANCOM_detected <- comparison_test$W.taxa$otu.names[comparison_test$W.taxa$detected_0.6]
ancom_detected <- unlist(lapply(strsplit(as.character(ANCOM_detected), split = ""), function(x){
  paste0(x[5:length(x)], collapse = "")
}))

ancom_detected <- as.numeric(ancom_detected)
length(ancom_detected) 
ancom_precision <- sum((ancom_detected %in% truth_DA))/ (length(ancom_detected))
ancom_recall <- sum((ancom_detected %in% truth_DA))/ (length(truth_DA))
ancom_F1_score <- 2 * (ancom_precision * ancom_recall) / (ancom_precision + ancom_recall)

ancom_metrics <- c(ancom_precision, ancom_recall, ancom_F1_score)
ancom_metrics


#mbimpute
OTUdat1 <- 10^imputed_mat - 1.01
taxanames <- unlist( lapply(1:dim(imputed_mat)[2], FUN = function(x){
  paste("taxa", x, sep = "")
}))
colnames(OTUdat1) <- taxanames
OTUdat1 <- cbind(1:100, OTUdat1)
colnames(OTUdat1)[1] = "Sample.ID"

OTUdat1 <- as.data.frame(OTUdat1)

comparison_test1=ANCOM.main(OTUdat1,
                            Vardat,
                            adjusted=F,
                            repeated=F,
                            main.var="condition",
                            adj.formula=NULL,
                            repeat.var=NULL,
                            multcorr=2,
                            sig=0.9,
                            prev.cut=0.90,
                            longitudinal = F)

ANCOM_detected1 <- comparison_test1$W.taxa$otu.names[comparison_test1$W.taxa$detected_0.6]
ancom_detected1 <- unlist(lapply(strsplit(as.character(ANCOM_detected1), split = ""), function(x){
  paste0(x[5:length(x)], collapse = "")
}))

ancom_detected1 <- as.numeric(ancom_detected1)
length(ancom_detected1) 
ancom_mbimpute_precision <- sum((ancom_detected1 %in% truth_DA))/ (length(ancom_detected1))
ancom_mbimpute_recall <- sum((ancom_detected1 %in% truth_DA))/ (length(truth_DA))
ancom_mbimpute_F1_score <- 2 * (ancom_mbimpute_precision * ancom_mbimpute_recall) / (ancom_mbimpute_precision + ancom_mbimpute_recall)
ancom_mbimpute_metrics <- c(ancom_mbimpute_precision, ancom_mbimpute_recall, ancom_mbimpute_F1_score)

ancom_mbimpute_metrics

#TphPMF
OTUdat2 <- 10^bhmpf_pre_data -1.01
taxanames <- unlist( lapply(1:dim(bhmpf_pre_data)[2], FUN = function(x){
  paste("taxa", x, sep = "")
}))
colnames(OTUdat2) <- taxanames
OTUdat2 <- cbind(1:100, OTUdat2)
colnames(OTUdat2)[1] = "Sample.ID"

OTUdat2 <- as.data.frame(OTUdat2)

comparison_test2=ANCOM.main(OTUdat2,
                            Vardat,
                            adjusted=F,
                            repeated=F,
                            main.var="condition",
                            adj.formula=NULL,
                            repeat.var=NULL,
                            multcorr=2,
                            sig=0.9,
                            prev.cut=0.90,
                            longitudinal = F)

ANCOM_detected2 <- comparison_test2$W.taxa$otu.names[comparison_test2$W.taxa$detected_0.6]
ancom_detected2 <- unlist(lapply(strsplit(as.character(ANCOM_detected2), split = ""), function(x){
  paste0(x[5:length(x)], collapse = "")
}))

ancom_detected2 <- as.numeric(ancom_detected2)
length(ancom_detected2) 


ancom_precison_recall_BHPMF <- c( sum((ancom_detected2 %in% truth_DA))/ (length(ancom_detected2)),
                                  sum((ancom_detected2 %in% truth_DA))/ (length(truth_DA)) )



print(ancom_precison_recall_BHPMF)
ancom_BHPMF_precision <- sum((ancom_detected2 %in% truth_DA))/ (length(ancom_detected2))
ancom_BHPMF_recall <- sum((ancom_detected2 %in% truth_DA))/ (length(truth_DA))
ancom_BHPMF_F1_score <- 2 * (ancom_BHPMF_precision * ancom_BHPMF_recall) / (ancom_BHPMF_precision + ancom_BHPMF_recall)
ancom_BHPMF_metrics <- c(ancom_BHPMF_precision, ancom_BHPMF_recall, ancom_BHPMF_F1_score)

ancom_BHPMF_metrics




library(metagenomeSeq)

Vardat <- cbind(rnorm(100,0,1), condition)
colnames(Vardat) <- c("rand_cov", "condition")

Vardat <- as.data.frame(Vardat)
sim_tab1 <- as.matrix(apply(t(sim_tab), 2, as.numeric))
sim_tab1 <- as.data.frame(sim_tab1)
colnames(sim_tab1) <- 1:100

mr <- newMRexperiment(counts = sim_tab1, phenoData = AnnotatedDataFrame(Vardat))
mr_norm <- cumNorm(mr)
conditions <- factor(Vardat$condition)
fit <- fitZig(obj = mr_norm, mod = model.matrix(~ conditions),useCSSoffset = FALSE,
              zeroInflation = FALSE,
              distribution = "gaussian")

summary(fit)

library(limma)
fit.eb <- eBayes(fit@fit)
pvalues <- fit.eb$p.value
p.values <- pvalues[, "conditions1"]
adjusted_pvalues <- p.adjust(p.values, method = "BH")
significant_indices <- which(adjusted_pvalues < 0.3)
significant_indices_unnamed <- unname(significant_indices)

metagenomeSeq_precision <- (sum(significant_indices_unnamed %in% truth_DA) ) / (length(significant_indices_unnamed) )
metagenomeSeq_recall <- (sum(significant_indices_unnamed %in% truth_DA) ) / (length(truth_DA))
metagenomeSeq_F1_score <- 2 * (metagenomeSeq_precision * metagenomeSeq_recall) / (metagenomeSeq_precision + metagenomeSeq_recall)

metagenomeSeq_metrics <- c(metagenomeSeq_precision, metagenomeSeq_recall, metagenomeSeq_F1_score)
metagenomeSeq_metrics



#mbimpute
sim_tab2 <- as.matrix(apply(t(imputed_mat), 2, as.numeric))
sim_tab2 <- as.data.frame(sim_tab2)
colnames(sim_tab2) <- 1:100

mr2 <- newMRexperiment(counts = sim_tab2, phenoData = AnnotatedDataFrame(Vardat))
mr_norm2 <- cumNorm(mr2)
conditions <- factor(Vardat$condition)
fit2 <- fitZig(obj = mr_norm2, mod = model.matrix(~ conditions),useCSSoffset = FALSE,
               zeroInflation = FALSE,
               distribution = "gaussian")

summary(fit2)
library(limma)

fit.eb2 <- eBayes(fit2@fit)
pvalues2 <- fit.eb2$p.value
p.values2 <- pvalues2[, "conditions1"]
adjusted_pvalues2 <- p.adjust(p.values2, method = "BH")
significant_indices2 <- which(adjusted_pvalues2 < 0.3)
significant_indices_unnamed2 <- unname(significant_indices2)
metagenomeSeq_mbImpute_precision <- (sum(significant_indices_unnamed2 %in% truth_DA) ) / (length(significant_indices_unnamed2) )
metagenomeSeq_mbImpute_recall <- (sum(significant_indices_unnamed2 %in% truth_DA) ) / (length(truth_DA))
metagenomeSeq_mbImpute_F1_score <- 2 * (metagenomeSeq_mbImpute_precision * metagenomeSeq_mbImpute_recall) / (metagenomeSeq_mbImpute_precision + metagenomeSeq_mbImpute_recall)

metagenomeSeq_mbImpute_metrics <- c(metagenomeSeq_mbImpute_precision, metagenomeSeq_mbImpute_recall, metagenomeSeq_mbImpute_F1_score)
metagenomeSeq_mbImpute_metrics


#TphPMF
sim_tab3 <- as.matrix(apply(t(bhmpf_pre_data), 2, as.numeric))
sim_tab3 <- as.data.frame(sim_tab3)
colnames(sim_tab3) <- 1:100
mr3 <- newMRexperiment(counts = sim_tab3, phenoData = AnnotatedDataFrame(Vardat))
mr_norm3 <- cumNorm(mr3)
conditions <- factor(Vardat$condition)
fit3 <- fitZig(obj = mr_norm3, mod = model.matrix(~ conditions),useCSSoffset = FALSE,
               zeroInflation = FALSE,
               distribution = "gaussian")

summary(fit2)

library(limma)
fit.eb3 <- eBayes(fit3@fit)
pvalues3 <- fit.eb3$p.value
p.values3 <- pvalues3[, "conditions1"]
adjusted_pvalues3 <- p.adjust(p.values3, method = "BH")
significant_indices3 <- which(adjusted_pvalues3 < 0.3)
significant_indices_unnamed3 <- unname(significant_indices3)
metagenomeSeq_BHPMF_precision <- (sum(significant_indices_unnamed3 %in% truth_DA) ) / (length(significant_indices_unnamed3) )
metagenomeSeq_BHPMF_recall <- (sum(significant_indices_unnamed3 %in% truth_DA) ) / (length(truth_DA))
metagenomeSeq_BHPMF_F1_score <- 2 * (metagenomeSeq_BHPMF_precision * metagenomeSeq_BHPMF_recall) / (metagenomeSeq_BHPMF_precision + metagenomeSeq_BHPMF_recall)

metagenomeSeq_BHPMF_metrics <- c(metagenomeSeq_BHPMF_precision, metagenomeSeq_BHPMF_recall, metagenomeSeq_BHPMF_F1_score)
metagenomeSeq_BHPMF_metrics


#Omnibus

library(mbzinb)
Vardat <- cbind(rnorm(100,0,1), condition)
colnames(Vardat) <- c("rand_cov", "condition")
mbzinb.dataset <-
  function(count, sample, taxon=NULL) {
    #Return some warnings if data type is unexpected.
    #For now, I don't have anything to do with the taxonomy table so it doesn't matter.
    if (!inherits(count, "matrix")) {
      stop("Count table should have class matrix \n")
    }
    if (!inherits(sample, "data.frame")) {
      stop("Sample information should have class data.frame \n")
    }
    #Check for null names
    if (is.null(colnames(count))) {
      stop("Count matrix must have column names matching sample row names \n")
    }
    if (is.null(rownames(sample))) {
      stop("Sample data frame must have column names \n")
    }
    if (!is.null(taxon) & is.null(rownames(taxon))) {
      stop("Taxonomy table must have row names \n")
    }
    #Reorder samples to match names
    sample.match <- match(colnames(count), rownames(sample))
    if (any(is.na(sample.match))) {
      cat(paste(sum(is.na(sample.match)), "samples in count table trimmed due to missing sample information \n"))
      new.count <- count[, !is.na(sample.match)]
    } else {
      new.count <- count
    }
    new.sample <- sample[colnames(new.count), , drop=FALSE]
    #If taxonomy data supplied, reorder to match names
    if (!is.null(taxon)) {
      taxon.match <- match(rownames(count), rownames(taxon))
      if (any(is.na(taxon.match))) {
        cat(paste(sum(is.na(taxon.match)), "taxa in count table trimmed due to missing taxonomy information \n"))
        new.count <- new.count[!is.na(taxon.match), ]
      } else {
        new.count <- new.count
      }
      new.taxon <- taxon[rownames(new.count), ]
    } else {
      new.taxon <- NULL
    }
    l <- list(count=new.count, sample=new.sample, taxon=new.taxon, filtered=FALSE)
    class(l) <- "mbzinb"
    return(l)
  }



Vardat <- as.data.frame(Vardat)
OTUdat <- 100000000000000000^sim_tab
sim_tabb<-t(OTUdat)
colnames(sim_tabb) <- rownames(Vardat)
mbzinb_data <- mbzinb.dataset(sim_tabb, Vardat)
mbzinb_data_filtered <- filter.dataset(mbzinb_data, min.prev=0.1, min.reads=50, niter=1)
mbzinb_test_result <- mbzinb.test(mbzinb_data_filtered, group = "condition")
p_values <- mbzinb_test_result$results$PValue
p_values_adjusted <- p.adjust(p_values, method = "BH")
print(p_values_adjusted)
significant_indices_Omnibus <- which(p_values_adjusted < 0.05)
significant_indices_Omnibus_unnamed <- unname(significant_indices_Omnibus )
Omnibus_precision <- (sum(significant_indices_Omnibus_unnamed %in% truth_DA) ) / (length(significant_indices_Omnibus_unnamed) )
Omnibus_recall <- (sum(significant_indices_Omnibus_unnamed %in% truth_DA) ) / (length(truth_DA))
Omnibus_F1_score <- 2 * (Omnibus_precision * Omnibus_recall) / (Omnibus_precision + Omnibus_recall)

Omnibus_metrics <- c(Omnibus_precision, Omnibus_recall, Omnibus_F1_score)
Omnibus_metrics

#mbImpute
OTUdat1 <- 100000000000000000^imputed_mat
sim_tabb1<-t(OTUdat1)
colnames(sim_tabb1) <- rownames(Vardat)
mbzinb_data1 <- mbzinb.dataset(sim_tabb1, Vardat)
mbzinb_data_filtered1 <- filter.dataset(mbzinb_data1, min.prev=0.1, min.reads=50, niter=1)
mbzinb_test_result1 <- mbzinb.test(mbzinb_data_filtered1, group = "condition")
p_values1 <- mbzinb_test_result1$results$PValue
p_values_adjusted1 <- p.adjust(p_values1, method = "BH")
print(p_values_adjusted1)
significant_indices_Omnibus1 <- which(p_values_adjusted1 < 0.05)
significant_indices_Omnibus_unnamed1 <- unname(significant_indices_Omnibus1 )
Omnibus_mbImpute_precision <- (sum(significant_indices_Omnibus_unnamed1 %in% truth_DA) ) / (length(significant_indices_Omnibus_unnamed1) )
Omnibus_mbImpute_recall <- (sum(significant_indices_Omnibus_unnamed1 %in% truth_DA) ) / (length(truth_DA))
Omnibus_mbImpute_F1_score <- 2 * (Omnibus_mbImpute_precision * Omnibus_mbImpute_recall) / (Omnibus_mbImpute_precision + Omnibus_mbImpute_recall)

Omnibus_mbImpute_metrics <- c(Omnibus_mbImpute_precision, Omnibus_mbImpute_recall, Omnibus_mbImpute_F1_score)
Omnibus_mbImpute_metrics

##TphPMF

OTUdat2 <- 100000000000000000000000000000000000000000000000000000000^bhmpf_pre_data+10
sim_tabb2<-t(OTUdat2)
colnames(sim_tabb2) <- rownames(Vardat)
mbzinb_data2 <- mbzinb.dataset(sim_tabb2, Vardat)
mbzinb_data_filtered2 <- filter.dataset(mbzinb_data2, min.prev=0.1, min.reads=0.000005, niter=1)
mbzinb_test_result2 <- mbzinb.test(mbzinb_data_filtered2, group = "condition")
p_values2 <- mbzinb_test_result2$results$PValue
print(p_values2)
p_values_adjusted2 <- p.adjust(p_values2, method = "BH")
print(p_values_adjusted2)
significant_indices_Omnibus2 <- which(p_values_adjusted2 < 0.05)
print(significant_indices_Omnibus2)
significant_indices_Omnibus_unnamed2 <- unname(significant_indices_Omnibus2 )

Omnibus_BHPMF_precision <- (sum(significant_indices_Omnibus_unnamed2 %in% truth_DA) ) / (length(significant_indices_Omnibus_unnamed2) )
Omnibus_BHPMF_recall <- (sum(significant_indices_Omnibus_unnamed2 %in% truth_DA) ) / (length(truth_DA))
Omnibus_BHPMF_F1_score <- 2 * (Omnibus_BHPMF_precision * Omnibus_BHPMF_recall) / (Omnibus_BHPMF_precision + Omnibus_BHPMF_recall)

Omnibus_BHPMF_metrics <- c(Omnibus_BHPMF_precision, Omnibus_BHPMF_recall, Omnibus_BHPMF_F1_score)



########################DESeq2-phyloseq#################3
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("joey711/phyloseq")
library(phyloseq)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(DESeq2)

Vardat <- cbind(rnorm(100,0,1), condition)
colnames(Vardat) <- c("rand_cov", "condition")
Vardat <- as.data.frame(Vardat)
sim_tabbb<-t(sim_tab)
OTUdat <- floor(100000000000000^sim_tabbb +1)
physeq2 <- phyloseq(otu_table(OTUdat, taxa_are_rows = TRUE), 
                    sample_data(Vardat))
sample_data(physeq2)$condition <- factor(sample_data(physeq2)$condition)
Deseq2_obj <- phyloseq_to_deseq2(physeq2, ~ condition)
results <- DESeq(Deseq2_obj, test="Wald", fitType="parametric")
dds_results <- results(results)
yyresults <- dds_results$padj

print(yyresults)

significant_indices_DESeq2 <- which(yyresults < 0.05)
print(significant_indices_DESeq2)
significant_indices_DESeq2_unnamed2 <- unname(significant_indices_DESeq2 )

DESeq2_precision <- (sum(significant_indices_DESeq2_unnamed2 %in% truth_DA) ) / (length(significant_indices_DESeq2_unnamed2) )
DESeq2_recall <- (sum(significant_indices_DESeq2_unnamed2 %in% truth_DA) ) / (length(truth_DA))
DESeq2_F1_score <- 2 * (DESeq2_precision * DESeq2_recall) / (DESeq2_precision + DESeq2_recall)
DESeq2_metrics <- c(DESeq2_precision, DESeq2_recall, DESeq2_F1_score)





##mbImpute
sim_tabbb1<-t(imputed_mat)
physeq2 <- phyloseq(otu_table(sim_tabbb1, taxa_are_rows = TRUE), 
                    sample_data(Vardat))
sample_data(physeq2)$condition <- factor(sample_data(physeq2)$condition)
Deseq2_obj <- phyloseq_to_deseq2(physeq2, ~ condition)
results <- DESeq(Deseq2_obj, test="Wald", fitType="parametric")
dds_results <- results(results)
yyresults <- dds_results$padj

print(yyresults)
significant_indices_DESeq2 <- which(yyresults < 0.05)
print(significant_indices_DESeq2)
significant_indices_DESeq2_unnamed2 <- unname(significant_indices_DESeq2 )
DESeq2_mbimpute_precision <- (sum(significant_indices_DESeq2_unnamed2 %in% truth_DA) ) / (length(significant_indices_DESeq2_unnamed2) )
DESeq2_mbimpute_recall <- (sum(significant_indices_DESeq2_unnamed2 %in% truth_DA) ) / (length(truth_DA))
DESeq2_mbimpute_F1_score <- 2 * (DESeq2_mbimpute_precision * DESeq2_mbimpute_recall) / (DESeq2_mbimpute_precision + DESeq2_mbimpute_recall)

DESeq2_mbimpute_metrics <- c(DESeq2_mbimpute_precision, DESeq2_mbimpute_recall, DESeq2_mbimpute_F1_score)



##TphPMF

sim_tabbb1<-t(bhmpf_pre_data)
OTUdat <- floor(10000000000000000000000000000000000000000000000000000000000000000^sim_tabbb1 +1)

colnames(OTUdat) <- rownames(Vardat)
otu_table_object <- otu_table(OTUdat, taxa_are_rows = TRUE)
otu_sample_names <- colnames(otu_table_object)
sample_data_sample_names <- rownames(Vardat)
setequal(otu_sample_names, sample_data_sample_names)

otu_sample_names <- colnames(otu_table_object)
sample_data_sample_names <- rownames(Vardat)

sample_data(physeq2)$condition <- factor(sample_data(physeq2)$condition)
Deseq2_obj <- phyloseq_to_deseq2(physeq2, ~ condition)
results <- DESeq(Deseq2_obj, test="Wald", fitType="parametric")
dds_results <- results(results)
yyresults <- dds_results$padj

print(yyresults)
significant_indices_DESeq2 <- which(yyresults < 0.05)
print(significant_indices_DESeq2)
significant_indices_DESeq2_unnamed2 <- unname(significant_indices_DESeq2 )

DESeq2_BHPMF_precision <- (sum(significant_indices_DESeq2_unnamed2 %in% truth_DA) ) / (length(significant_indices_DESeq2_unnamed2) )
DESeq2_BHPMF_recall <- (sum(significant_indices_DESeq2_unnamed2 %in% truth_DA) ) / (length(truth_DA))
DESeq2_BHPMF_F1_score <- 2 * (DESeq2_BHPMF_precision * DESeq2_BHPMF_recall) / (DESeq2_BHPMF_precision + DESeq2_BHPMF_recall)

DESeq2_BHPMF_metrics <- c(DESeq2_BHPMF_precision, DESeq2_BHPMF_recall, DESeq2_BHPMF_F1_score)


library(ggplot2)
values <- c(
  Wilcox_precision, Wilcox_recall, Wilcox_F1_score,
  Wilcox_mbImpute_precision, Wilcox_mbImpute_recall, Wilcox_mbImpute_F1_score,
  Wilcox_BHPMF_precision, Wilcox_BHPMF_recall, Wilcox_BHPMF_F1_score,
  ancom_precision, ancom_recall, ancom_F1_score,
  ancom_mbimpute_precision, ancom_mbimpute_recall, ancom_mbimpute_F1_score,
  ancom_BHPMF_precision, ancom_BHPMF_recall, ancom_BHPMF_F1_score,
  metagenomeSeq_precision, metagenomeSeq_recall, metagenomeSeq_F1_score,
  metagenomeSeq_mbImpute_precision, metagenomeSeq_mbImpute_recall, metagenomeSeq_mbImpute_F1_score,
  metagenomeSeq_BHPMF_precision, metagenomeSeq_BHPMF_recall, metagenomeSeq_BHPMF_F1_score,
  DESeq2_precision, DESeq2_recall, DESeq2_F1_score,
  DESeq2_mbimpute_precision, DESeq2_mbimpute_recall, DESeq2_mbimpute_F1_score,
  DESeq2_BHPMF_precision, DESeq2_BHPMF_recall, DESeq2_BHPMF_F1_score,
  Omnibus_precision, Omnibus_recall, Omnibus_F1_score,
  Omnibus_mbImpute_precision, Omnibus_mbImpute_recall, Omnibus_mbImpute_F1_score,
  Omnibus_BHPMF_precision, Omnibus_BHPMF_recall, Omnibus_BHPMF_F1_score)
print(values)


results_df <- data.frame(
  Method = rep(c("Wilcoxon", "ANCOM", "metagenomeSeq", "DESeq2-phyloseq", "Omnibus"), each = 9),
  Metric = rep(c("Precision", "Recall", "F1_Score"), times = 5 * 3),
  Value = values,
  Imputation = rep(c("DA method", "mbImpute + DA method", "TphPMF + DA method"), times = 5)
)

ggplot(results_df, aes(x = Metric, y = Value, fill = Imputation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~ Method, scales = "free_y") +
  scale_fill_manual(values = c("DA method" = "blue", "mbImpute + DA method" = "red", "TphPMF + DA method" = "green")) +
  labs(x = "", y = "Scores") +  
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = NULL))  
results_df <- data.frame(
  Method = rep(c("Wilcoxon", "ANCOM", "metagenomeSeq", "DESeq2-phyloseq", "Omnibus"), each = 9),
  Metric = rep(c("Precision", "Recall", "F1_Score"), times = 5*3),
  Value = values,
  Imputation = rep(c("DA method", "mbImpute + DA method", "TphPMF + DA method"), times = 15)
)




ggplot(results_df, aes(x = interaction(Metric, Imputation), y = Value, fill = Imputation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.1)) +
  facet_wrap(~ Method, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "") +  
  ylim(0, 1.00) + 
  theme_minimal() +
  theme(
    legend.position = "bottom",
  )


