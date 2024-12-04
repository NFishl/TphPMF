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


##mbDenoise
library(mbDenoise)
sim_tab_zi_matrix <- as.matrix(sim_tab)
result <- ZIPPCApn(sim_tab_zi_matrix , family = "negative.binomial", n.factors = 2, rank = TRUE)
mbDenoise_imputed <- result$muz
less_pval_mbDenoise <- apply(mbDenoise_imputed, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
greater_pval_mbDenoise <- apply(mbDenoise_imputed, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
less_identified_mbDenoise <- which(less_pval_mbDenoise < 0.1)
greater_identified_mbDenoise <- which(greater_pval_mbDenoise < 0.1)
Wilcox_mbDenoise_precision <- (sum(less_identified_mbDenoise %in% truth_DA) +
sum(greater_identified_mbDenoise %in% truth_DA)) / (length(less_identified_mbDenoise) + length(greater_identified_mbDenoise))
Wilcox_mbDenoise_recall <- (sum(less_identified_mbDenoise %in% truth_DA) +
sum(greater_identified_mbDenoise %in% truth_DA)) / length(truth_DA)
Wilcox_mbDenoise_F1_score <- 2 * (Wilcox_mbDenoise_precision * Wilcox_mbDenoise_recall) / (Wilcox_mbDenoise_precision + Wilcox_mbDenoise_recall)

Wilcox_mbDenoise_metrics <- c(Wilcox_mbDenoise_precision, Wilcox_mbDenoise_recall, Wilcox_mbDenoise_F1_score)


############ ANCOM-BC2 ##############
library(devtools)
library(ANCOMBC)
Vardat <- cbind(rnorm(100, 0, 1), condition)
colnames(Vardat) <- c("rand_cov", "condition")
Vardat <- cbind(1:100, Vardat)
colnames(Vardat)[1] <- "Sample.ID"
OTUdat <- 10^sim_tab 
taxanames <- unlist(lapply(1:dim(sim_tab)[2], FUN = function(x) {
paste("taxa", x, sep = "")
}))
colnames(OTUdat) <- taxanames
OTUdat <- cbind(1:100, OTUdat)
colnames(OTUdat)[1] = "Sample.ID"
Vardat <- as.data.frame(Vardat)
OTUdat <- as.data.frame(OTUdat)
library(phyloseq)
OTU_table <- otu_table(OTUdat, taxa_are_rows = FALSE)
sample_data <- sample_data(Vardat)
physeq <- phyloseq(OTU_table, sample_data)
OTUdat <- as(otu_table(physeq), "matrix")
Vardat <- as(sample_data(physeq), "data.frame")
comparison_test <- ancombc2(
data = OTUdat,
taxa_are_rows = FALSE,
meta_data = Vardat,
fix_formula = "condition",
group = "condition",
p_adj_method = "BH",
prv_cut = 0.1,
lib_cut = 0,
struc_zero = FALSE,
neg_lb = FALSE,
pseudo_sens = FALSE,            
alpha = 0.05,
global = TRUE,
trend = FALSE
)
str(comparison_test)
ancom_detected <- comparison_test$res$taxon[comparison_test$res$diff_condition]
ancom_detected <- as.numeric(gsub("taxa", "", ancom_detected))
ancom_precision <- sum((ancom_detected %in% truth_DA)) / length(ancom_detected)
ancom_recall <- sum((ancom_detected %in% truth_DA)) / length(truth_DA)
ancom_F1_score <- 2 * (ancom_precision * ancom_recall) / (ancom_precision + ancom_recall)
ancom_metrics <- c(ancom_precision, ancom_recall, ancom_F1_score)
ancom_metrics


###mbImpute
OTUdat1 <- imputed_mat
OTUdat1_log<-10^imputed_mat
taxanames <- unlist(lapply(1:dim(imputed_mat)[2], FUN = function(x) {
paste("taxa", x, sep = "")
}))
colnames(OTUdat1_log) <- taxanames
OTUdat1_log <- cbind(1:100, OTUdat1_log)
colnames(OTUdat1_log)[1] = "Sample.ID"
OTUdat1_log <- as.data.frame(OTUdat1_log)
Vardat <- as.data.frame(Vardat)
Vardat$condition <- as.factor(Vardat$condition)
comparison_test1 <- ancombc2(
data = OTUdat1_log,
taxa_are_rows = FALSE,
meta_data = Vardat,
fix_formula = "condition",
group = "condition",
p_adj_method = "BH",
prv_cut = 0.1,
lib_cut = 0,
struc_zero = FALSE,
neg_lb = FALSE,
pseudo_sens = FALSE,            
alpha = 0.05,
global = TRUE,
trend = FALSE
)
str(comparison_test1)
ancom_detected1 <- comparison_test1$res$taxon[comparison_test1$res$diff_condition]
ancom_detected1 <- as.numeric(gsub("taxa", "", ancom_detected1))
ancom_precision_mbImpute <- sum((ancom_detected1 %in% truth_DA)) / length(ancom_detected1)
ancom_recall_mbImpute <- sum((ancom_detected1 %in% truth_DA)) / length(truth_DA)
ancom_F1_score_mbImpute <- 2 * (ancom_precision_mbImpute * ancom_recall_mbImpute) / (ancom_precision_mbImpute + ancom_recall_mbImpute)
ancom_metrics_mbImpute <- c(ancom_precision_mbImpute, ancom_recall_mbImpute, ancom_F1_score_mbImpute)
ancom_metrics_mbImpute



###TphPMF
OTUdat2 <- 10^(bhmpf_pre_data)
taxanames <- unlist(lapply(1:dim(bhmpf_pre_data)[2], FUN = function(x) {
paste("taxa", x, sep = "")
}))
colnames(OTUdat2) <- taxanames
OTUdat2 <- cbind(1:100, OTUdat2)
colnames(OTUdat2)[1] = "Sample.ID"
OTUdat2 <- as.data.frame(OTUdat2)
comparison_test2 <- ancombc2(
data = OTUdat2,
taxa_are_rows = FALSE,
meta_data = Vardat,
fix_formula = "condition",
group = "condition",
p_adj_method = "BH",
prv_cut = 0.1,
lib_cut = 0,
struc_zero = FALSE,
neg_lb = FALSE,
pseudo_sens = FALSE,            
alpha = 0.05,
global = TRUE,
trend = FALSE
)
ancom_detected2 <- comparison_test2$res$taxon[comparison_test2$res$diff_condition]
ancom_detected2 <- as.numeric(gsub("taxa", "", ancom_detected2))
ancom_precision_TphPMF <- sum((ancom_detected2 %in% truth_DA)) / length(ancom_detected2)
ancom_recall_TphPMF <- sum((ancom_detected2 %in% truth_DA)) / length(truth_DA)
ancom_F1_score_TphPMF <- 2 * (ancom_precision_TphPMF * ancom_recall_TphPMF) / (ancom_precision_TphPMF + ancom_recall_TphPMF)
ancom_metrics_TphPMF <- c(ancom_precision_TphPMF, ancom_recall_TphPMF, ancom_F1_score_TphPMF)
ancom_metrics_TphPMF


###mbDenoise
OTUdat3<-10^mbDenoise_imputed
taxanames <- unlist(lapply(1:dim(imputed_mat)[2], FUN = function(x) {
paste("taxa", x, sep = "")
}))
OTUdat3_log <- cbind(1:100, OTUdat3)
colnames(OTUdat3_log)[1] = "Sample.ID"
OTUdat3_log <- as.data.frame(OTUdat3_log)
Vardat <- as.data.frame(Vardat)
Vardat$condition <- as.factor(Vardat$condition)
comparison_test3 <- ancombc2(
data = OTUdat3_log,
taxa_are_rows = FALSE,
meta_data = Vardat,
fix_formula = "condition",
group = "condition",
p_adj_method = "BH",
prv_cut = 0.1,
lib_cut = 0,
struc_zero = FALSE,
neg_lb = FALSE,
pseudo_sens = FALSE,           
alpha = 0.05,
global = TRUE,
trend = FALSE
)
str(comparison_test3)
ancom_detected3 <- comparison_test3$res$taxon[comparison_test3$res$diff_condition]
ancom_precision3 <- sum((ancom_detected3 %in% truth_DA)) / length(ancom_detected3)
ancom_recall3 <- sum((ancom_detected3 %in% truth_DA)) / length(truth_DA)
ancom_F1_score3 <- 2 * (ancom_precision3 * ancom_recall3) / (ancom_precision3 + ancom_recall3)
ancom_metrics_mbDenoise <- c(ancom_precision3, ancom_recall3, ancom_F1_score3)
ancom_metrics_mbDenoise


######LOCOM#####
library(LOCOM)
Vardat <- cbind(rnorm(100, 0, 1), condition)  
colnames(Vardat) <- c("rand_cov", "condition")  
Vardat <- cbind(1:100, Vardat)  
colnames(Vardat)[1] <- "Sample.ID"  
OTUdat <- 10^sim_tab  
taxanames <- unlist(lapply(1:dim(sim_tab)[2], FUN = function(x) {
paste("taxa", x, sep = "")
}))  
colnames(OTUdat) <- taxanames  
OTUdat <- cbind(1:100, OTUdat)  
colnames(OTUdat)[1] <- "Sample.ID" 
Vardat <- as.data.frame(Vardat)
OTUdat <- as.data.frame(OTUdat)

Y <- Vardat$condition  
C <- Vardat$rand_cov  
result <- locom(
otu.table = OTUdat[, -1],  
Y = Y,                     
C = C,                   
seed = 1,                 
n.cores = 4,               
n.perm.max = 1000          
)
locom_detected<-result$detected.otu
locom_detected <- as.numeric(gsub("taxa", "", locom_detected))
locom_precision <- sum((locom_detected %in% truth_DA)) / length(locom_detected)
locom_recall <- sum((locom_detected %in% truth_DA)) / length(truth_DA)
locom_F1_score <- 2 * (locom_precision * locom_recall) / (locom_precision + locom_recall)
locom_metrics <- c(locom_precision, locom_recall, locom_F1_score)
locom_metrics

###mbImpute
Vardat <- cbind(rnorm(100, 0, 1), condition) 
colnames(Vardat) <- c("rand_cov", "condition") 
Vardat <- cbind(1:100, Vardat)  
colnames(Vardat)[1] <- "Sample.ID"  
OTUdat1 <- imputed_mat  
taxanames <- unlist(lapply(1:dim(sim_tab)[2], FUN = function(x) {
paste("taxa", x, sep = "")
})) 
colnames(OTUdat1) <- taxanames 
OTUdat1 <- cbind(1:100, OTUdat1)  
colnames(OTUdat1)[1] <- "Sample.ID" 
Vardat <- as.data.frame(Vardat)
OTUdat1 <- as.data.frame(OTUdat1)

Y <- Vardat$condition  
C <- Vardat$rand_cov   
result <- locom(
otu.table = OTUdat1[, -1],  
Y = Y,                    
C = C,                    
seed = 1,                  
n.cores = 4,               
n.perm.max = 1000          
)
locom_detected1<-result$detected.otu
locom_detected1 <- as.numeric(gsub("taxa", "", locom_detected1))
locom_precision_mbImpute <- sum((locom_detected1 %in% truth_DA)) / length(locom_detected1)
locom_recall_mbImpute <- sum((locom_detected1 %in% truth_DA)) / length(truth_DA)
locom_F1_score_mbImpute <- 2 * (locom_precision_mbImpute * locom_recall_mbImpute) / (locom_precision_mbImpute + locom_recall_mbImpute)
locom_metrics_mbImpute <- c(locom_precision_mbImpute, locom_recall_mbImpute, locom_F1_score_mbImpute)
locom_metrics_mbImpute

###Tphpmf
Vardat <- cbind(rnorm(100, 0, 1), condition)  
colnames(Vardat) <- c("rand_cov", "condition")  
Vardat <- cbind(1:100, Vardat)  
colnames(Vardat)[1] <- "Sample.ID"  
OTUdat2 <- bhmpf_pre_data   
taxanames <- unlist(lapply(1:dim(sim_tab)[2], FUN = function(x) {
paste("taxa", x, sep = "")
})) 
colnames(OTUdat2) <- taxanames  
OTUdat2 <- cbind(1:100, OTUdat2)  
colnames(OTUdat2)[1] <- "Sample.ID"  
Vardat <- as.data.frame(Vardat)
OTUdat2 <- as.data.frame(OTUdat2)

Y <- Vardat$condition  
C <- Vardat$rand_cov  
result <- locom(
otu.table = OTUdat2[, -1],  
Y = Y,                    
C = C,                     
seed = 1,                  
n.cores = 4,               
n.perm.max = 1000          
)
print(result$detected.otu)
locom_detected2<-result$detected.otu
locom_detected2 <- as.numeric(gsub("taxa", "", locom_detected2))
locom_precision_TphPMF <- sum((locom_detected2 %in% truth_DA)) / length(locom_detected2)
locom_recall_TphPMF <- sum((locom_detected2 %in% truth_DA)) / length(truth_DA)
locom_F1_score_TphPMF <- 2 * (locom_precision_TphPMF * locom_recall_TphPMF) / (locom_precision_TphPMF + locom_recall_TphPMF)
locom_metrics_TphPMF <- c(locom_precision_TphPMF, locom_recall_TphPMF, locom_F1_score_TphPMF)
locom_metrics_TphPMF

###mbDenoise
Vardat <- cbind(rnorm(100, 0, 1), condition)  
colnames(Vardat) <- c("rand_cov", "condition")  
Vardat <- cbind(1:100, Vardat)  
colnames(Vardat)[1] <- "Sample.ID" 
OTUdat3 <- 10^mbDenoise_imputed
taxanames <- unlist(lapply(1:dim(sim_tab)[2], FUN = function(x) {
paste("taxa", x, sep = "")
}))  
colnames(OTUdat3) <- taxanames  
OTUdat3 <- cbind(1:100, OTUdat3) 
colnames(OTUdat3)[1] <- "Sample.ID" 
Vardat <- as.data.frame(Vardat)
OTUdat3 <- as.data.frame(OTUdat3)
Y <- Vardat$condition  
C <- Vardat$rand_cov   
result1 <- locom(
otu.table = OTUdat3[, -1],  
Y = Y,                     
C = C,                    
seed = 1,                  
n.cores = 4,               
n.perm.max = 1000          
)
print(result1$detected.otu)
locom_detected3<-result1$detected.otu
locom_detected3 <- as.numeric(gsub("taxa", "", locom_detected3))
locom_precision_mbDenoise <- sum((locom_detected3 %in% truth_DA)) / length(locom_detected3)
locom_recall_mbDenoise <- sum((locom_detected3 %in% truth_DA)) / length(truth_DA)
locom_F1_score_mbDenoise <- 2 * (locom_precision_mbDenoise * locom_recall_mbDenoise) / (locom_precision_mbDenoise + locom_recall_mbDenoise)
locom_metrics_mbDenoise <- c(locom_precision_mbDenoise, locom_recall_mbDenoise, locom_F1_score_mbDenoise)
locom_metrics_mbDenoise


######Linda######
library(LinDA)
Vardat <- cbind(rnorm(100, 0, 1), condition)  
colnames(Vardat) <- c("rand_cov", "condition") 
Vardat <- cbind(1:100, Vardat)  
colnames(Vardat)[1] <- "Sample.ID"  
OTUdat <- 10^sim_tab   
taxanames <- unlist(lapply(1:dim(sim_tab)[2], FUN = function(x) {
paste("taxa", x, sep = "")
}))  
colnames(OTUdat) <- taxanames  
OTUdat <- cbind(1:100, OTUdat)  
colnames(OTUdat)[1] <- "Sample.ID" 
Vardat <- as.data.frame(Vardat)
OTUdat <- as.data.frame(OTUdat)
otu.tab <- OTUdat[, -1]
otu.tab<-t(otu.tab)
meta <- Vardat
colnames(meta)[2] <- "rand_cov" 
colnames(meta)[3] <- "condition"  
linda.obj <- linda(
otu.tab = otu.tab,        
meta = meta,             
formula = '~condition+rand_cov',  
alpha = 0.05,            
prev.cut = 0.001,          
lib.cut = 0.01,           
winsor.quan = NULL        
)
significant_otus <- rownames(linda.obj$output[[1]])[which(linda.obj$output[[1]]$reject)]
print(significant_otus)
linda_detected<-significant_otus
linda_detected <- as.numeric(gsub("taxa", "", linda_detected))
linda_precision <- sum((linda_detected %in% truth_DA)) / length(linda_detected)
linda_recall <- sum((linda_detected %in% truth_DA)) / length(truth_DA)
linda_F1_score <- 2 * (linda_precision * linda_recall) / (linda_precision + linda_recall)
linda_metrics <- c(linda_precision, linda_recall, linda_F1_score)
linda_metrics

###mbimpute
Vardat <- cbind(rnorm(100, 0, 1), condition) 
colnames(Vardat) <- c("rand_cov", "condition") 
Vardat <- cbind(1:100, Vardat)  
colnames(Vardat)[1] <- "Sample.ID" 

OTUdat1 <- imputed_mat   
taxanames <- unlist(lapply(1:dim(sim_tab)[2], FUN = function(x) {
paste("taxa", x, sep = "")
}))  
print(taxanames)
colnames(OTUdat1) <- taxanames  
OTUdat1 <- cbind(1:100, OTUdat1)  
colnames(OTUdat1)[1] <- "Sample.ID"  
Vardat <- as.data.frame(Vardat)
OTUdat <- as.data.frame(OTUdat1)
otu.tab <- OTUdat1[, -1]
otu.tab<-t(otu.tab)
meta <- Vardat
colnames(meta)[2] <- "rand_cov"  
colnames(meta)[3] <- "condition"  
linda.obj1 <- linda(
otu.tab = otu.tab,        
meta = meta,             
formula = '~condition+rand_cov', 
alpha = 0.05,             
prev.cut = 0.001,          
lib.cut = 0.01,           
winsor.quan = NULL        
  )
str(linda.obj1)
significant_otus1 <- rownames(linda.obj1$output[[1]])[which(linda.obj1$output[[1]]$reject)]
linda_detected1<-significant_otus1
linda_detected1 <- as.numeric(gsub("taxa", "", linda_detected1))
linda_precision1 <- sum((linda_detected1 %in% truth_DA)) / length(linda_detected1)
linda_recall1 <- sum((linda_detected1 %in% truth_DA)) / length(truth_DA)
linda_F1_score1 <- 2 * (linda_precision1 * linda_recall1) / (linda_precision1 + linda_recall1)
linda_metrics_mbImpute <- c(linda_precision1, linda_recall1, linda_F1_score1)
linda_metrics_mbImpute


###Tphpmf
Vardat <- cbind(rnorm(100, 0, 1), condition) 
colnames(Vardat) <- c("rand_cov", "condition") 
Vardat <- cbind(1:100, Vardat) 
colnames(Vardat)[1] <- "Sample.ID"  
OTUdat2 <- 10^(bhmpf_pre_data) 
taxanames <- unlist(lapply(1:dim(sim_tab)[2], FUN = function(x) {
paste("taxa", x, sep = "")
}))  
print(taxanames)
colnames(OTUdat2) <- taxanames 
OTUdat2 <- cbind(1:100, OTUdat2)  
colnames(OTUdat2)[1] <- "Sample.ID"  
Vardat <- as.data.frame(Vardat)
OTUdat2 <- as.data.frame(OTUdat2)
otu.tab <- OTUdat2[, -1]
otu.tab<-t(otu.tab)
meta <- Vardat
colnames(meta)[2] <- "rand_cov"  
colnames(meta)[3] <- "condition"  
linda.obj2 <- linda(
otu.tab = otu.tab,        
meta = meta,              
formula = '~condition+rand_cov',  
alpha = 0.05,             
prev.cut = 0.001,           
lib.cut = 0.01,          
winsor.quan = NULL       
)
str(linda.obj2)
significant_otus2 <- rownames(linda.obj2$output[[1]])[which(linda.obj2$output[[1]]$reject)]
linda_detected2<-significant_otus2
linda_detected2 <- as.numeric(gsub("taxa", "", linda_detected2))
linda_precision2 <- sum((linda_detected2 %in% truth_DA)) / length(linda_detected2)
linda_recall2 <- sum((linda_detected2 %in% truth_DA)) / length(truth_DA)
linda_F1_score2 <- 2 * (linda_precision2 * linda_recall2) / (linda_precision2 + linda_recall2)
linda_metrics_TphPMF <- c(linda_precision2, linda_recall2, linda_F1_score2)
linda_metrics_TphPMF



###mbDenoise
Vardat <- cbind(rnorm(100, 0, 1), condition)  
colnames(Vardat) <- c("rand_cov", "condition")  
Vardat <- cbind(1:100, Vardat) 
colnames(Vardat)[1] <- "Sample.ID"  

OTUdat3 <- 10^(mbDenoise_imputed) 
taxanames <- unlist(lapply(1:dim(sim_tab)[2], FUN = function(x) {
  paste("taxa", x, sep = "")
})) 
print(taxanames)
colnames(OTUdat3) <- taxanames  
OTUdat3 <- cbind(1:100, OTUdat3)  
colnames(OTUdat3)[1] <- "Sample.ID"  
Vardat <- as.data.frame(Vardat)
OTUdat3 <- as.data.frame(OTUdat3)
otu.tab <- OTUdat3[, -1]
otu.tab<-t(otu.tab)
meta <- Vardat
colnames(meta)[2] <- "rand_cov"  
colnames(meta)[3] <- "condition"  
linda.obj3 <- linda(
  otu.tab = otu.tab,        
  meta = meta,              
  formula = '~condition+rand_cov', 
  alpha = 0.05,            
  prev.cut = 0.001,           
  lib.cut = 0.01,          
  winsor.quan = NULL        
)
str(linda.obj3)
significant_otus3 <- rownames(linda.obj3$output[[1]])[which(linda.obj3$output[[1]]$reject)]
print(significant_otus3)
linda_detected3<-significant_otus3
linda_detected3 <- as.numeric(gsub("taxa", "", linda_detected3))
linda_precision3 <- sum((linda_detected3 %in% truth_DA)) / length(linda_detected3)
linda_recall3 <- sum((linda_detected3 %in% truth_DA)) / length(truth_DA)
linda_F1_score3 <- 2 * (linda_precision3 * linda_recall3) / (linda_precision3 + linda_recall3)
linda_metrics_mbDenoise <- c(linda_precision3, linda_recall3, linda_F1_score3)
linda_metrics_mbDenoise


################metagenomeSeq#########################
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
significant_indices <- which(adjusted_pvalues < 0.05)
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
significant_indices2 <- which(adjusted_pvalues2 < 0.05)
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
significant_indices3 <- which(adjusted_pvalues3 < 0.05)
significant_indices_unnamed3 <- unname(significant_indices3)
metagenomeSeq_BHPMF_precision <- (sum(significant_indices_unnamed3 %in% truth_DA) ) / (length(significant_indices_unnamed3) )
metagenomeSeq_BHPMF_recall <- (sum(significant_indices_unnamed3 %in% truth_DA) ) / (length(truth_DA))
metagenomeSeq_BHPMF_F1_score <- 2 * (metagenomeSeq_BHPMF_precision * metagenomeSeq_BHPMF_recall) / (metagenomeSeq_BHPMF_precision + metagenomeSeq_BHPMF_recall)

metagenomeSeq_BHPMF_metrics <- c(metagenomeSeq_BHPMF_precision, metagenomeSeq_BHPMF_recall, metagenomeSeq_BHPMF_F1_score)
metagenomeSeq_BHPMF_metrics


###mbDenoise
sim_tab4 <- as.matrix(apply(t(mbDenoise_imputed), 2, as.numeric))
sim_tab4 <- as.data.frame(sim_tab4)

colnames(sim_tab4) <- 1:100
mr4 <- newMRexperiment(counts = sim_tab4, phenoData = AnnotatedDataFrame(Vardat))
mr_norm4 <- cumNorm(mr4)
conditions <- factor(Vardat$condition)
fit4 <- fitZig(obj = mr_norm4, mod = model.matrix(~ conditions),useCSSoffset = FALSE,
               zeroInflation = FALSE,
               distribution = "gaussian")
summary(fit4)
library(limma)
fit.eb4 <- eBayes(fit4@fit)
pvalues4 <- fit.eb4$p.value
p.values4 <- pvalues4[, "conditions1"]
print(pvalues4)
adjusted_pvalues4 <- p.adjust(p.values4, method = "BH")
significant_indices4 <- which(p.values4 < 0.05)
significant_indices_unnamed4 <- unname(significant_indices4)
metagenomeSeq_precision4 <- (sum(significant_indices_unnamed4 %in% truth_DA) ) / (length(significant_indices_unnamed4) )
metagenomeSeq_recall4 <- (sum(significant_indices_unnamed4 %in% truth_DA) ) / (length(truth_DA))
metagenomeSeq_F1_score4 <- 2 * (metagenomeSeq_precision4 * metagenomeSeq_recall4) / (metagenomeSeq_precision4 + metagenomeSeq_recall4)
metagenomeSeq_mbDenoise_metrics <- c(metagenomeSeq_precision4, metagenomeSeq_recall4, metagenomeSeq_F1_score4)
metagenomeSeq_mbDenoise_metrics



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
OTUdat <- floor(10^sim_tabbb +1.01)
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
OTUdat <- floor(10^sim_tabbb1 +1.01)

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

###mbDenoise

sim_tabbb3 <- t(mbDenoise_imputed)
sim_tabbb3 <- floor(10^sim_tabbb3 +1.01)
physeq3 <- phyloseq(otu_table(sim_tabbb3, taxa_are_rows = TRUE),
sample_data(Vardat))
sample_data(physeq3)$condition <- factor(sample_data(physeq3)$condition)
Deseq2_obj <- phyloseq_to_deseq2(physeq3, ~ condition)
results <- DESeq(Deseq2_obj, test="Wald", fitType="parametric")
dds_results <- results(results)
yyresults <- dds_results$pvalue
significant_indices_DESeq2 <- which(yyresults < 0.05)
significant_indices_DESeq2_unnamed3 <- unname(significant_indices_DESeq2 )
DESeq2_precision3 <- (sum(significant_indices_DESeq2_unnamed3 %in% truth_DA) ) / (length(significant_indices_DESeq2_unnamed3) )
DESeq2_recall3 <- (sum(significant_indices_DESeq2_unnamed3 %in% truth_DA) ) / (length(truth_DA))
DESeq2_F1_score3 <- 2 * (DESeq2_precision3 * DESeq2_recall3) / (DESeq2_precision3 + DESeq2_recall3)
DESeq2_mbDenoise_metrics <- c(DESeq2_precision3, DESeq2_recall3, DESeq2_F1_score3)
DESeq2_mbDenoise_metrics



library(ggplot2)
data <- data.frame(
  Method = rep(c('Wilcoxon', 'ANCOM-BC2', 'metagenomeSeq', 'DESeq2-phyloseq', 'LOCOM', 'LinDA'), each = 12),
  Metric = factor(rep(c('Precision', 'Recall', 'F1_Score'), times = 24),
                  levels = c('Precision', 'Recall', 'F1_Score')),
  Value = c(
    Wilcox_precision, Wilcox_recall, Wilcox_F1_score,
    Wilcox_mbImpute_precision, Wilcox_mbImpute_recall, Wilcox_mbImpute_F1_score,
    Wilcox_BHPMF_precision, Wilcox_BHPMF_recall, Wilcox_BHPMF_F1_score,
    Wilcox_mbDenoise_precision, Wilcox_mbDenoise_recall, Wilcox_mbDenoise_F1_score,
    ancom_precision, ancom_recall, ancom_F1_score,
    ancom_precision_mbImpute, ancom_recall_mbImpute, ancom_F1_score_mbImpute,
    ancom_precision_TphPMF, ancom_recall_TphPMF, ancom_F1_score_TphPMF,
    ancom_precision3, ancom_recall3, ancom_F1_score3,
    metagenomeSeq_precision, metagenomeSeq_recall, metagenomeSeq_F1_score,
    metagenomeSeq_mbImpute_precision, metagenomeSeq_mbImpute_recall, metagenomeSeq_mbImpute_F1_score,
    metagenomeSeq_BHPMF_precision, metagenomeSeq_BHPMF_recall, metagenomeSeq_BHPMF_F1_score,
    metagenomeSeq_precision4, metagenomeSeq_recall4, metagenomeSeq_F1_score4,
    DESeq2_precision, DESeq2_recall, DESeq2_F1_score,
    DESeq2_mbimpute_precision, DESeq2_mbimpute_recall, DESeq2_mbimpute_F1_score,
    DESeq2_BHPMF_precision, DESeq2_BHPMF_recall, DESeq2_BHPMF_F1_score,
    DESeq2_precision3, DESeq2_recall3, DESeq2_F1_score3,
    locom_precision, locom_recall, locom_F1_score,
    locom_precision_mbImpute, locom_recall_mbImpute, locom_F1_score_mbImpute,
    locom_precision_TphPMF, locom_recall_TphPMF, locom_F1_score_TphPMF,
    locom_precision_mbDenoise, locom_recall_mbDenoise, locom_F1_score_mbDenoise,
    linda_precision, linda_recall, linda_F1_score,
    linda_precision1, linda_recall1, linda_F1_score1,
    linda_precision2, linda_recall2, linda_F1_score2,
    linda_precision3, linda_recall3, linda_F1_score3
  ),
  Imputation = factor(rep(c('DA method', 'mbImpute + DA method', 'TphPMF + DA method', 'mbDenoise + DA method'), each = 3, times = 6),
                      levels = c('DA method', 'mbImpute + DA method', 'TphPMF + DA method', 'mbDenoise + DA method'))
)
pdf("/Users/hanxinyu/Desktop/mbimpute1/DA_plots.pdf", width = 8, height = 6)
for (method in unique(data$Method)) {
  method_data <- subset(data, Method == method)
  p <- ggplot(method_data, aes(x = Metric, y = Value, fill = Imputation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.88)) +
    scale_y_continuous(limits = c(0, 1.00), breaks = seq(0, 1, 0.25)) +
    scale_fill_manual(values = c("DA method" = "lightblue",
                                 "mbImpute + DA method" = "dodgerblue",
                                 "TphPMF + DA method" = "steelblue",
                                 "mbDenoise + DA method" = "darkblue")) +
    ggtitle(method) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.text = element_text(size = 16),
      axis.text.y = element_text(size = 18), 
      axis.title = element_text(size = 18),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10)
    ) +
    labs(x = "", y = "", fill = "Imputation Method")
  print(p)
}

dev.off()


