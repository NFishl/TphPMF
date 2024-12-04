T2D_real <- read.csv("/Users/hanxinyu/Desktop/otu_real_data_T2D-3.csv", row.names = "X")
control_real <- read.csv("/Users/hanxinyu/Desktop/otu_real_data_control-2.csv", row.names = "X")
D <- read.csv("/Users/hanxinyu/Desktop/D.csv", row.names = "X")
D [] <- lapply(D , as.numeric)
meta_data_T2D <- read.csv("/Users/hanxinyu/Desktop/meta_data_T2D-2.csv", row.names = "X")
meta_data_control <- read.csv("/Users/hanxinyu/Desktop/meta_data_control-2.csv", row.names = "X")


real_data <- rbind(T2D_real, control_real)
real_data[] <- lapply(real_data, as.numeric)
real_data_zi_rate <- apply(real_data, 2, FUN = function(x){
  sum(x <(log10(1.01)+0.01))/length(x)
})
chosen_taxa <- which(real_data_zi_rate <0.85)
print(length(chosen_taxa))

meta_real_data <- rbind(meta_data_T2D, meta_data_control)
meta_real_data[] <- lapply(meta_real_data, as.numeric)


filtered_data <- real_data[, chosen_taxa]
non_zero_counts <- apply(filtered_data, 1, function(row) sum(row > (log10(1.01) + 0.01)))
non_zero_proportions <- non_zero_counts / length(chosen_taxa)

selected_samples <- which(non_zero_proportions >= 0.5)
print(length(selected_samples))

selected_data <- filtered_data[selected_samples, ]
non_zero_values <- selected_data[selected_data > (log10(1.01)+ 0.01)]
mean_value <- mean(non_zero_values)
sd_value <- sd(non_zero_values)
summary(non_zero_values)

replace_zeros <- function(data, mean_val, sd_val) {
  apply(data, 1, function(row) {
    row[row <= (log10(1.01) + 0.01)] <- rnorm(sum(row <= (log10(1.01) + 0.01)), 
                                              mean = mean_val, 
                                              sd = sd_val)
    return(row)
  })
}

filled_data <- replace_zeros(selected_data, mean_value, sd_value)
simulated2<-t(filled_data)
meta_simulated2 <-meta_real_data[selected_samples, ]
D_sim <- D[chosen_taxa, chosen_taxa]

dist_obj <- as.dist(D_sim)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 2)
genus_clusters <- cutree(hc, h = 8)
family_clusters <- cutree(hc, h = 15)

microbe_data <- data.frame(
  plant_id = 1:146,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)

microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
hierarchy.info<-microbe_data


# introduce zeros based on real data.
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
      return(list("d" = y < log10(1.01) + 10^(-3) ))
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

otable <- real_data
meta_tab <- meta_real_data
# read in the gamma_norm_mix function
mean_record <- c()
percentage_record <- c()
D_vale_record <- c()
beta_record <- list()
for(j in 1:dim(otable)[2]){
  result <- gamma_norm_mix(otable[,j], data.matrix(meta_tab))
  beta_record[[j]] = result$cov_par
  mean_record <- c(mean_record, mean(otable[which(result$d < 0.5),j]))
  percentage_record <- c(percentage_record, sum(result$d > 0.5)/dim(otable)[1])
  D_vale_record <- c(D_vale_record, result$d)
}

missing_rate <- function(mean, emp_mean, emp_miss){
  win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/5
  mean_up <- mean + win_len 
  mean_lo <- mean - win_len
  sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
}
# missing_rate(1, mean_record, percentage_record)

y_sim <- simulated2
y_preserve <- y_sim
col_mean <- colMeans(y_sim)
zero_rate <- unlist(lapply(col_mean, FUN = function(x){
  print(x)
  return(missing_rate(x, mean_record, percentage_record))
}))
zero_rate[zero_rate > 0.9] = 0.9
#zero_rate[zero_rate > 0.8] = 0.8
n = dim(y_sim)[1]
m = dim(y_sim)[2]
zero_mat <- matrix(NA, nrow = n, ncol = m)
for(i in 1:m){
  zero_mat[,i] = rbinom(n,1, 1-zero_rate[i])
}
sim_tab_zi = y_sim * zero_mat


###mbImpute
library(mbImpute)
library(glmnet)
library(Matrix)
sim_tab_zi_mb <- sim_tab_zi
otu_tab <- sim_tab_zi_mb
otu_tab[1:6, 1:6]
D <- D_sim
D[1:6, 1:6]
meta_tab <- meta_simulated2
meta_tab[, 1] <- 0  
matching_rows <- which(rownames(meta_tab) %in% rownames(meta_data_T2D))
meta_tab[matching_rows, 1] <- 1

meta_data <- meta_tab

study_condition = meta_data[,1]
meta_data <- as.data.frame(unclass(meta_data))
meta_data <- meta_data[,-1]
imputed_count_mat_list <- mbImpute(condition = study_condition, otu_tab = otu_tab, metadata = meta_data, D = D)
mbImpute_mat <- imputed_count_mat_list$imp_count_mat_lognorm

###mbDenoise
library(mbDenoise)
sim_tab_zi_matrix <- as.matrix(sim_tab_zi)
result <- ZIPPCApn(sim_tab_zi_matrix , family = "negative.binomial", n.factors = 2, rank = TRUE)
mbDenoise_imputed <- result$muz

##TphPMF
sim_tab_zi[sim_tab_zi == 0] <- NA
sim_tab_zi_matrix <- as.matrix(sim_tab_zi)
sim_tab_zi_transposed <- t(sim_tab_zi_matrix)
rownames(sim_tab_zi_transposed) <- as.character(1:146)

trait.info<- sim_tab_zi_transposed
head(trait.info)
library(BHPMF)
setwd("/Users/hanxinyu/Desktop/mbimpute1/Simulation3_K_1")
tmp.dir<-dirname("/Users/hanxinyu/Desktop/mbimpute1/Simulation3_K_1/tmp/")
nrow(hierarchy.info)==nrow(trait.info)

GapFilling(trait.info,hierarchy.info,num.samples=10, burn=3, gaps=2,
           num.latent.feats=1,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir)
data <- read.table("/Users/hanxinyu/Desktop/mbimpute1/Simulation3_K_1/mean_gap_filled.txt", header = TRUE, sep = "\t")
bhmpf_pre <- as.matrix(data)
bhmpf_pre_data<-t(bhmpf_pre)

sqrt(sum((sim_tab_zi - y_sim)^2))
sim_tab_zi_mbImpute <- sim_tab_zi

write.csv(sim_tab_zi, "simulated_zi_matrix_mbImpute.csv")
write.csv(meta_simulated2, "simulated_meta_data_mbImpute.csv")
write.csv(D_sim, "D_sim.csv")

dim(read.csv("simulated_zi_matrix_mbImpute.csv", row.names = "X"))

sim_tab_zi <- 10^(sim_tab_zi) - 1.01
sim_tab_zi_trans <- t(sim_tab_zi)
sim_tab_zi_trans[sim_tab_zi_trans<0]<-0
write.csv(sim_tab_zi_trans, "simulated_zi_matrix.csv")

write.csv(t(sim_tab_zi_trans), "simulated_zi_matrix_Magic.csv")

y_sim_rec <- y_sim
sim_tab_zi_rec <- sim_tab_zi

write.csv(y_sim_rec, "truth.csv")
# alra
randomized.svd <- function(A,K, q, method = 'rsvd', mkl.seed = -1) {
  out <- setNames(vector("list", 3), c("u", "d", "v"))
  if (method == 'rsvd') {
    library(rsvd)
    out <- rsvd(A,K,q=q)
  }else if (method == 'rsvd-mkl') {
    library(fastRPCA)
    fastPCAOut <- fastPCA(A, k=K, its=q, l=(K+10), seed=mkl.seed)
    out$u <- fastPCAOut$U
    out$v <- fastPCAOut$V
    out$d <- diag(fastPCAOut$S)
  }else{
    stop('Method not recognized')
  }
  return(out)
}

normalize_data <- function (A) {
  #  Simple convenience function to library and log normalize a matrix
  
  totalUMIPerCell <- rowSums(A);
  if (any(totalUMIPerCell == 0)) {
    toRemove <- which(totalUMIPerCell == 0)
    A <- A[-toRemove,]
    totalUMIPerCell <- totalUMIPerCell[-toRemove]
    cat(sprintf("Removed %d cells which did not express any genes\n", length(toRemove)))
  }
  
  A_norm <- sweep(A, 1, totalUMIPerCell, '/');
  A_norm <- A_norm * 10E3
  A_norm <- log(A_norm +1);
}

choose_k <- function (A_norm,K=100, thresh=6, noise_start=80,q=2,use.mkl=F, mkl.seed =-1) {
  #  Heuristic for choosing rank k for the low rank approximation based on
  #  statistics of the spacings between consecutive singular values. Finds
  #  the smallest singular value \sigma_i such that $\sigma_i - \sigma_{i-1}
  #  is significantly different than spacings in the tail of the singular values.
  # 
  #
  # Args:
  #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
  #   K: Number of singular values to compute. Must be less than the
  #   smallest dimension of the matrix.
  #   thresh: Number of standard deviations away from the ``noise'' singular
  #   values which you consider to be signal
  #   noise_start : Index for which all smaller singular values are
  #   considered noise
  #   q : Number of additional power iterations
  #   use.mkl : Use the Intel MKL based implementation of SVD. Needs to be
  #             installed from https://github.com/KlugerLab/rpca-mkl.
  #   mkl.seed : Only relevant if use.mkl=T. Set the seed for the random
  #   generator for the Intel MKL implementation of SVD. Any number <0 will
  #   use the current timestamp. If use.mkl=F, set the seed using
  #   set.seed() function as usual.
  # Returns:
  #   A list with three items
  #       1) Chosen k
  #       2) P values of each possible k 
  #       3) Singular values of the matrix A_norm
  
  if (K > min(dim(A_norm))) {
    stop("For an m by n matrix, K must be smaller than the min(m,n).\n")
  }
  if (noise_start >K-5) {
    stop("There need to be at least 5 singular values considered noise.\n")
  }
  noise_svals <- noise_start:K
  if (!use.mkl) {
    rsvd_out <- randomized.svd(A_norm,K,q=q)
  }else {    
    rsvd_out <- randomized.svd(A_norm,K,q=q, method='rsvd-mkl', mkl.seed=mkl.seed)
  }
  diffs <- rsvd_out$d[1:(length(rsvd_out$d)-1)] - rsvd_out$d[2:length(rsvd_out$d)]
  mu <- mean(diffs[noise_svals-1])
  sigma <- sd(diffs[noise_svals-1])
  num_of_sds <- (diffs-mu)/sigma
  k <- max (which(num_of_sds > thresh))
  return (list( k=k,num_of_sds = num_of_sds,d=rsvd_out$d))
}

alra <- function( A_norm, k=0,q=10, quantile.prob = 0.001, use.mkl = F, mkl.seed=-1) {
  cat(sprintf("Read matrix with %d cells and %d genes\n", nrow(A_norm), ncol(A_norm)))
  #if (class(A_norm) != 'matrix') {
  #stop(sprintf("A_norm is of class %s, but it should be of class matrix. Did you forget to run as.matrix()?",class(A_norm)))
  #}
  
  if (k ==0 ) {
    k_choice <- choose_k(A_norm)
    k <-  k_choice$k
    cat(sprintf("Chose k=%d\n",k))
  }
  
  cat("Getting nonzeros\n")
  originally_nonzero <- A_norm >0 
  
  cat("Randomized SVD\n")
  if (!use.mkl) {
    fastDecomp_noc <- randomized.svd(A_norm,k,q=q)
  }else {
    fastDecomp_noc <- randomized.svd(A_norm,k,q=q, method='rsvd-mkl', mkl.seed=mkl.seed)
  }
  A_norm_rank_k <- fastDecomp_noc$u[,1:k]%*%diag(fastDecomp_noc$d[1:k])%*% t(fastDecomp_noc$v[,1:k])
  
  
  cat(sprintf("Find the %f quantile of each gene\n", quantile.prob))
  #A_norm_rank_k_mins <- abs(apply(A_norm_rank_k,2,min))
  A_norm_rank_k_mins <- abs(apply(A_norm_rank_k,2,FUN=function(x) quantile(x,quantile.prob)))
  cat("Sweep\n")
  A_norm_rank_k_cor <- replace(A_norm_rank_k, A_norm_rank_k <= A_norm_rank_k_mins[col(A_norm_rank_k)], 0)
  
  
  sd_nonzero <- function(x) sd(x[!x == 0])
  sigma_1 <- apply(A_norm_rank_k_cor, 2, sd_nonzero)
  sigma_2 <- apply(A_norm, 2, sd_nonzero)
  mu_1 <- colSums(A_norm_rank_k_cor)/colSums(!!A_norm_rank_k_cor)
  mu_2 <- colSums(A_norm)/colSums(!!A_norm)
  
  toscale <- !is.na(sigma_1) & !is.na(sigma_2) & !(sigma_1 == 0 & sigma_2 == 0) & !(sigma_1 == 0)
  cat(sprintf("Scaling all except for %d columns\n", sum(!toscale)))
  
  sigma_1_2 <- sigma_2/sigma_1
  toadd  <- -1*mu_1*sigma_2/sigma_1 + mu_2
  
  A_norm_rank_k_temp <- A_norm_rank_k_cor[,toscale]
  A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp,2, sigma_1_2[toscale],FUN = "*")
  A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp,2, toadd[toscale],FUN = "+")
  
  A_norm_rank_k_cor_sc <- A_norm_rank_k_cor
  A_norm_rank_k_cor_sc[,toscale] <- A_norm_rank_k_temp
  A_norm_rank_k_cor_sc[A_norm_rank_k_cor==0] = 0
  
  lt0 <- A_norm_rank_k_cor_sc  <0
  A_norm_rank_k_cor_sc[lt0] <- 0 
  cat(sprintf("%.2f%% of the values became negative in the scaling process and were set to zero\n", 100*sum(lt0)/(nrow(A_norm)*ncol(A_norm))))
  
  A_norm_rank_k_cor_sc[originally_nonzero & A_norm_rank_k_cor_sc ==0] <- A_norm[originally_nonzero & A_norm_rank_k_cor_sc ==0]
  
  colnames(A_norm_rank_k_cor) <- colnames(A_norm)
  colnames(A_norm_rank_k_cor_sc) <- colnames(A_norm)
  colnames(A_norm_rank_k) <- colnames(A_norm)
  
  original_nz <- sum(A_norm>0)/(nrow(A_norm)*ncol(A_norm))
  completed_nz <- sum(A_norm_rank_k_cor_sc>0)/(nrow(A_norm)*ncol(A_norm))
  cat(sprintf("The matrix went from %.2f%% nonzero to %.2f%% nonzero\n", 100*original_nz, 100*completed_nz))
  
  list(A_norm_rank_k=A_norm_rank_k,A_norm_rank_k_cor =A_norm_rank_k_cor, A_norm_rank_k_cor_sc=A_norm_rank_k_cor_sc)
}
si_impute <- function(sim_tab){
  set.seed(12345)
  y_sim = sim_tab
  # add filter
  filter <- lapply(X = 1:ncol(y_sim), FUN = function(col_i){
    y = y_sim[,col_i]
    n = length(y)
    nz <- sum(y <= (log10(1.01) + 1e-6))
    pz = 1 - nz/n
    test = pz - 1.96 * sqrt(pz * (1-pz)/n)
    if(nz == n || test <= 0){
      return(0)
    }else{
      return(1)
    }
  })
  y_imp <- y_sim
  #perform imputation on the rest
  filter_vec <- which(unlist(filter) == 1)
  y_sim = y_sim[, filter_vec]
  
  na_idx <- y_sim < log10(1.01) + 0.001
  y_sim[y_sim < log10(1.01) + 0.001] = NA
  
  # y_sim_cv <- unlist(y_sim)
  # y_sim_validate <- matrix(y_sim_cv, nrow(y_sim), ncol = ncol(y_sim))
  # identical(y_sim_validate, y_sim)
  y_sim_cv <- unlist(y_sim)
  na_intro <- sample(which(!is.na(y_sim_cv)), floor(sum(!is.na(y_sim_cv))/10))
  y_sim_cv_intro <- y_sim_cv
  y_sim_cv_intro[na_intro] = NA
  y_sim_cv_intro <- matrix(y_sim_cv_intro, nrow = nrow(y_sim), ncol = ncol(y_sim))
  
  j = 1
  se = 1e10
  for(i in 1:5){
    si_cv_1 <- softImpute(y_sim_cv_intro, rank.max = i, lambda = 0)
    y_imp_cv <- complete(y_sim_cv_intro, si_cv_1)
    y_sim_vali <- as.vector(y_imp_cv)
    se2 <- sum((y_sim_cv[na_intro] - y_sim_vali[na_intro])^2)
    print(se2)
    if(se2 < se){
      se = se2
      j = i
    }
  }
  print(j)
  si1 <- softImpute(y_sim_cv_intro, rank.max = 10, lambda = 0, trace.it = TRUE)
  impute_mat <- complete(y_sim_cv_intro, si1)
  
  y_imp[, filter_vec] = impute_mat
  
  return(y_imp)
}

saver_mat <- saver(real_data, ncores = 1, estimates.only = TRUE)
saver_mat_eval <- log10(saver_mat + 1.01)
saver_mat_eval <- t(saver_mat_eval)
for(j in 1:dim(saver_mat_eval)[2]){
  saver_mat_eval[,j] <- saver_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(saver_mat_eval[,j])
}
sqrt(sum((y_sim_rec- saver_mat_eval)^2))

scimpute(count_path = "simulated_zi_matrix.csv", Kcluster = 1, out_dir = "sim_imp")
scImpute_mat <- read.csv("sim_impscimpute_count.csv", row.names = "X")
dim(scImpute_mat)
dim(sim_tab_zi_mbImpute)
scImpute_mat_eval <- log10(scImpute_mat + 1.01)
scImpute_mat_eval <- t(scImpute_mat_eval)
for(j in 1:dim(scImpute_mat_eval)[2]){
  scImpute_mat_eval[,j] <- scImpute_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(scImpute_mat_eval[,j])
}
sqrt(sum((y_sim_rec- scImpute_mat_eval)^2))

T2D_mat <- normalize_data(sim_tab_zi)
k_chosen <- choose_k(T2D_mat, K = 49, noise_start = 44)$k
T2D_mat <- as.matrix(T2D_mat)
alra_mat_eval <- alra(T2D_mat, k = k_chosen)$A_norm_rank_k_cor_sc
for(j in 1:dim(alra_mat_eval)[2]){
  alra_mat_eval[,j] <- alra_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(alra_mat_eval[,j])
}
sqrt(sum((y_sim_rec - alra_mat_eval)^2)) 

# softImpute
softImpute_mat = si_impute(sim_tab_zi)
softImpute_mat[softImpute_mat < 0] = 0
softImpute_mat_eval <- log10(softImpute_mat+1.01)
for(j in 1:dim(softImpute_mat_eval)[2]){
  softImpute_mat_eval[,j] <- softImpute_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(softImpute_mat_eval[,j])
}
sqrt(sum((y_sim_rec - softImpute_mat_eval)^2))

bt_mse <- matrix(ncol = 800, nrow = 8)
set.seed(1234)
for(i in 1:800){
  bt_idx <-  sample(1:146, 146, replace = TRUE)
  bt_mse[,i] <- c(sum((y_sim_rec[,bt_idx] - bhmpf_pre_data[,bt_idx])^2),
                  sum((y_sim_rec[,bt_idx] - mbImpute_mat[,bt_idx])^2),
                  sum((y_sim_rec[,bt_idx] - mbDenoise_imputed[,bt_idx])^2),
                  sum((y_sim_rec[,bt_idx] - scImpute_mat_eval[,bt_idx])^2),
                  sum((y_sim_rec[,bt_idx] - alra_mat_eval[,bt_idx])^2),
                  sum((y_sim_rec[,bt_idx] - saver_mat_eval[,bt_idx])^2),
                  sum((y_sim_rec[,bt_idx] - softImpute_mat_eval[,bt_idx])^2),
                  sum((y_sim_rec[,bt_idx] - sim_tab_zi_mbImpute[,bt_idx])^2)) / (146 * 60)
}
apply(bt_mse, 1, mean)
mse11_vector<-apply(bt_mse, 1, mean)
mse_sd <- apply(bt_mse, 1, sd)
mse_se <- mse_sd / sqrt(800)
df <- data.frame(
  MSE = mse11_vector,
  method = c("TphPMF", "mbImpute", "mbDenoise", "scImpute", "ALRA", "SAVER", "softImpute", "No imputation"),
  SE = mse_se
)

df$method <- factor(df$method, levels = c("TphPMF", "mbImpute", "mbDenoise", "scImpute", "ALRA", "SAVER", "softImpute", "No imputation"))
pdf("/Users/hanxinyu/Desktop/mbimpute1/mse_comparision3.pdf")

ggplot(df, aes(x=method, y=MSE, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.9) +  
  geom_errorbar(aes(ymin=MSE - SE, ymax=MSE + SE), width=0.2, position=position_dodge(0.7)) +  
  ylab("MSE") +
  xlab("Imputation method") +
  ggtitle("Simulation 3") +  
  scale_y_continuous(limits = c(0, 7.5)) +  
  scale_fill_manual("method", values = c(
    "No imputation" = "#fb9a99", 
    "ALRA" = "#ece2f0", 
    "TphPMF" = "#a6cee3", 
    "scImpute" = "#1f78b4", 
    "SAVER" = "#b2df8a", 
    "softImpute" = "#6A95B6", 
    "mbImpute" = "#1D9E78", 
    "mbDenoise" = "#FF6347"  
  )) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20), 
    axis.text.x = element_text(size = 7),    
    axis.text.y = element_text(size = 30)    
  )

dev.off()


cor <- matrix(NA, nrow = 8, ncol = dim(mbImpute_mat)[2])
rownames(cor) <- c("TphPMF", "mbImpute", "mbDenoise","scImpute", "ALRA", "SAVER", "softImpute", "No imputation")
for(i in 1:dim(mbImpute_mat)[2]){
  cor[1,i] <- cor(y_sim_rec[,i], bhmpf_pre_data[,i])
  cor[2,i] <- cor(y_sim_rec[,i], mbImpute_mat[,i])
  cor[3,i] <- cor(y_sim_rec[,i], mbDenoise_imputed[,i])
  cor[4,i] <- cor(y_sim_rec[,i], scImpute_mat_eval[,i])
  cor[5,i] <- cor(y_sim_rec[,i], alra_mat_eval[,i])
  cor[6,i] <- cor(y_sim_rec[,i], saver_mat_eval[,i])
  cor[7,i] <- cor(y_sim_rec[,i], softImpute_mat_eval[,i])
  cor[8,i] <- cor(y_sim_rec[,i], sim_tab_zi_mbImpute[,i])
  
}
order <- sort(cor[1,], decreasing = TRUE, index.return = TRUE)$ix

bt_cor <- matrix(ncol = 800, nrow = 8)
set.seed(1234)
for(i in 1:800){
  bt_idx <-  sample(1:146, 146, replace = TRUE)
  bt_cor[,i] <- rowMeans(cor[,bt_idx])
}
pearson_cor_vector<-apply(bt_cor, 1, mean)
print(pearson_cor_vector)
apply(bt_cor, 1, sd)
pearson_cor_mean <- apply(bt_cor, 1, mean)
pearson_cor_sd <- apply(bt_cor, 1, sd)
df <- data.frame(
  method = c("TphPMF", "mbImpute", "mbDenoise", "scImpute", "ALRA", "SAVER", "softImpute", "No imputation"),
  correlation = pearson_cor_mean,
  SE = pearson_cor_se
)

df$method <- factor(df$method, levels = c("TphPMF", "mbImpute", "mbDenoise", "scImpute", "ALRA", "SAVER", "softImpute", "No imputation"))

library(ggplot2)

pdf("/Users/hanxinyu/Desktop/mbimpute1/correlation_analysis3.pdf")

ggplot(df, aes(x=method, y=correlation, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.9) +  
  geom_errorbar(aes(ymin=correlation - SE, ymax=correlation + SE), width=0.2, position=position_dodge(0.7)) +  
  ylab("Correlation") +
  xlab("Imputation method") +
  ggtitle("Simulation 3") +  
  scale_y_continuous(limits = c(0, 0.5)) +  
  scale_fill_manual("method", values = c(
    "No imputation" = "#fb9a99", 
    "ALRA" = "#ece2f0", 
    "TphPMF" = "#a6cee3", 
    "scImpute" = "#1f78b4", 
    "SAVER" = "#b2df8a", 
    "softImpute" = "#6A95B6", 
    "mbImpute" = "#1D9E78", 
    "mbDenoise" = "#FF6347"  
  )) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    axis.text.x = element_text(size = 6),    
    axis.text.y = element_text(size = 30)    
  )

dev.off()
