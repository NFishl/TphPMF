control_coeff <- readRDS("/Users/hanxinyu/Desktop/control_dat1_sim_add_filter_coef.rds")
T2D_coeff <- readRDS("/Users/hanxinyu/Desktop/T2D_dat2_sim_add_filter_coef.rds")
otu_real_data_T2D <- read.csv("/Users/hanxinyu/Desktop/otu_real_data_T2D-3.csv", row.names = "X")
otu_real_data_control <- read.csv("/Users/hanxinyu/Desktop/otu_real_data_control-2.csv", row.names = "X")
D <- read.csv("/Users/hanxinyu/Desktop/D.csv", row.names = "X")
meta_data_T2D <- read.csv("/Users/hanxinyu/Desktop/meta_data_T2D-2.csv", row.names = "X")
meta_data_control <- read.csv("/Users/hanxinyu/Desktop/meta_data_control-2.csv", row.names = "X")
set.seed(12345)
y_sim = otu_real_data_T2D
x = meta_data_T2D
k = 5
c1 = T2D_coeff
c1 = c1*5
c1[1] = T2D_coeff[1]
c1[length(c1)-6] = 0.2
c1[length(c1)-11] = -0.3

#loading used functions
design_mat_row_gen2 <- function(count_mat, covariate_mat, row_index, col_index, close_taxa){
  n = dim(count_mat)[1]
  m = dim(count_mat)[2]
  k = length(close_taxa[[1]])
  if(is.vector(covariate_mat)){
    p = 1
    i = row_index
    j = col_index
    
    close_taxa_set <- close_taxa[[j]]
    #generate a row including response and a row of design matrix.
    row_gen <- rep(0, m*k + (n-1) * n + n*p)
    
    row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
    row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
  }else{
    p = dim(covariate_mat)[2]
    i = row_index
    j = col_index
    
    close_taxa_set <- close_taxa[[j]]
    #generate a row including response and a row of design matrix.
    row_gen <- rep(0, m*k + (n-1) * n + n*p)
    
    row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
    row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
  }
  return(row_gen)
}

#identify the low abundance taxa by binomial test
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

filter_vec <- which(unlist(filter) == 1)
y_sim = y_sim[, filter_vec]
D = D[filter_vec,filter_vec]
dist_obj <- as.dist(D)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 2)
genus_clusters <- cutree(hc, h = 10)
family_clusters <- cutree(hc, h = 15)
microbe_data <- data.frame(
  plant_id = 1:193,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)
microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
hierarchy.info<-microbe_data
y_preserve <- y_sim
m = dim(y_sim)[2]
n = dim(y_sim)[1]

X <- as.matrix(x)
#generate close set
close_taxa <- list()
for(j in 1:m){
  close_dist <- D[D[,j] %in% sort(D[,j])[2:(k+1)],j]
  close_taxa_vec = which(D[,j] %in% close_dist[close_dist != max(close_dist)])
  if(length(close_taxa_vec) < k){
    close_taxa_vec <- c(close_taxa_vec, which(D[,j] == max(close_dist))[1:(k-length(close_taxa_vec))])
  }
  close_taxa[[j]] = close_taxa_vec
}
#generate design matrix
p = dim(x)[2]
if(is.null(p)){
  p = 1
}
#generate a row including response and a row of design matrix.
row_length <- m * k + (n-1) * n + n*p
idx_set <- matrix(NA, m*n, 2)
for(i in 1:n){
  for(j in 1:m){
    idx_set[(i-1)*m+j, ] <- c(i, j)
  }
}


design_mat_gen <- matrix(0, nrow = dim(idx_set)[1], ncol = row_length)
for(i in 1:dim(idx_set)[1]){
  design_mat_gen[i,] <- design_mat_row_gen2(y_sim, x[1:n,], idx_set[i,1], idx_set[i,2], close_taxa)
}
design_mat_gen <- cbind(rep(1, dim(design_mat_gen)[1]), design_mat_gen)
imputed_value <- design_mat_gen %*% c1
impute_mat <- y_sim
for(i in 1:dim(idx_set)[1]){
  impute_mat[idx_set[i,1], idx_set[i,2]] = max(imputed_value[i], log10(1.01))
}
print(sqrt(sum((impute_mat - y_sim)^2)))
y_sim <- impute_mat

y_sim <- y_sim - 1.5
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

otable <- y_preserve
meta_tab <- meta_data_T2D
D

dim(otable)
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
#filter <- which(is.nan(mean_record) | mean_record < log10(1.01) + 0.001)

# filter out the nan values.
plot(mean_record, percentage_record)

# build a map between the mean and percentage of missing
missing_rate <- function(mean, emp_mean, emp_miss){
  win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/3
  mean_up <- mean + win_len
  mean_lo <- mean - win_len
  sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
}
# missing_rate(1, mean_record, percentage_record)

col_mean <- colMeans(y_sim)
zero_rate <- unlist(lapply(col_mean, FUN = function(x){
  print(x)
  return(missing_rate(x, mean_record, percentage_record))
}))
zero_mat <- matrix(NA, nrow = n, ncol = m)
for(i in 1:m){
  zero_mat[,i] = rbinom(n,1, 1-zero_rate[i])
}
sim_tab_zi = y_sim * zero_mat
sim_tab_zi[sim_tab_zi < log10(1.01)+1e-6] = log10(1.01)


###mbImpute
library(mbImpute)
library(glmnet)
library(Matrix)
sim_tab_zi_mb <- sim_tab_zi
sim_tab_zi_mb[sim_tab_zi_mb == log10(1.01)] = 0
otu_tab <- sim_tab_zi_mb
otu_tab[1:6, 1:6]
D[1:6, 1:6]
meta_data <- meta_tab
meta_data[1:6, 1:6]
study_condition = meta_data[,1]
meta_data <- as.data.frame(unclass(meta_data))
meta_data <- meta_data[,-1]
imputed_count_mat_list <- mbImpute(condition = study_condition, otu_tab = otu_tab, metadata = meta_data, D = D)
mbImpute_mat <- imputed_count_mat_list$imp_count_mat_lognorm


sqrt(sum((sim_tab_zi - y_preserve)^2))
sim_tab_zi_mbImpute <- sim_tab_zi

write.csv(sim_tab_zi, "simulated_zi_matrix_mbImpute.csv")

sim_tab_zi <- 10^(sim_tab_zi) - 1.01
sim_tab_zi_trans <- t(sim_tab_zi)
write.csv(sim_tab_zi_trans, "simulated_zi_matrix.csv")

y_sim_rec <- y_sim
sim_tab_zi_rec <- sim_tab_zi

write.csv(y_sim_rec, "truth.csv")

print(max(colMeans(y_sim)))
print(max(colMeans(sim_tab_zi)))
print(max(rowMeans(sim_tab_zi_trans)))

library(usethis)
library(devtools)
library(scImpute)
library(ggplot2)
library(gplots)
library(SAVER)
library(softImpute)
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

print(max(colMeans(y_sim)))
print(max(colMeans(sim_tab_zi)))
print(max(rowMeans(sim_tab_zi_trans)))

eval_omit <- c(26, 34, 42, 117, 123, 158, 191, 192)

real_data <- data.matrix(sim_tab_zi_trans)
dim(real_data)
saver_mat <- saver(real_data, ncores = 1, estimates.only = TRUE)
saver_mat_eval <- log10(saver_mat + 1.01)
saver_mat_eval <- t(saver_mat_eval)
for(j in 1:dim(saver_mat_eval)[2]){
  saver_mat_eval[,j] <- saver_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(saver_mat_eval[,j])
}
sqrt(sum((y_sim_rec[,-eval_omit] - saver_mat_eval[,-eval_omit])^2))

scimpute(count_path = "simulated_zi_matrix.csv", Kcluster = 1, out_dir = "sim_imp")
scImpute_mat <- read.csv("sim_impscimpute_count.csv", row.names = "X")
dim(scImpute_mat)
scImpute_mat_eval <- log10(scImpute_mat + 1.01)
scImpute_mat_eval <- t(scImpute_mat_eval)
for(j in 1:dim(scImpute_mat_eval)[2]){
  scImpute_mat_eval[,j] <- scImpute_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(scImpute_mat_eval[,j])
}
sqrt(sum((y_sim_rec[,-eval_omit] - scImpute_mat_eval[,-eval_omit])^2))

T2D_mat <- normalize_data(sim_tab_zi)
k_chosen <- choose_k(T2D_mat, K = 49, noise_start = 44)$k
T2D_mat <- as.matrix(T2D_mat)
alra_mat_eval <- alra(T2D_mat, k = k_chosen)$A_norm_rank_k_cor_sc
for(j in 1:dim(alra_mat_eval)[2]){
  alra_mat_eval[,j] <- alra_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(alra_mat_eval[,j])
}
sqrt(sum((y_sim_rec[,-eval_omit] - alra_mat_eval[,-eval_omit])^2))

# softImpute
softImpute_mat = si_impute(sim_tab_zi)
softImpute_mat[softImpute_mat < 0] = 0
softImpute_mat_eval <- log10(softImpute_mat+1.01)
for(j in 1:dim(softImpute_mat_eval)[2]){
  softImpute_mat_eval[,j] <- softImpute_mat_eval[,j] * max(sim_tab_zi_mbImpute[,j]) / max(softImpute_mat_eval[,j])
}
sqrt(sum((y_sim_rec[,-eval_omit] - softImpute_mat_eval[,-eval_omit])^2))

##TphPMF
sim_tab_zi[sim_tab_zi == log10(1.01)] <- NA
sim_tab_zi_matrix <- as.matrix(sim_tab_zi)
sim_tab_zi_transposed <- t(sim_tab_zi_matrix)
rownames(sim_tab_zi_transposed) <- as.character(1:193)

trait.info<- sim_tab_zi_transposed
library(usethis)
library(devtools)
library(Matrix)

library(BHPMF)
setwd("/Users/hanxinyu/Desktop/bbhpmf_sim1/")
tmp.dir<-dirname("/Users/hanxinyu/Desktop/bbhpmf_sim1/tmp/")
GapFilling(trait.info,hierarchy.info,num.samples=20, burn=5, gaps=2,
           num.latent.feats=11.5,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)

complete_data_matrix <- as.matrix(y_sim)
complete_data_transposed <- t(complete_data_matrix)

real_data<-t(complete_data_transposed)
rownames(complete_data_transposed) <- 1:193

data <- read.table("/Users/hanxinyu/Desktop/bbhpmf_sim1/mean_gap_filled.txt", header = TRUE, sep = "\t")
bhmpf_pre <- as.matrix(data)
bhmpf_pre_data<-t(bhmpf_pre)


bt_mse <- matrix(ncol = 800, nrow = 7)
set.seed(1234)
for(i in 1:800){
  bt_idx <-  sample(1:185, 185, replace = TRUE)
  bt_mse[,i] <- c(sum((y_sim_rec[,-eval_omit][,bt_idx] - bhmpf_pre_data[,-eval_omit][,bt_idx])^2),
                  sum((y_sim_rec[,-eval_omit][,bt_idx] - mbImpute_mat[,-eval_omit][,bt_idx])^2),
                  sum((y_sim_rec[,-eval_omit][,bt_idx] - scImpute_mat_eval[,-eval_omit][,bt_idx])^2),
                  sum((y_sim_rec[,-eval_omit][,bt_idx] - alra_mat_eval[,-eval_omit][,bt_idx])^2),
                  sum((y_sim_rec[,-eval_omit][,bt_idx] - saver_mat_eval[,-eval_omit][,bt_idx])^2),
                  sum((y_sim_rec[,-eval_omit][,bt_idx] - softImpute_mat_eval[,-eval_omit][,bt_idx])^2),
                  sum((y_sim_rec[,-eval_omit][,bt_idx] - sim_tab_zi_mbImpute[,-eval_omit][,bt_idx])^2)) / (185 * 53)
}
apply(bt_mse, 1, mean)
mse1_vector<-apply(bt_mse, 1, mean)
print(mse1_vector)
apply(bt_mse, 1, sd)

pdf("/Users/hanxinyu/Desktop/mbimpute/mse1_comparision1.pdf")
df <- data.frame(cbind(mse1_vector, c("TphPMF", "mbImpute", "scImpute", "ALRA", "SAVER", "softImpute", "No imputation")))
colnames(df) <- c("MSE", "method")
df$MSE <- as.numeric(as.character(df$MSE))
df$method <- factor(df$method, levels = c("TphPMF", "mbImpute", "scImpute", "ALRA", "SAVER", "softImpute", "No imputation"))
ggplot()+
  ylab("MSE") +
  xlab("Imputation method")+
  ggtitle("Simulation 2") +  
  geom_bar(data = df, aes(x=method, y=MSE, fill = method), stat="identity", position=position_dodge())+
  scale_y_continuous(limits = c(0,5)) + 
  scale_fill_manual("method", values = c("No imputation" = "#fb9a99", "ALRA" = "#ece2f0", "TphPMF" = "#a6cee3", "scImpute" = "#1f78b4", "SAVER" = "#b2df8a", "softImpute" = "#6A95B6", "mbImpute" = "#1D9E78")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none" , 
    plot.title = element_text(hjust = 0.5) )
dev.off()

cor <- matrix(NA, nrow = 7, ncol = dim(mbImpute_mat[,-eval_omit])[2])
rownames(cor) <- c("TphPMF","mbImpute", "scImpute", "ALRA", "SAVER", "softImpute", "No imputation")
idx_eval <- which(!(1:193 %in% eval_omit))
for(i in 1:dim(mbImpute_mat[,-eval_omit])[2]){
  cor[1,i] <- cor(y_sim_rec[,idx_eval[i]], bhmpf_pre_data[,idx_eval[i]])
  cor[2,i] <- cor(y_sim_rec[,idx_eval[i]], mbImpute_mat[,idx_eval[i]])
  cor[3,i] <- cor(y_sim_rec[,idx_eval[i]], scImpute_mat_eval[,idx_eval[i]])
  cor[4,i] <- cor(y_sim_rec[,idx_eval[i]], alra_mat_eval[,idx_eval[i]])
  cor[5,i] <- cor(y_sim_rec[,idx_eval[i]], saver_mat_eval[,idx_eval[i]])
  cor[6,i] <- cor(y_sim_rec[,idx_eval[i]], softImpute_mat_eval[,idx_eval[i]])
  cor[7,i] <- cor(y_sim_rec[,idx_eval[i]], sim_tab_zi_mbImpute[,idx_eval[i]])
}
which(cor[3,] > 0.95)
hist(cor[1,])
hist(cor[2,])
hist(cor[3,])
hist(cor[4,])
hist(cor[5,])
rowMeans(cor)

bt_cor <- matrix(ncol = 800, nrow = 7)
set.seed(1234)
for(i in 1:800){
  bt_idx <-  sample(1:185, 185, replace = TRUE)
  bt_cor[,i] <- rowMeans(cor[,bt_idx])
}
apply(bt_cor, 1, mean)
pearson_cor_vector<-apply(bt_cor, 1, mean)
print(pearson_cor_vector)
apply(bt_cor, 1, sd)

pdf("/Users/hanxinyu/Desktop/mbimpute/correlation1_anlaysis1.pdf")
df <- data.frame(cbind(pearson_cor_vector, c("TphPMF", "mbImpute", "scImpute", "ALRA", "SAVER", "softImpute", "No imputation")))
colnames(df) <- c("correlation", "method")
df$correlation <- as.numeric(as.character(df$correlation))
df$method <- factor(df$method, levels = c("No imputation", "SAVER", "ALRA", "scImpute", "softImpute", "mbImpute","TphPMF"))
ggplot()+
  ylab("Correlation") +
  xlab("Imputation method")+
  ggtitle("Simulation 2") +
  geom_bar(data = df, aes(x=method, y=correlation, fill = method), stat="identity", position=position_dodge())+
  scale_y_continuous(limits = c(0,0.6)) + 
  scale_fill_manual("method", values = c("No imputation" = "#fb9a99", "ALRA" = "#ece2f0",  "scImpute" = "#1f78b4", "SAVER" = "#b2df8a", "softImpute" = "#6A95B6", "mbImpute" = "#1D9E78", "TphPMF" = "#2B35C4")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none" , 
    plot.title = element_text(hjust = 0.5)  
  )
dev.off()

library(transport)
wasserstein1d_values<-c(wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(sim_tab_zi_mbImpute[,-eval_omit])/apply(sim_tab_zi_mbImpute[,-eval_omit], 2, sd)),
                        wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(alra_mat_eval[,-eval_omit])/apply(alra_mat_eval[,-eval_omit], 2, sd)),
                        wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(bhmpf_pre_data[,-eval_omit])/apply(bhmpf_pre_data[,-eval_omit], 2, sd)),
                        wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(saver_mat_eval[,-eval_omit])/apply(saver_mat_eval[,-eval_omit], 2, sd)),
                        wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(scImpute_mat_eval[,-eval_omit])/apply(scImpute_mat_eval[,-eval_omit], 2, sd)),
                        wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(softImpute_mat_eval[,-eval_omit])/apply(softImpute_mat_eval[,-eval_omit], 2, sd)),
                        wasserstein1d(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd), colMeans(mbImpute_mat[,-eval_omit])/apply(mbImpute_mat[,-eval_omit], 2, sd)))

names(wasserstein1d_values) <- c("No imputation", "ALRA", "TphPMF", "SAVER", "scImpute", "softImpute", "mbImpute")

result_df <- as.data.frame(t(wasserstein1d_values))

pdf("/Users/hanxinyu/Desktop/mbimpute/bmean_divide_sd_comparison1.pdf")

# Set up the layout for a 3x3 grid of plots
par(mfrow = c(3, 3))

# Define a function to add Wasserstein distance text to a plot
add_wasserstein_value <- function(wasserstein_value, x_pos, y_pos ,col = "black") {
  # Add the Wasserstein distance value to the plot at the specified position
  text(x_pos, y_pos, labels = paste("d = ", round(wasserstein_value, 3)), col = col,
       adj = c(1, 1), cex = 1.2, xpd = TRUE)
}

# Histogram for the complete data (no Wasserstein distance needed)
hist(colMeans(y_sim_rec[,-eval_omit])/apply(y_sim_rec[,-eval_omit], 2, sd),
     col = "#fb9a99", main = "Complete", xlab = "mean/sd")

# No imputation histogram with Wasserstein distance
hist(colMeans(sim_tab_zi_mbImpute[,-eval_omit])/apply(sim_tab_zi_mbImpute[,-eval_omit], 2, sd),
     col = "#BC6C36", main = "No imputation", xlab = "mean/sd")
add_wasserstein_value(wasserstein1d_values["No imputation"], x_pos = 14, y_pos = 130)

# Empty plot for placeholder
plot.new()

# ALRA histogram with Wasserstein distance
hist(colMeans(alra_mat_eval[,-eval_omit])/apply(alra_mat_eval[,-eval_omit], 2, sd),
     col = "#CE95CE", main = "ALRA", xlab = "mean/sd")
add_wasserstein_value(wasserstein1d_values["ALRA"], x_pos = 20, y_pos = 130)

# TphPMF histogram with Wasserstein distance
hist(colMeans(bhmpf_pre_data[,-eval_omit])/apply(bhmpf_pre_data[,-eval_omit], 2, sd),
     col = "#B1CEE3", main = "TphPMF", xlab = "mean/sd")
add_wasserstein_value(wasserstein1d_values["TphPMF"], x_pos = 15, y_pos = 35)

# SAVER histogram with Wasserstein distance
hist(colMeans(saver_mat_eval[,-eval_omit])/apply(saver_mat_eval[,-eval_omit], 2, sd),
     col = "#b2df8a", main = "SAVER", xlab = "mean/sd")
add_wasserstein_value(wasserstein1d_values["SAVER"], x_pos = 40, y_pos = 135)

# scImpute histogram with Wasserstein distance
hist(colMeans(scImpute_mat_eval[,-eval_omit])/apply(scImpute_mat_eval[,-eval_omit], 2, sd),
     col = "#1f78b4", main = "scImpute", xlab = "mean/sd")
add_wasserstein_value(wasserstein1d_values["scImpute"], x_pos = 14, y_pos = 85)

# softImpute histogram with Wasserstein distance
hist(colMeans(softImpute_mat_eval[,-eval_omit])/apply(softImpute_mat_eval[,-eval_omit], 2, sd),
     col = "#6BAED6", main = "softImpute", xlab = "mean/sd")
add_wasserstein_value(wasserstein1d_values["softImpute"], x_pos = 10, y_pos = 50)

# mbImpute histogram with Wasserstein distance
hist(colMeans(mbImpute_mat[,-eval_omit])/apply(mbImpute_mat[,-eval_omit], 2, sd),
     col = "#36A02C", main = "mbImpute", xlab = "mean/sd")
add_wasserstein_value(wasserstein1d_values["mbImpute"], x_pos = 12, y_pos = 38)

# Add an overall title
mtext("Simulation 2", outer = TRUE, cex = 1.2, line = -1.3)

# Close the PDF device
dev.off()


library(scales)
pdf("/Users/hanxinyu/Desktop/mbimpute/bmean_sd_scatter1.pdf")
par(mfrow = c(3,3))
plot(y = colMeans(y_sim_rec[,-eval_omit]), x = apply(y_sim_rec[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "Complete", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(sim_tab_zi_mbImpute[,-eval_omit]), x = apply(sim_tab_zi_mbImpute[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "No imputation", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot.new()
plot(y = colMeans(alra_mat_eval[,-eval_omit]), x = apply(alra_mat_eval[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "ALRA", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(bhmpf_pre_data[,-eval_omit]), x = apply(bhmpf_pre_data[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "TphPMF", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(saver_mat_eval[,-eval_omit]), x = apply(saver_mat_eval[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "SAVER", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(scImpute_mat_eval[,-eval_omit]), x = apply(scImpute_mat_eval[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "scImpute", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(softImpute_mat_eval[,-eval_omit]), x = apply(softImpute_mat_eval[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "softImpute", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
plot(y = colMeans(mbImpute_mat[,-eval_omit]), x = apply(mbImpute_mat[,-eval_omit], 2, sd), col = alpha("#BC6C36", 0.4), main = "mbImpute", xlim = c(0, 2.5), ylim = c(0,6), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1)
# Add a main title across the top of all plots
mtext("Simulation 2", outer = TRUE, cex = 1.2, line = -1.3)
dev.off()


