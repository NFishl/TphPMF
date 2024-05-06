library(HMP16SData)
V13_stool <-
  V13() %>%
  subset(select = HMP_BODY_SUBSITE == "Stool")

V13_oral <-
  V13() %>%
  subset(select = HMP_BODY_SUBSITE == "Tongue Dorsum")
library(SummarizedExperiment)
oral_16S <- assay(V13_oral, "16SrRNA")
dim(oral_16S)
sum(rowSums(oral_16S) != 0 )

taxa_map <- rowData(V13_oral)

taxa_map[1,7]
rownames(taxa_map)[1]

# Process the data and find the genus that overlaps with the Metagenomic data.

genus_16S <- rownames(oral_16S)

genus_16s_trans <- unlist(lapply(genus_16S, FUN = function(x){
  return(taxa_map[x,7])
}))

oral_16S <- t(oral_16S)
oral_16S_processed <- oral_16S / rowSums(oral_16S) * 10^6
oral_16S_processed <- log10(oral_16S_processed+1.01)
colnames(oral_16S_processed) <- genus_16s_trans
sum(colnames(oral_16S) %in% V13_oral@metadata$phylogeneticTree$tip.label)
edges <- V13_oral@metadata$phylogeneticTree$edge

m = length(V13_oral@metadata$phylogeneticTree$tip.label)
k = 50
nd_mat <- matrix(rep(1, k*m), k, m)
l <- rep(1,k)

n_taxa = m
tree_height = k

#generate the node set
for(i in 1:n_taxa){
  if(i %% 1000 == 0){
    print(i)
  }
  l <- rep(1,tree_height+1)
  l[1] = i
  for(j in 2:(tree_height+1)){
    if(sum(edges[,2] %in% l[j-1]) != 0){
      l[j] = edges[edges[,2] %in% l[j-1], 1]
    }
    else{
      l[j] = NA
    }
  }
  nd_mat[,i] = l[2:(tree_height+1)]
}

saveRDS(nd_mat, "nd_mat.RData")
nd_mat <- readRDS("nd_mat.RData")

library(devtools)
library(HMP16SData)

V13_stool <-
  V13() %>%
  subset(select = HMP_BODY_SUBSITE == "Stool")
library(SummarizedExperiment)
stool_16S <- assay(V13_stool, "16SrRNA")
dim(stool_16S)
sum(rowSums(stool_16S) != 0 )

taxa_map <- rowData(V13_stool)

taxa_map[1,7]
rownames(taxa_map)[1]


# Process the data and find the genus that overlaps with the Metagenomic data.

nz_prop <- apply(stool_16S, 1, FUN = function(x){
  sum(x == 0) / length(x)
})

set.seed(1234)
selected_taxa <- sample(which(nz_prop < 0.7), 300)

stool_16S <- stool_16S[selected_taxa,]

stool_16S <- t(stool_16S)

sum(colnames(stool_16S) %in% V13_stool@metadata$phylogeneticTree$tip.label)

phy_tree_taxa <- which(V13_stool@metadata$phylogeneticTree$tip.label %in% colnames(stool_16S))


nd_mat <- readRDS("nd_mat.RData")

m = length( colnames(stool_16S) )
d1_mat <- matrix(0, nrow = m, ncol = m)

#generate the distance matrix
for(i in 1:m){
  for(j in 1:m){
    int_sc <- intersect(nd_mat[, phy_tree_taxa[i]], nd_mat[,phy_tree_taxa[j]])
    leni <- sum(!is.na(int_sc))
    len1 <- sum(!is.na(nd_mat[,phy_tree_taxa[i]]))
    len2 <- sum(!is.na(nd_mat[,phy_tree_taxa[j]]))
    d1_mat[i, j] = len1 - leni + 1 + len2 - leni + 1
  }
}
diag(d1_mat) = 0

colnames(d1_mat) <- V13_stool@metadata$phylogeneticTree$tip.label[phy_tree_taxa]

data_16S_phy_tree <- d1_mat

stool_16S <- stool_16S[,V13_stool@metadata$phylogeneticTree$tip.label[phy_tree_taxa]]

selected_sample <- which(rowSums(stool_16S) >= 600 & rowSums(stool_16S) < 1500)

stool_16S <- stool_16S[selected_sample,]


meta_data_16S <- V13_stool@colData
#get a reasonable subset of meta features
head(meta_data_16S)
gender <- as.integer(as.factor(meta_data_16S$SEX))
total_reads <- rowSums(stool_16S)
total_reads <- (total_reads - min(total_reads))/sd(total_reads)
meta_data_16S <- data.frame("gender" = gender[selected_sample], "total_reads" = total_reads)




sequencing_depth <- 1000
stool_16S_normalized <- stool_16S / rowSums(stool_16S) * 10^3
stool_16S_processed <- log10(stool_16S_normalized+1.01)

# learn zero introduction parameter k from the benchmark truth.
zp = apply(stool_16S, 2, FUN = function(x){
  return( sum(x == 0)/length(x) )
})
mean_val = apply(stool_16S_processed, 2, FUN = function(x){
  return( mean(x) )
})
zp_map <- cbind(zp, mean_val)
benchmark_truth <- stool_16S_normalized
for(j in 1:dim(stool_16S_normalized)[2]){
  abundance_truth <- log10(stool_16S_normalized[,j] + 1.01)
  # change the zero values to normal distribution
  abundance_truth[abundance_truth < log10(1.02)] = rnorm(sum(abundance_truth < log10(1.02)), mean = mean(abundance_truth[abundance_truth > log10(1.02)]), sd = sd(abundance_truth[abundance_truth > log10(1.02)]) ) 
  benchmark_truth[,j] = floor(10^abundance_truth - 1.01)
  benchmark_truth[benchmark_truth < 0] = 0
}

sum(benchmark_truth == 0) / (dim(benchmark_truth)[1] * dim(benchmark_truth)[2])
benchmark_truth_normalized <- log10( benchmark_truth / rowSums(benchmark_truth) * 10^3 + 1.01)





###TphPMF
dist_obj <- as.dist(data_16S_phy_tree)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 30)
genus_clusters <- cutree(hc, h = 60)
family_clusters <- cutree(hc, h = 90)
microbe_data <- data.frame(
  plant_id = 1:300,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)
microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
hierarchy.info<-microbe_data


# adjust sequencing depth
run_simulation <- function(benchmark_truth, sequencing_depth, meta_data_16S, data_16S_phy_tree, zp_map){
  simulated_truth <- benchmark_truth
  for(i in 1:dim(simulated_truth)[1]){
    simulated_truth[i,] = rmultinom(1, size = sequencing_depth, prob = benchmark_truth[i,] / sum(benchmark_truth[i,]))
  }
  simulated_truth <- floor(simulated_truth)
  
  simulated_truth_normalized <- simulated_truth / rowSums(simulated_truth) * 10^3
  simulated_truth_processed <- log10(simulated_truth_normalized+1.01)
  
  simulated_ZI_data <- simulated_truth_processed
  for(j in 1:dim(simulated_truth_processed)[2]){
    # introduce zeros to the simulated truth processed
    new_mean <- mean(simulated_truth_processed[,j])
    zp_idx <- which( zp_map[,2] + 0.5 > new_mean & zp_map[,2] - 0.5 < new_mean )
    zp_intro <- sample(zp_map[zp_idx,1], 1)
    new_abundance <-  rbinom(length(simulated_truth_processed[,j]), 1, prob = 1 - zp_intro) * simulated_ZI_data[,j]
    new_abundance[new_abundance == 0] = log10(1.01)
    simulated_ZI_data[,j] <- new_abundance
  }
  
  mse_ori <- sum( (simulated_ZI_data  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  #run TphPMF
  simulated_ZI_data[simulated_ZI_data == log10(1.01)] <- NA
  y_sim <- t(simulated_ZI_data)
  trait.info<- y_sim
  setwd("/Users/hanxinyu/Desktop/sim5_10000/")
  tmp.dir<-dirname("/Users/hanxinyu/Desktop/sim5_10000/tmp/")
  GapFilling(trait.info,hierarchy.info,num.samples=30, burn=5, gaps=2,
             num.latent.feats=1,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)
  data <- read.table("/Users/hanxinyu/Desktop/sim5_10000/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat<-t(bhmpf_pre)
  
  mse_imputed <- sum( (imputed_mat  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  mse_benchmark_ori <- sum( (simulated_ZI_data  - benchmark_truth_normalized )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  mse_benchmark_imp <- sum( (imputed_mat  - benchmark_truth_normalized )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  return(c(mse_ori, mse_imputed, mse_benchmark_ori, mse_benchmark_imp))
}

set.seed(1234)
parallel = FALSE

MSE_rec_1000 <- c(0, 0, 0, 0)
for(i in 1:5){
  MSE_rec_1000 = rbind(MSE_rec_1000, run_simulation(benchmark_truth, 1000, meta_data_16S, data_16S_phy_tree, zp_map))
}

MSE_rec_2000 <- c(0, 0, 0, 0)
for(i in 1:5){
  MSE_rec_2000 = rbind(MSE_rec_2000, run_simulation(benchmark_truth, 2000, meta_data_16S, data_16S_phy_tree, zp_map))
}

MSE_rec_5000 <- c(0, 0, 0, 0)
for(i in 1:5){
  MSE_rec_5000 = rbind(MSE_rec_5000, run_simulation(benchmark_truth, 5000, meta_data_16S, data_16S_phy_tree, zp_map))
}

MSE_rec_10000 <- c(0, 0, 0, 0)
for(i in 1:5){
  MSE_rec_10000 = rbind(MSE_rec_10000, run_simulation(benchmark_truth, 10000, meta_data_16S, data_16S_phy_tree, zp_map))
}


MSE_rec_1000 <- MSE_rec_1000[-1,]
print(MSE_rec_1000)
MSE_rec_2000 <- MSE_rec_2000[-1,]
print(MSE_rec_2000)
MSE_rec_5000 <- MSE_rec_5000[-1,]
print(MSE_rec_5000)
MSE_rec_10000 <- MSE_rec_10000[-1,]
print(MSE_rec_10000)



library(ggplot2)

head(ToothGrowth)

df <- cbind( c( MSE_rec_10000[,1], MSE_rec_10000[,2], MSE_rec_5000[,1], MSE_rec_5000[,2], MSE_rec_2000[,1], MSE_rec_2000[,2], MSE_rec_1000[,1], MSE_rec_1000[,2] ),
             rep( c( rep("zi", 5), rep("imp", 5) ), 4),
             c( rep(10000, 10), rep(5000, 10), rep(2000, 10), rep(1000, 10) ) )
colnames(df) <- c("mse", "status", "seq_depth")
df <- data.frame(df)
df$mse <- as.numeric(as.character(df$mse))
df$status <- factor(df$status, levels = c("zi", "imp"))
df$seq_depth <- factor(df$seq_depth, levels = c("1000", "2000", "5000", "10000"))

pdf("/Users/hanxinyu/Desktop/mbimpute/sequencing_depth_adjust2.pdf")
ggplot(df, aes(x=seq_depth, y=mse, fill=status)) +
  geom_boxplot()
dev.off()


####add outlier samples
run_simulation_oulier_samples <- function(benchmark_truth, sequencing_depth, meta_data_16S, data_16S_phy_tree, zp_map){
  simulated_truth <- benchmark_truth
  for(i in 1:dim(simulated_truth)[1]){
    simulated_truth[i,] = rmultinom(1, size = sequencing_depth, prob = benchmark_truth[i,] / sum(benchmark_truth[i,]))
  }
  simulated_truth <- floor(simulated_truth)
  
  simulated_truth_normalized <- simulated_truth / rowSums(simulated_truth) * 10^3
  simulated_truth_processed <- log10(simulated_truth_normalized+1.01)
  
  simulated_ZI_data <- simulated_truth_processed
  for(j in 1:dim(simulated_truth_processed)[2]){
    # introduce zeros to the simulated truth processed
    new_mean <- mean(simulated_truth_processed[,j])
    zp_idx <- which( zp_map[,2] + 0.5 > new_mean & zp_map[,2] - 0.5 < new_mean )
    zp_intro <- sample(zp_map[zp_idx,1], 1)
    new_abundance <-  rbinom(length(simulated_truth_processed[,j]), 1, prob = 1 - zp_intro) * simulated_ZI_data[,j]
    new_abundance[new_abundance == 0] = log10(1.01)
    simulated_ZI_data[,j] <- new_abundance
  }
  
  # create outlier samples
  sim_summary <- apply(simulated_ZI_data, 2, FUN = function(x){
    return( c(max(x), mean(x)) )
  })
  mean_rank <- rank(sim_summary[2,])
  max_rank <- rank(sim_summary[1,])
  nz_number <- apply(simulated_ZI_data, 2, FUN = function(x){
    return( sum(x > log10(1.01)) )
  })
  
  new_sample <- rep(0, 300)
  new_sample[which(mean_rank <= 150 & nz_number > 10)] = sample(sim_summary[1,][which(max_rank > 200)], length(which(mean_rank <= 150 & nz_number > 10)), replace = TRUE)
  new_sample[new_sample == 0] = log10(1.01)
  simulated_ZI_data_outlier <- rbind(simulated_ZI_data, new_sample)
  
  mse_ori <- sum( (simulated_ZI_data  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  
  ##run TphPMF
  simulated_ZI_data[simulated_ZI_data == log10(1.01)] <- NA
  y_sim <- t(simulated_ZI_data)
  trait.info<- y_sim
  setwd("/Users/hanxinyu/Desktop/sim5_outlier/")
  tmp.dir<-dirname("/Users/hanxinyu/Desktop/sim5_outlier/tmp/")
  GapFilling(trait.info,hierarchy.info,num.samples=30, burn=5, gaps=2,
             num.latent.feats=1,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)
  data <- read.table("/Users/hanxinyu/Desktop/sim5_outlier/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat<-t(bhmpf_pre)
  #
  simulated_ZI_data_outlier[simulated_ZI_data_outlier == log10(1.01)] <- NA
  y_sim <- t(simulated_ZI_data_outlier)
  trait.info<- y_sim
  setwd("/Users/hanxinyu/Desktop/sim5_outlier1/")
  tmp.dir<-dirname("/Users/hanxinyu/Desktop/sim5_outlier1/tmp/")
  GapFilling(trait.info,hierarchy.info,num.samples=30, burn=5, gaps=2,
             num.latent.feats=1,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)
  data <- read.table("/Users/hanxinyu/Desktop/sim5_outlier1/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat_outlier1<-t(bhmpf_pre)

  new_sample <- rep(0, 300)
  new_sample[which(mean_rank <= 150 & nz_number > 10)] = sample(sim_summary[1,][which(max_rank > 200)], length(which(mean_rank <= 150 & nz_number > 10)), replace = TRUE)
  new_sample[new_sample == 0] = log10(1.01)
  simulated_ZI_data_outlier2 <- rbind(simulated_ZI_data_outlier, new_sample)
  #
  simulated_ZI_data_outlier2[simulated_ZI_data_outlier2 == log10(1.01)] <- NA
  y_sim <- t(simulated_ZI_data_outlier2)
  trait.info<- y_sim
  setwd("/Users/hanxinyu/Desktop/sim5_outlier2/")
  tmp.dir<-dirname("/Users/hanxinyu/Desktop/sim5_outlier2/tmp/")
  GapFilling(trait.info,hierarchy.info,num.samples=30, burn=5, gaps=2,
             num.latent.feats=1,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)
  data <- read.table("/Users/hanxinyu/Desktop/sim5_outlier2/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat_outlier2<-t(bhmpf_pre)
  
  new_sample <- rep(0, 300)
  new_sample[which(mean_rank <= 150 & nz_number > 10)] = sample(sim_summary[1,][which(max_rank > 200)], length(which(mean_rank <= 150 & nz_number > 10)), replace = TRUE)
  new_sample[new_sample == 0] = log10(1.01)
  simulated_ZI_data_outlier3 <- rbind(simulated_ZI_data_outlier2, new_sample)
  
  #
  simulated_ZI_data_outlier3[simulated_ZI_data_outlier3 == log10(1.01)] <- NA
  y_sim <- t(simulated_ZI_data_outlier3)
  trait.info<- y_sim
  setwd("/Users/hanxinyu/Desktop/sim5_outlier3/")
  tmp.dir<-dirname("/Users/hanxinyu/Desktop/sim5_outlier3/tmp/")
  GapFilling(trait.info,hierarchy.info,num.samples=30, burn=5, gaps=2,
             num.latent.feats=1,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)
  data <- read.table("/Users/hanxinyu/Desktop/sim5_outlier3/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat_outlier3<-t(bhmpf_pre)
  
  
  mse_imputed <- sum( (imputed_mat  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  mse_outlier1 <- sum( (imputed_mat_outlier1[-55,]  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  mse_outlier2 <- sum( (imputed_mat_outlier2[-c(55,56),]  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  mse_outlier3 <- sum( (imputed_mat_outlier3[-c(55,56,57),]  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  
  
  return(c(mse_ori, mse_imputed, mse_outlier1, mse_outlier2, mse_outlier3))
}

MSE_rec_outlier <- c(0, 0, 0, 0)
for(i in 1:10){
  MSE_rec_outlier = rbind(MSE_rec_outlier, run_simulation_oulier_samples(benchmark_truth, 2000, meta_data_16S, data_16S_phy_tree, zp_map))
}

MSE_rec_outlier <- MSE_rec_outlier[-1,]
print(MSE_rec_outlier)

df <- cbind( c( MSE_rec_outlier[,1], MSE_rec_outlier[,2], MSE_rec_outlier[,3], MSE_rec_outlier[,4] ),
             c( rep("no imputation", 10), rep("0 outlier", 10), rep("1 outlier", 10),  rep("2 outliers", 10)))
colnames(df) <- c("mse", "status")
df <- data.frame(df)
df$mse <- as.numeric(as.character(df$mse))
df$status <- factor(df$status, levels = c("no imputation", "0 outlier", "1 outlier", "2 outliers"))

pdf("/Users/hanxinyu/Desktop/mbimpute/outlier_mse1.pdf")
ggplot(df, aes(x=status, y=mse)) +
  geom_boxplot()
dev.off()


run_simulation_oulier_samples <- function(benchmark_truth, sequencing_depth, meta_data_16S, data_16S_phy_tree, zp_map){
  simulated_truth <- benchmark_truth
  for(i in 1:dim(simulated_truth)[1]){
    simulated_truth[i,] = rmultinom(1, size = sequencing_depth, prob = benchmark_truth[i,] / sum(benchmark_truth[i,]))
  }
  simulated_truth <- floor(simulated_truth)
  
  simulated_truth_normalized <- simulated_truth / rowSums(simulated_truth) * 10^3
  simulated_truth_processed <- log10(simulated_truth_normalized+1.01)
  
  simulated_ZI_data <- simulated_truth_processed
  for(j in 1:dim(simulated_truth_processed)[2]){
    # introduce zeros to the simulated truth processed
    new_mean <- mean(simulated_truth_processed[,j])
    zp_idx <- which( zp_map[,2] + 0.5 > new_mean & zp_map[,2] - 0.5 < new_mean )
    zp_intro <- sample(zp_map[zp_idx,1], 1)
    new_abundance <-  rbinom(length(simulated_truth_processed[,j]), 1, prob = 1 - zp_intro) * simulated_ZI_data[,j]
    new_abundance[new_abundance == 0] = log10(1.01)
    simulated_ZI_data[,j] <- new_abundance
  }
  
  # create outlier samples
  sim_summary <- apply(simulated_ZI_data, 2, FUN = function(x){
    return( c(max(x), mean(x)) )
  })
  mean_rank <- rank(sim_summary[2,])
  max_rank <- rank(sim_summary[1,])
  nz_number <- apply(simulated_ZI_data, 2, FUN = function(x){
    return( sum(x > log10(1.01)) )
  })
  
  new_sample <- rep(0, 300)
  new_sample[which(mean_rank <= 150 & nz_number > 10)] = sample(sim_summary[1,][which(max_rank > 200)], length(which(mean_rank <= 150 & nz_number > 10)), replace = TRUE)
  new_sample[new_sample == 0] = log10(1.01)
  simulated_ZI_data_outlier <- rbind(simulated_ZI_data, new_sample)
  
  mse_ori <- sum( (simulated_ZI_data  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  
  ##run TphPMF
  
  data <- read.table("/Users/hanxinyu/Desktop/sim5_outlier/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat<-t(bhmpf_pre)
  #
  
  data <- read.table("/Users/hanxinyu/Desktop/sim5_outlier1/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat_outlier1<-t(bhmpf_pre)
  
  new_sample <- rep(0, 300)
  new_sample[which(mean_rank <= 150 & nz_number > 10)] = sample(sim_summary[1,][which(max_rank > 200)], length(which(mean_rank <= 150 & nz_number > 10)), replace = TRUE)
  new_sample[new_sample == 0] = log10(1.01)
  simulated_ZI_data_outlier2 <- rbind(simulated_ZI_data_outlier, new_sample)
  #
  
  data <- read.table("/Users/hanxinyu/Desktop/sim5_outlier2/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat_outlier2<-t(bhmpf_pre)
  
  new_sample <- rep(0, 300)
  new_sample[which(mean_rank <= 150 & nz_number > 10)] = sample(sim_summary[1,][which(max_rank > 200)], length(which(mean_rank <= 150 & nz_number > 10)), replace = TRUE)
  new_sample[new_sample == 0] = log10(1.01)
  simulated_ZI_data_outlier3 <- rbind(simulated_ZI_data_outlier2, new_sample)
  
  #
  
  data <- read.table("/Users/hanxinyu/Desktop/sim5_outlier3/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat_outlier3<-t(bhmpf_pre)
  
  
  mse_imputed <- sum( (imputed_mat  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  mse_outlier1 <- sum( (imputed_mat_outlier1[-55,]  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  mse_outlier2 <- sum( (imputed_mat_outlier2[-c(55,56),]  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  mse_outlier3 <- sum( (imputed_mat_outlier3[-c(55,56,57),]  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  
  
  return(list(
    simulated_truth_processed = simulated_truth_processed,
    simulated_ZI_data = simulated_ZI_data,
    imputed_mat = imputed_mat,
    simulated_ZI_data_outlier = simulated_ZI_data_outlier,
    imputed_mat_outlier1 = imputed_mat_outlier1,
    simulated_ZI_data_outlier2 = simulated_ZI_data_outlier2,
    imputed_mat_outlier2 = imputed_mat_outlier2,
    simulated_ZI_data_outlier3 = simulated_ZI_data_outlier3,
    imputed_mat_outlier3 = imputed_mat_outlier3
  ))
}

results <- run_simulation_oulier_samples(benchmark_truth, 2000, meta_data_16S, data_16S_phy_tree, zp_map)
simulated_truth_processed <- results$simulated_truth_processed
simulated_ZI_data <- results$simulated_ZI_data
imputed_mat <- results$imputed_mat
simulated_ZI_data_outlier <- results$simulated_ZI_data_outlier
imputed_mat_outlier1 <- results$imputed_mat_outlier1
simulated_ZI_data_outlier2 <- results$simulated_ZI_data_outlier2
imputed_mat_outlier2 <- results$imputed_mat_outlier2
simulated_ZI_data_outlier3 <- results$simulated_ZI_data_outlier3
imputed_mat_outlier3 <- results$imputed_mat_outlier3


i = 38
pdf("/Users/hanxinyu/Desktop/mbimpute/outlier_example1.pdf")
par(mfrow = c(2, 3))
hist(simulated_ZI_data[,i], xlim = c(0,2), breaks = seq(0, 2, 0.1))
hist( simulated_ZI_data_outlier[,i], xlim = c(0,2), breaks = seq(0, 2, 0.1))
hist( simulated_ZI_data_outlier2[,i], xlim = c(0,2), breaks = seq(0, 2, 0.1))

hist(imputed_mat[,i], xlim = c(-3,2), breaks = seq(-3, 2, 0.1))
hist(imputed_mat_outlier1[,i], xlim = c(-3,2), breaks = seq(-3, 2, 0.1))
hist(imputed_mat_outlier2[,i], xlim = c(-3,2), breaks = seq(-3, 2, 0.1))
dev.off()



#########Sim_seq_depth_imp_methods
##mbimpute
library(Matrix)
library(glmnet)
data_fit2 <- function(y_sim, x, D, k){
  #loading used functions
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
  # design_mat_row_gen2(y_sim, x[1:n], confidence_set[i,1], confidence_set[i,2], close_taxa)
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
      # row_gen <- rep(0, m*k + (n-1) * n + n*p)
      # 
      # row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
      # row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
      # row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i])
      
      non_zero_idx <- c(((j-1)*k+1):(j*k), (m*k + (i-1)*(n-1)+1):(k*m + i*(n-1)),(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p) )
      non_zero_values <- c( as.numeric(count_mat[i,close_taxa[[j]]]), as.numeric(count_mat[-i,j]), as.numeric(covariate_mat[i]))
    }else{
      p = dim(covariate_mat)[2]
      i = row_index
      j = col_index
      
      close_taxa_set <- close_taxa[[j]]
      #generate a row including response and a row of design matrix.
      # row_gen <- rep(0, m*k + (n-1) * n + n*p)
      # 
      # row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
      # row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
      # row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i,])
      
      non_zero_idx <- c(((j-1)*k+1):(j*k), (m*k + (i-1)*(n-1)+1):(k*m + i*(n-1)),(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p) )
      non_zero_values <- c( as.numeric(count_mat[i,close_taxa[[j]]]), as.numeric(count_mat[-i,j]), as.numeric(covariate_mat[i,]))
    }
    return(list("nz_idx" = non_zero_idx, "nz_val" = non_zero_values))
  }
  design_mat_row_gen2_imp <- function(count_mat, covariate_mat, row_index, col_index, close_taxa){
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
      row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i])
    }else{
      p = dim(covariate_mat)[2]
      i = row_index
      j = col_index
      
      close_taxa_set <- close_taxa[[j]]
      #generate a row including response and a row of design matrix.
      row_gen <- rep(0, m*k + (n-1) * n + n*p)
      
      row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
      row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
      row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i,])
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
  
  #keep the original matrix
  y_imp <- y_sim
  # print(dim(y_imp))
  # #keep the vectors that we neither impute nor borrow information
  # remain_vec <- which(unlist(filter) == 0)
  # y_rem <- y_sim[, remain_vec]
  
  #perform imputation on the rest
  filter_vec <- which(unlist(filter) == 1)
  y_sim = y_sim[, filter_vec]
  D = D[filter_vec,filter_vec]
  
  #apply the imputation method on the simulated matrix
  m = dim(y_sim)[2]
  n = dim(y_sim)[1]
  
  #x[1:50,]
  X <- cbind(rep(1,n), x[1:n,])
  X <- as.matrix(X)
  
  #identifying the group needs imputation and group doesn't need imputation
  
  result <- lapply(X = 1:ncol(y_sim), FUN = function(col_i){
    y = y_sim[,col_i]
    return(gamma_norm_mix(y, X)$d)
  })
  
  hd_result <- result[1:6]
  confidence_set <- c(0,0)
  impute_set <- c(0,0)
  for(i in 1:length(result)){
    d_set <- result[[i]]
    cf_set <- which(d_set < 0.5)
    im_set <- which(d_set > 0.5)
    n_c = length(cf_set)
    confidence_set <- rbind(confidence_set, cbind(cf_set, rep(i, n_c)))
    impute_set <- rbind(impute_set, cbind(im_set, rep(i, n-n_c)))
  }
  #print(dim(confidence_set)[1]/(dim(confidence_set)[1] + dim(impute_set)[1]))
  
  #find k closest taxa
  # k = 10
  #generate close set
  close_taxa <- list()
  for(j in 1:m){
    #close_taxa[[j-1]] = which(D[,j] %in% sort(unique(D[,j]))[2:(k+1)])
    close_dist <- D[D[,j] %in% sort(D[,j])[2:(k+1)],j]
    close_taxa_vec = which(D[,j] %in% close_dist[close_dist != max(close_dist)])
    if(length(close_taxa_vec) < k){
      close_taxa_vec <- c(close_taxa_vec, which(D[,j] == max(close_dist))[1:(k-length(close_taxa_vec))])
    }
    close_taxa[[j]] = close_taxa_vec
  }
  
  
  #generate design matrix
  p = dim(X)[2]
  if(is.null(p)){
    p = 1
  }
  #generate a row including response and a row of design matrix.
  row_length <- m * k + (n-1) * n + n*p
  
  size = 5000
  if(dim(confidence_set)[1] - 1 < size){
    mat_num = ceiling( (dim(confidence_set)[1] - 1) / size)
    remaining = (dim(confidence_set)[1] - 1) %% size
    design_mat_fit =  sparseMatrix(i = 1, j =1, x = 0, dims = c(remaining, row_length))
    track = ( (mat_num-1)*size+1):((mat_num - 1)*size+remaining)
    for(i in 1:remaining){
      if(is.vector(X)){
        result <- design_mat_row_gen2(y_sim, X[1:n], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
        design_mat_fit[i,result$nz_idx] <- result$nz_val
      }
      else{
        result <- design_mat_row_gen2(y_sim, X[1:n,], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
        design_mat_fit[i,result$nz_idx] <- result$nz_val
      }
    }
  }else{
    mat_num = ceiling( (dim(confidence_set)[1] - 1) / size)
    mat_list = list()
    if(!parallel){
      for(mat_new in 1:(mat_num-1)){
        print(mat_num-1-mat_new)
        design_mat_fit = sparseMatrix(i = 1, j =1, x = 0, dims = c(size, row_length))
        track = ((mat_new-1)*size+1):(mat_new*size)
        for(i in 1:size){
          if(is.vector(X)){
            result <- design_mat_row_gen2(y_sim, X[1:n], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
            design_mat_fit[i,result$nz_idx] <- result$nz_val
          }
          else{
            result <- design_mat_row_gen2(y_sim, X[1:n,], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
            design_mat_fit[i,result$nz_idx] <- result$nz_val
          }
        }
        mat_list[[mat_new]] = design_mat_fit
      }
    }else{
      no_cores <- max(ncores, detectCores() - 1)  
      registerDoParallel(cores=no_cores)  
      cl <- makeCluster(no_cores, "fork")  
      f <- function(mat_new){
        size = 2000
        design_mat_fit = sparseMatrix(i = 1, j =1, x = 0, dims = c(size, row_length))
        track = ((mat_new-1)*size+1):(mat_new*size)
        for(i in 1:size){
          if(is.vector(X)){
            result <- design_mat_row_gen2(y_sim, X[1:n], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
            design_mat_fit[i,result$nz_idx] <- result$nz_val
          }
          else{
            result <- design_mat_row_gen2(y_sim, X[1:n,], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
            design_mat_fit[i,result$nz_idx] <- result$nz_val
          }
        }
        return(design_mat_fit)
      }
      result <- parLapply(cl, 1:(mat_num-1), f)
    }
    remaining = (dim(confidence_set)[1] - 1) %% size
    design_mat_fit =  sparseMatrix(i = 1, j =1, x = 0, dims = c(remaining, row_length))
    track = ( (mat_num-1)*size+1):((mat_num - 1)*size+remaining)
    for(i in 1:remaining){
      if(is.vector(X)){
        result <- design_mat_row_gen2(y_sim, X[1:n], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
        design_mat_fit[i,result$nz_idx] <- result$nz_val
      }
      else{
        result <- design_mat_row_gen2(y_sim, X[1:n,], confidence_set[track[i]+1,1], confidence_set[track[i]+1,2], close_taxa)
        design_mat_fit[i,result$nz_idx] <- result$nz_val
      }
    }
    for(j in length(mat_list):1){
      design_mat_fit = rbind(mat_list[[j]], design_mat_fit)
    }
  }
  
  # y_sim <- as.matrix(y_sim@.Data)
  #generate covariate matrix for values need imputation
  # design_mat_impute <- sparseMatrix(i = 1, j =1, x = 0, dims = c(dim(impute_set)[1] - 1, row_length))
  # for(i in 2:dim(impute_set)[1]){
  #   if(i %% 10000 == 0){
  #     print(i)
  #   }
  #   if(is.vector(X)){
  #     result <- design_mat_row_gen2(y_sim, X[1:n], impute_set[i,1], impute_set[i,2], close_taxa)
  #     design_mat_impute[i-1,result$nz_idx] <- result$nz_val
  #   }
  #   else{
  #     result <- design_mat_row_gen2(y_sim, X[1:n,], impute_set[i,1], impute_set[i,2], close_taxa)
  #     design_mat_impute[i-1,result$nz_idx] <- result$nz_val
  #   }
  # }
  
  mse <- rep(5,4)
  psi <- c(0.1, 0.5, 1, 2)
  c1 <- list()
  for(i in 1:4){
    weights_pen <- D^psi[i]
    penalized_weights <- c()
    for(j in 1:m){
      penalized_weights <- c(penalized_weights, weights_pen[close_taxa[[j]],j])
    }
    #penalized_weights <- c(penalized_weights, rep(1, n*(n-1)), rep(c(0, rep(1, p-1)), n))
    penalized_weights <- c(penalized_weights, rep(1, n*(n-1)), rep(1, n*p))
    #print(length(penalized_weights))
    #print(dim(design_mat_fit))
    response <- y_sim[confidence_set]
    set.seed(1)
    print(dim(design_mat_fit))
    cv.result <- cv.glmnet(x = design_mat_fit, y = response, family = "gaussian", penalty.factor = penalized_weights, intercept = TRUE)
    c1[[i]] <- coef(cv.result, s = cv.result$lambda.min)
    mse[i] = min(cv.result$cvm)
    print(mse)
  }
  c1 <- c1[[which.min(mse)]]
  impute_mat = y_sim
  for(i in 2:dim(impute_set)[1]){
    if(i %% 10000 == 0){
      print(i)
    }
    if(is.vector(X)){
      result <- design_mat_row_gen2_imp(y_sim, X[1:n], impute_set[i,1], impute_set[i,2], close_taxa)
      impute_mat[impute_set[i,1], impute_set[i,2]] = min(max(y_sim[,impute_set[i,2]]), max(c(1, result) %*% c1, log10(1.01)))
    }else{
      result <- design_mat_row_gen2_imp(y_sim, X[1:n,], impute_set[i,1], impute_set[i,2], close_taxa)
      impute_mat[impute_set[i,1], impute_set[i,2]] <- min(max(y_sim[,impute_set[i,2]]), max(c(1, result) %*% c1, log10(1.01)))
    }
  }
  
  y_imp[, filter_vec] = impute_mat
  return(y_imp)
}


library(devtools)
library(scImpute)
library(ggplot2)
library(gplots)
library(SAVER)
library(softImpute)
library(reticulate)

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
  # Computes the k-rank approximation to A_norm and adjusts it according to the
  # error distribution learned from the negative values.
  #
  # Args:
  #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
  #   k : the rank of the rank-k approximation. Set to 0 for automated choice of k.
  #   q : the number of additional power iterations in randomized SVD
  #   use.mkl : Use the Intel MKL based implementation of SVD. Needs to be
  #             installed from https://github.com/KlugerLab/rpca-mkl.
  #   mkl.seed : Only relevant if use.mkl=T. Set the seed for the random
  #   generator for the Intel MKL implementation of SVD. Any number <0 will
  #   use the current timestamp. If use.mkl=F, set the seed using
  #   set.seed() function as usual.
  #
  # Returns:
  #   A list with three items
  #       1) The rank k approximation of A_norm.
  #       2) The rank k approximation of A_norm, adaptively thresholded
  #       3) The rank k approximation of A_norm, adaptively thresholded and
  #       with the first two moments of the non-zero values matched to the
  #       first two moments of the non-zeros of A_norm. This is the completed
  #       matrix most people will want to work with
  # Example:
  #     result.completed <- adjusted_svd(A_norm,15)
  #     A_norm_rank15 <- result.completed[[1]]     # The low rank approximation for reference purposes...not suggested for matrix completion
  #     A_norm_rank15_cor <- result.completed[[3]] # The actual adjusted, completed matrix
  
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



stool_16S <- read.csv("/Users/hanxinyu/Desktop/stool_16S.csv", row.names = 1)
meta_data_16S <- read.csv("/Users/hanxinyu/Desktop/meta_data_16S.csv", row.names = 1)
data_16S_phy_tree <- read.csv("/Users/hanxinyu/Desktop/data_16S_phy_tree.csv", row.names = 1)


library(usethis)
library(devtools)
library(Matrix)
library(BHPMF)

sequencing_depth <- 1000
stool_16S_normalized <- stool_16S / rowSums(stool_16S) * 10^3
stool_16S_processed <- log10(stool_16S_normalized+1.01)

# learn zero introduction parameter k from the benchmark truth.
zp = apply(stool_16S, 2, FUN = function(x){
  return( sum(x == 0)/length(x) )
})
mean_val = apply(stool_16S_processed, 2, FUN = function(x){
  return( mean(x) )
})
zp_map <- cbind(zp, mean_val)

benchmark_truth <- stool_16S_normalized
for(j in 1:dim(stool_16S_normalized)[2]){
  abundance_truth <- log10(stool_16S_normalized[,j] + 1.01)
  # change the zero values to normal distribution
  abundance_truth[abundance_truth < log10(1.02)] = rnorm(sum(abundance_truth < log10(1.02)), mean = mean(abundance_truth[abundance_truth > log10(1.02)]), sd = sd(abundance_truth[abundance_truth > log10(1.02)]) ) 
  benchmark_truth[,j] = floor(10^abundance_truth - 1.01)
  benchmark_truth[benchmark_truth < 0] = 0
}

sum(benchmark_truth == 0) / (dim(benchmark_truth)[1] * dim(benchmark_truth)[2])
benchmark_truth_normalized <- log10( benchmark_truth / rowSums(benchmark_truth) * 10^3 + 1.01)



# adjust sequencing depth
run_simulation <- function(benchmark_truth, sequencing_depth, meta_data_16S, data_16S_phy_tree, zp_map){
  simulated_truth <- benchmark_truth
  for(i in 1:dim(simulated_truth)[1]){
    simulated_truth[i,] = rmultinom(1, size = sequencing_depth, prob = benchmark_truth[i,] / sum(benchmark_truth[i,]))
  }
  simulated_truth <- floor(simulated_truth)
  
  simulated_truth_normalized <- simulated_truth / rowSums(simulated_truth) * 10^3
  simulated_truth_processed <- log10(simulated_truth_normalized+1.01)
  
  simulated_ZI_data <- simulated_truth_processed
  for(j in 1:dim(simulated_truth_processed)[2]){
    # introduce zeros to the simulated truth processed
    new_mean <- mean(simulated_truth_processed[,j])
    zp_idx <- which( zp_map[,2] + 0.5 > new_mean & zp_map[,2] - 0.5 < new_mean )
    zp_intro <- sample(zp_map[zp_idx,1], 1)
    new_abundance <-  rbinom(length(simulated_truth_processed[,j]), 1, prob = 1 - zp_intro) * simulated_ZI_data[,j]
    new_abundance[new_abundance == 0] = log10(1.01)
    simulated_ZI_data[,j] <- new_abundance
  }
  
  mse_ori <- sum( (simulated_ZI_data  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  
  # run mbImpute
  D = data_16S_phy_tree
  x = meta_data_16S
  k = 5
  y_sim = simulated_ZI_data
  
  imputed_mat = data_fit2(y_sim, x, D, k)
  mse_imputed <- sum( (imputed_mat  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  mse_benchmark_ori <- sum( (simulated_ZI_data  - benchmark_truth_normalized )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  mse_benchmark_imp <- sum( (imputed_mat  - benchmark_truth_normalized )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  registerDoSEQ()
  Saver_imp <- saver(t(y_sim), ncores = 1, estimates.only = TRUE)
  # sum( (t(Saver_imp)/rescale_factor  - benchmark_truth_normalized )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  Saver_mse <- sum( (t(Saver_imp)  - simulated_truth_processed )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  write.csv(t(y_sim), "y_sim_trans.csv")
  scimpute(count_path = "y_sim_trans.csv", Kcluster = 1, out_dir = "y_sim_imp", ncores = 1)
  scImpute_imp <- read.csv("y_sim_impscimpute_count.csv")
  scImpute_imp <- scImpute_imp[,-1]
  # sum( (t(scImpute_imp)/rescale_factor  - benchmark_truth_normalized )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  scImpute_mse <- sum( (t(scImpute_imp)  - simulated_truth_processed )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  y_sim_norm <- normalize_data(10^y_sim - 1.01)
  k_chosen <- choose_k(y_sim_norm, K = 49, noise_start = 44)$k
  y_sim_norm <- as.matrix(y_sim_norm)
  y_sim_norm <- alra(y_sim_norm, k = 2)$A_norm_rank_k_cor_sc
  y_sim_norm <- log10(y_sim_norm + 1.01)
  ALRA_mse <- sum( ( y_sim_norm  - simulated_truth_processed )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  # ALRA_mse <- sum( (y_sim_norm/rescale_factor  - benchmark_truth_normalized )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])

  y_sim_count = 10^y_sim-1.01
  y_sim_count[y_sim_count == 0] <- NA
  y_sim_count <- as.matrix(y_sim_count)
  print(y_sim_count)
  softImpute_y_imp <- softImpute(y_sim_count)
  softImpute_y_imp <- complete(y_sim_count, softImpute_y_imp)
  softImpute_y_imp[softImpute_y_imp <= 0] = 0
  softImpute_mse <- sum( ( log10(softImpute_y_imp  + 1.01 )  - simulated_truth_processed )^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  data <- read.table("/Users/hanxinyu/Desktop/sim5_BHPMF_10000/mean_gap_filled.txt", header = TRUE, sep = "\t")
  bhmpf_pre <- as.matrix(data)
  imputed_mat<-t(bhmpf_pre)
  BHPMF_mse <- sum( (imputed_mat  - simulated_truth_processed)^2 ) / (dim(simulated_ZI_data)[1] * dim(simulated_ZI_data)[2])
  
  return(c(mse_ori, mse_imputed, mse_benchmark_ori, mse_benchmark_imp,  Saver_mse, scImpute_mse, ALRA_mse, softImpute_mse,BHPMF_mse))
}



set.seed(1234)
parallel = FALSE
MSE_rec_1000 <- c(0, 0, 0, 0, 0, 0, 0)
for(i in 1:5){
  MSE_rec_1000 = rbind(MSE_rec_1000, run_simulation(benchmark_truth, 1000, meta_data_16S, data_16S_phy_tree, zp_map))
}
print(MSE_rec_1000)

set.seed(1234)
parallel = FALSE
MSE_rec_2000 <- c(0, 0, 0, 0, 0, 0, 0)
for(i in 1:5){
  MSE_rec_2000 = rbind(MSE_rec_2000, run_simulation(benchmark_truth, 2000, meta_data_16S, data_16S_phy_tree, zp_map))
}
print(MSE_rec_2000)


set.seed(1234)
parallel = FALSE
MSE_rec_5000 <- c(0, 0, 0, 0, 0, 0, 0)
for(i in 1:5){
  MSE_rec_5000 = rbind(MSE_rec_5000, run_simulation(benchmark_truth, 5000, meta_data_16S, data_16S_phy_tree, zp_map))
}
print(MSE_rec_5000)


set.seed(1234)
parallel = FALSE
MSE_rec_10000 <- c(0, 0, 0, 0, 0, 0, 0)
for(i in 1:5){
  MSE_rec_10000 = rbind(MSE_rec_10000, run_simulation(benchmark_truth, 10000, meta_data_16S, data_16S_phy_tree, zp_map))
}
print(MSE_rec_10000)


df <- cbind( c( mean(MSE_rec_10000[,1]), mean(MSE_rec_10000[,2]), mean(MSE_rec_10000[,5]), mean(MSE_rec_10000[,6]), mean(MSE_rec_10000[,7]), mean(MSE_rec_10000[,8]), mean(MSE_rec_10000[,9]), mean(MSE_rec_5000[,1]), mean(MSE_rec_5000[,2]), mean(MSE_rec_5000[,5]), mean(MSE_rec_5000[,6]), mean(MSE_rec_5000[,7]), mean(MSE_rec_5000[,8]), mean(MSE_rec_5000[,9]), mean(MSE_rec_2000[,1]), mean(MSE_rec_2000[,2]), mean(MSE_rec_2000[,5]), mean(MSE_rec_2000[,6]), mean(MSE_rec_2000[,7]), mean(MSE_rec_2000[,8]), mean(MSE_rec_2000[,9]), mean(MSE_rec_1000[,1]), mean(MSE_rec_1000[,2]), mean(MSE_rec_1000[,5]), mean(MSE_rec_1000[,6]), mean(MSE_rec_1000[,7]), mean(MSE_rec_1000[,8]), mean(MSE_rec_1000[,9])),
             c( sd(MSE_rec_10000[,1]), sd(MSE_rec_10000[,2]), sd(MSE_rec_10000[,5]), sd(MSE_rec_10000[,6]), sd(MSE_rec_10000[,7]), sd(MSE_rec_10000[,8]), sd(MSE_rec_10000[,9]), sd(MSE_rec_5000[,1]), sd(MSE_rec_5000[,2]), sd(MSE_rec_5000[,5]), sd(MSE_rec_5000[,6]), sd(MSE_rec_5000[,7]), sd(MSE_rec_5000[,8]), sd(MSE_rec_5000[,9]), sd(MSE_rec_2000[,1]), sd(MSE_rec_2000[,2]), sd(MSE_rec_2000[,5]), sd(MSE_rec_2000[,6]), sd(MSE_rec_2000[,7]), sd(MSE_rec_2000[,8]), sd(MSE_rec_2000[,9]), sd(MSE_rec_1000[,1]), sd(MSE_rec_1000[,2]), sd(MSE_rec_1000[,5]), sd(MSE_rec_1000[,6]), sd(MSE_rec_1000[,7]), sd(MSE_rec_1000[,8]), sd(MSE_rec_1000[,9])),
             
             rep( c( "No imputation", "mbImpute",  "SAVER", "scImpute", "ALRA", "SoftImpute","TphPMF"), 4),
             c( rep(10000, 7), rep(5000, 7), rep(2000, 7), rep(1000, 7) ) )
colnames(df) <- c("MSE_mean", "MSE_sd", "Method", "Depth")

df <- as.data.frame(df)
df$Method <- factor(df$Method)
df$Depth <- factor(df$Depth, levels = c("1000", "2000", "5000", "10000"))
df$MSE_mean <- as.numeric(as.character(df$MSE_mean))
df$MSE_sd <- as.numeric(as.character(df$MSE_sd))

gp1 <- ggplot(df, aes(x = Depth, y = MSE_mean, fill = Method)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) + 
  geom_errorbar(aes(ymin = MSE_mean - MSE_sd, ymax = MSE_mean + MSE_sd), 
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Sequencing Depth", y = "Mean Squared Error (MSE)", fill = "Imputation Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Paired") 
ggsave("/Users/hanxinyu/Desktop/mbimpute/sequencing_depth_adjust_Multiple_imp4_1.pdf", plot = gp1, width = 20, height = 10)


