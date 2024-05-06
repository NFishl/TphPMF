library(curatedMetagenomicData)
###T2D data
##read in Qin data
curatedMetagenomicData("QinJ_2012.relative_abundance", dryrun = FALSE, counts = TRUE, rownames = "short")

Qin_cd <- curatedMetagenomicData("QinJ_2012.relative_abundance", dryrun = FALSE, counts = TRUE)
str(Qin_cd)
tree_sum_exp <- Qin_cd[["2021-10-14.QinJ_2012.relative_abundance"]]
assay_data <- assays(tree_sum_exp)
relative_abundance_matrix <- assay_data[["relative_abundance"]]
Qin_OTU <- relative_abundance_matrix

col_data <- colData(tree_sum_exp)
str(col_data)
all_list_data <- col_data@listData
Qin_meta <- data.frame(all_list_data)
phylo_list <- rowTree(tree_sum_exp)
edges <- phylo_list$edge

row_ranges <- rowRanges(tree_sum_exp)
element_metadata <- row_ranges@elementMetadata
list_data <- element_metadata@listData
species <- list_data$species
genus <- list_data$genus
family <- list_data$family
Qin_hierarchy.info <- data.frame(
  plant_id = 1:length(species),
  species = species,
  genus = genus,
  family = family
)


#
Qin_OTU <- t(Qin_OTU)
na_idx <- which( is.na(Qin_meta$study_condition))
Qin_meta <- Qin_meta[-na_idx,]
head(Qin_meta)
chosen_meta <- Qin_meta[,c("study_condition", "BMI", "number_reads")]
chosen_meta$BMI <- (chosen_meta$BMI-min(chosen_meta$BMI))/sqrt(var(chosen_meta$BMI))
chosen_meta$number_reads <- (chosen_meta$number_reads-min(chosen_meta$number_reads))/sqrt(var(chosen_meta$number_reads))

Qin_OTU <- Qin_OTU[-na_idx,]
OTU_tab <- Qin_OTU
Qin_condition <- chosen_meta$study_condition


m = dim(OTU_tab)[2]
n = dim(OTU_tab)[1]
k = 50
nd_mat <- matrix(rep(1, k*m), k, m)
l <- rep(1,k)

n_taxa = m
tree_height = k

#generate the node set
for(i in 1:n_taxa){
  print(i)
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


d1_mat <- matrix(0, nrow = m, ncol = m)
taxa_vec <- match(1:m, edges[,2])
#generate the distance matrix
for(i in 1:m){
  for(j in 1:m){
    int_sc <- intersect(nd_mat[,i], nd_mat[,j])
    leni <- sum(!is.na(int_sc))
    len1 <- sum(!is.na(nd_mat[,i]))
    len2 <- sum(!is.na(nd_mat[,j]))
    d1_mat[i, j] = len1 - leni + 1 + len2 - leni + 1
  }
}
diag(d1_mat) = 0

#process the otu_table
otu_tab_processed <- OTU_tab / rowSums(OTU_tab) * 10^6
otu_tab_processed <- log10(otu_tab_processed+1.01)

colnames(d1_mat) <- colnames(otu_tab_processed)
rownames(d1_mat) <- colnames(otu_tab_processed)

#construct x
head(chosen_meta)
x <- cbind(rep(1, dim(chosen_meta)[1]), chosen_meta)

Qin_otable <- otu_tab_processed
Qin_meta_data <- x
Qin_D_mat <- d1_mat

# filter out the taxa that have too few reads in one of the two conditions.
zero_taxa <- c()
for(j in 1:dim(otu_tab_processed)[2]){
  for(cond in unique(Qin_condition)){
    cond_samp <- otu_tab_processed[Qin_condition == cond,j]
    if( sum(cond_samp != log10(1.01)) / length(cond_samp) < 0.1 ){
      zero_taxa <- c(zero_taxa, j)
    }
  }
}
zero_taxa <- unique(zero_taxa)

Qin_otu_tab = Qin_otable[,-zero_taxa]
Qin_raw <- Qin_otu_tab
Qin_D_mat = Qin_D_mat[-zero_taxa , -zero_taxa ]
Qin_hierarchy.info = Qin_hierarchy.info[-zero_taxa , ]


##read in Karlsson data
Karlsson_cd <- curatedMetagenomicData("KarlssonFH_2013.relative_abundance", dryrun = FALSE, counts = TRUE)
tree_sum_exp <- Karlsson_cd[["2021-10-14.KarlssonFH_2013.relative_abundance"]]
assay_data <- assays(tree_sum_exp)
relative_abundance_matrix <- assay_data[["relative_abundance"]]
Karlsson_OTU <- relative_abundance_matrix
col_data <- colData(tree_sum_exp)
all_list_data <- col_data@listData
Karlsson_meta <- data.frame(all_list_data)
phylo_list <- rowTree(tree_sum_exp)
edges <- phylo_list$edge
row_ranges <- rowRanges(tree_sum_exp)
element_metadata <- row_ranges@elementMetadata
list_data <- element_metadata@listData
species <- list_data$species
genus <- list_data$genus
family <- list_data$family
Karlsson_hierarchy.info <- data.frame(
  plant_id = 1:length(species),
  species = species,
  genus = genus,
  family = family
)



Karlsson_OTU <- t(Karlsson_OTU)



#get a reasonable subset of meta features
head(Karlsson_meta)
########## CHANGE #############
apply(Karlsson_meta, 2, function(x){
  x <- as.numeric(x)
  var(x[!is.na(x)])/mean(x[!is.na(x)])
})
chosen_meta <- Karlsson_meta[,c("study_condition", "age", "number_reads", "triglycerides", "hba1c", "ldl", "c_peptide", "cholesterol", "glucose", "adiponectin", "hscrp", "leptin")]
chosen_meta$age <- (chosen_meta$age - min(chosen_meta$age)) / sqrt(var(chosen_meta$age))
chosen_meta$number_reads <- (chosen_meta$number_reads-min(chosen_meta$number_reads))/sqrt(var(chosen_meta$number_reads))
chosen_meta$triglycerides <- (chosen_meta$triglycerides - min(chosen_meta$triglycerides))/sqrt(var(chosen_meta$triglycerides))
chosen_meta$hba1c <- (chosen_meta$hba1c - min(chosen_meta$hba1c)) / sqrt(var(chosen_meta$hba1c))
chosen_meta$ldl <- (chosen_meta$ldl  - min(chosen_meta$ldl)) / sqrt(var(chosen_meta$ldl)) 
chosen_meta$c_peptide <- (chosen_meta$c_peptide - min(chosen_meta$c_peptide)) / sqrt(var(chosen_meta$c_peptide))
chosen_meta$cholesterol <- (chosen_meta$cholesterol - min(chosen_meta$cholesterol)) / sqrt(var(chosen_meta$cholesterol)) 
chosen_meta$glucose[is.na(chosen_meta$glucose)] = mean(chosen_meta$glucose[!is.na(chosen_meta$glucose)])
chosen_meta$glucose <- (chosen_meta$glucose - min(chosen_meta$glucose)) / sqrt(var(chosen_meta$glucose))
chosen_meta$adiponectin[is.na(chosen_meta$adiponectin)] <- mean(chosen_meta$adiponectin[!is.na(chosen_meta$adiponectin)])
chosen_meta$adiponectin <- (chosen_meta$adiponectin - min(chosen_meta$adiponectin)) / sqrt(var(chosen_meta$adiponectin))
chosen_meta$hscrp <- (chosen_meta$hscrp - min(chosen_meta$hscrp)) / sqrt( var(chosen_meta$hscrp))
chosen_meta$leptin[is.na(chosen_meta$leptin)] <- mean(chosen_meta$leptin[!is.na(chosen_meta$leptin)])
chosen_meta$leptin <- (chosen_meta$leptin - min(chosen_meta$leptin)) / sqrt( var(chosen_meta$leptin))

OTU_tab <- Karlsson_OTU

#get the distance matrix D base on phylogenetic tree
m = dim(OTU_tab)[2]
n = dim(OTU_tab)[1]
k = 50
nd_mat <- matrix(rep(1, k*m), k, m)
l <- rep(1,k)

n_taxa = m
tree_height = k

#generate the node set
for(i in 1:n_taxa){
  print(i)
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


d1_mat <- matrix(0, nrow = m, ncol = m)
taxa_vec <- match(1:m, edges[,2])
#generate the distance matrix
for(i in 1:m){
  for(j in 1:m){
    int_sc <- intersect(nd_mat[,i], nd_mat[,j])
    leni <- sum(!is.na(int_sc))
    len1 <- sum(!is.na(nd_mat[,i]))
    len2 <- sum(!is.na(nd_mat[,j]))
    d1_mat[i, j] = len1 - leni + 1 + len2 - leni + 1
  }
}
diag(d1_mat) = 0

#process the otu_table
otu_tab_processed <- OTU_tab / rowSums(OTU_tab) * 10^6
otu_tab_processed <- log10(otu_tab_processed+1.01)

#construct x
head(chosen_meta)
x <- cbind(rep(1, dim(chosen_meta)[1]), chosen_meta)

Karlsson_otable <- otu_tab_processed
Karlsson_meta_data <- x
Karlsson_D_mat <- d1_mat



Karlsson_condition = Karlsson_meta_data$study_condition

Karlsson_otable <- Karlsson_otable[Karlsson_condition != "IGT",]
Karlsson_meta_data <- Karlsson_meta_data[Karlsson_condition != "IGT",]
Karlsson_condition = Karlsson_condition [Karlsson_condition != "IGT"]

# filter out the taxa that have too few reads in one of the two conditions.
zero_taxa <- c()
for(j in 1:dim(Karlsson_otable)[2]){
  for(cond in unique(Karlsson_condition)){
    cond_samp <- Karlsson_otable[Karlsson_condition == cond,j]
    if( sum(cond_samp != log10(1.01)) / length(cond_samp) < 0.1 ){
      zero_taxa <- c(zero_taxa, j)
    }
  }
}
zero_taxa <- unique(zero_taxa)

otu_tab = Karlsson_otable[,-zero_taxa]
D = Karlsson_D_mat[-zero_taxa , -zero_taxa ]

Karlsson_otu_tab = otu_tab
y_sim = otu_tab@.Data
Karlsson_raw = y_sim
Karlsson_D_mat = D
Karlsson_hierarchy.info = Karlsson_hierarchy.info[-zero_taxa , ]



###CRC data
##read in Feng data
Feng_cd <- curatedMetagenomicData("FengQ_2015.relative_abundance", dryrun = FALSE, counts = TRUE)
tree_sum_exp <- Feng_cd[["2021-03-31.FengQ_2015.relative_abundance"]]
assay_data <- assays(tree_sum_exp)
relative_abundance_matrix <- assay_data[["relative_abundance"]]
Feng_OTU <- relative_abundance_matrix

col_data <- colData(tree_sum_exp)
str(col_data)
all_list_data <- col_data@listData
Feng_meta <- data.frame(all_list_data)
phylo_list <- rowTree(tree_sum_exp)

edges <- phylo_list$edge
row_ranges <- rowRanges(tree_sum_exp)
element_metadata <- row_ranges@elementMetadata
list_data <- element_metadata@listData
species <- list_data$species
genus <- list_data$genus
family <- list_data$family
Feng_hierarchy.info <- data.frame(
  plant_id = 1:length(species),
  species = species,
  genus = genus,
  family = family
)


#
Feng_OTU <- t(Feng_OTU)


#get a reasonable subset of meta features
head(Feng_meta)
chosen_meta <- Feng_meta[,c("study_condition", "age_category", "gender", "BMI", "number_reads")]
chosen_meta$age_category[chosen_meta$age_category == "adult"] = 1
chosen_meta$age_category[chosen_meta$age_category == "senior"] = 0
chosen_meta$gender[chosen_meta$gender == "male"] = 1
chosen_meta$gender[chosen_meta$gender == "female"] = 0
chosen_meta$BMI <- (chosen_meta$BMI-min(chosen_meta$BMI))/sqrt(var(chosen_meta$BMI))
chosen_meta$number_reads <- (chosen_meta$number_reads-min(chosen_meta$number_reads))/sqrt(var(chosen_meta$number_reads))

OTU_tab <- Feng_OTU

#get the distance matrix D base on phylogenetic tree
m = dim(OTU_tab)[2]
n = dim(OTU_tab)[1]
k = 50
nd_mat <- matrix(rep(1, k*m), k, m)
l <- rep(1,k)

n_taxa = m
tree_height = k

#generate the node set
for(i in 1:n_taxa){
  print(i)
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


d1_mat <- matrix(0, nrow = m, ncol = m)
taxa_vec <- match(1:m, edges[,2])
#generate the distance matrix
for(i in 1:m){
  for(j in 1:m){
    int_sc <- intersect(nd_mat[,i], nd_mat[,j])
    leni <- sum(!is.na(int_sc))
    len1 <- sum(!is.na(nd_mat[,i]))
    len2 <- sum(!is.na(nd_mat[,j]))
    d1_mat[i, j] = len1 - leni + 1 + len2 - leni + 1
  }
}
diag(d1_mat) = 0
#process the otu_table
otu_tab_processed <- OTU_tab / rowSums(OTU_tab) * 10^6
otu_tab_processed <- log10(otu_tab_processed+1.01)

#construct x
head(chosen_meta)
x <- cbind(rep(1, dim(chosen_meta)[1]), chosen_meta)

Feng_otable <- otu_tab_processed
Feng_meta_data <- x
Feng_D_mat <- d1_mat

Feng_condition = Feng_meta_data$study_condition
Feng_otable <- Feng_otable[Feng_condition != "adenoma",]
Feng_meta_data <- Feng_meta_data[Feng_condition != "adenoma",]
Feng_condition = Feng_condition [Feng_condition != "adenoma"]

# filter out the taxa that have too few reads in one of the two conditions.
zero_taxa <- c()
for(j in 1:dim(Feng_otable)[2]){
  for(cond in unique(Feng_condition)){
    cond_samp <- Feng_otable[Feng_condition == cond,j]
    if( sum(cond_samp != log10(1.01)) / length(cond_samp) < 0.1 ){
      zero_taxa <- c(zero_taxa, j)
    }
  }
}
zero_taxa <- unique(zero_taxa)

otu_tab = Feng_otable[,-zero_taxa]
D = Feng_D_mat[-zero_taxa , -zero_taxa ]
y_sim = otu_tab@.Data
Feng_raw = y_sim
Feng_otu_tab = otu_tab
Feng_D_mat = D
Feng_hierarchy.info = Feng_hierarchy.info[-zero_taxa , ]

##read in Vogtmann data
Vogtmann_cd <- curatedMetagenomicData("VogtmannE_2016.relative_abundance", dryrun = FALSE, counts = TRUE)
tree_sum_exp <- Vogtmann_cd[["2021-03-31.VogtmannE_2016.relative_abundance"]]
assay_data <- assays(tree_sum_exp)
relative_abundance_matrix <- assay_data[["relative_abundance"]]
Vogtmann_OTU <- relative_abundance_matrix

col_data <- colData(tree_sum_exp)
str(col_data)
all_list_data <- col_data@listData
Vogtmann_meta <- data.frame(all_list_data)

phylo_list <- rowTree(tree_sum_exp)
edges <- phylo_list$edge
row_ranges <- rowRanges(tree_sum_exp)
element_metadata <- row_ranges@elementMetadata
list_data <- element_metadata@listData
species <- list_data$species
genus <- list_data$genus
family <- list_data$family
Vogtmann_hierarchy.info <- data.frame(
  plant_id = 1:length(species),
  species = species,
  genus = genus,
  family = family
)


#
Vogtmann_OTU <- t(Vogtmann_OTU)


# leave out the missing conditons
na_idx <- which(is.na(Vogtmann_meta$study_condition))

Vogtmann_meta <- Vogtmann_meta[-na_idx,]


#get a reasonable subset of meta features
head(Vogtmann_meta)
chosen_meta <- Vogtmann_meta[,c("study_condition", "age_category", "gender", "BMI", "number_reads")]
# chosen_meta$study_condition[chosen_meta$study_condition == "adenoma"]=1
# chosen_meta$study_condition[chosen_meta$study_condition == "control"]=0
# chosen_meta$study_condition[chosen_meta$study_condition == "CRC"]=2
chosen_meta$age_category[chosen_meta$age_category == "adult"] = 1
chosen_meta$age_category[chosen_meta$age_category == "senior"] = 0
chosen_meta$gender[chosen_meta$gender == "male"] = 1
chosen_meta$gender[chosen_meta$gender == "female"] = 0
chosen_meta$BMI[is.na(chosen_meta$BMI)] <- mean(chosen_meta$BMI[!is.na(chosen_meta$BMI)])
chosen_meta$BMI <- (chosen_meta$BMI-min(chosen_meta$BMI))/sqrt(var(chosen_meta$BMI))
chosen_meta$number_reads <- (chosen_meta$number_reads-min(chosen_meta$number_reads))/sqrt(var(chosen_meta$number_reads))


OTU_tab <- Vogtmann_OTU
OTU_tab <- OTU_tab[-na_idx,]

#get the distance matrix D base on phylogenetic tree
m = dim(OTU_tab)[2]
n = dim(OTU_tab)[1]
k = 50
nd_mat <- matrix(rep(1, k*m), k, m)
l <- rep(1,k)

n_taxa = m
tree_height = k

#generate the node set
for(i in 1:n_taxa){
  print(i)
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


d1_mat <- matrix(0, nrow = m, ncol = m)
taxa_vec <- match(1:m, edges[,2])
#generate the distance matrix
for(i in 1:m){
  for(j in 1:m){
    int_sc <- intersect(nd_mat[,i], nd_mat[,j])
    leni <- sum(!is.na(int_sc))
    len1 <- sum(!is.na(nd_mat[,i]))
    len2 <- sum(!is.na(nd_mat[,j]))
    d1_mat[i, j] = len1 - leni + 1 + len2 - leni + 1
  }
}
diag(d1_mat) = 0

#process the otu_table
otu_tab_processed <- OTU_tab / rowSums(OTU_tab) * 10^6
otu_tab_processed <- log10(otu_tab_processed+1.01)

#construct x
head(chosen_meta)
x <- cbind(rep(1, dim(chosen_meta)[1]), chosen_meta)

Vogtmann_otable <- otu_tab_processed
Vogtmann_meta_data <- x
Vogtmann_D_mat <- d1_mat

Vogtmann_condition = Vogtmann_meta_data$study_condition

zero_taxa <- c()
for(j in 1:dim(Vogtmann_otable)[2]){
  for(cond in unique(Vogtmann_condition)){
    cond_samp <- Vogtmann_otable[Vogtmann_condition == cond,j]
    if( sum(cond_samp != log10(1.01)) / length(cond_samp) < 0.1 ){
      zero_taxa <- c(zero_taxa, j)
    }
  }
}
zero_taxa <- unique(zero_taxa)

otu_tab = Vogtmann_otable[,-zero_taxa]
D = Vogtmann_D_mat[-zero_taxa , -zero_taxa ]


y_sim = otu_tab@.Data
Vogtmann_raw = y_sim
Vogtmann_otu_tab = otu_tab
Vogtmann_D_mat = D
Vogtmann_hierarchy.info = Vogtmann_hierarchy.info[-zero_taxa , ]


##read in Yu data
Yu_cd <- curatedMetagenomicData("YuJ_2015.relative_abundance", dryrun = FALSE, counts = TRUE)
tree_sum_exp <- Yu_cd[["2021-03-31.YuJ_2015.relative_abundance"]]

assay_data <- assays(tree_sum_exp)
relative_abundance_matrix <- assay_data[["relative_abundance"]]
Yu_OTU <- relative_abundance_matrix
col_data <- colData(tree_sum_exp)
str(col_data)
all_list_data <- col_data@listData
Yu_meta <- data.frame(all_list_data)
phylo_list <- rowTree(tree_sum_exp)
edges <- phylo_list$edge
row_ranges <- rowRanges(tree_sum_exp)
element_metadata <- row_ranges@elementMetadata
list_data <- element_metadata@listData
species <- list_data$species
genus <- list_data$genus
family <- list_data$family
Yu_hierarchy.info <- data.frame(
  plant_id = 1:length(species),
  species = species,
  genus = genus,
  family = family
)



#
Yu_OTU <- t(Yu_OTU)


#get a reasonable subset of meta features
head(Yu_meta)
chosen_meta <- Yu_meta[,c("study_condition", "number_reads")]
# chosen_meta$study_condition[chosen_meta$study_condition == "adenoma"]=1
# chosen_meta$study_condition[chosen_meta$study_condition == "control"]=0
# chosen_meta$study_condition[chosen_meta$study_condition == "CRC"]=2
chosen_meta$number_reads <- (chosen_meta$number_reads-min(chosen_meta$number_reads))/sqrt(var(chosen_meta$number_reads))

OTU_tab <- Yu_OTU

#get the distance matrix D base on phylogenetic tree
m = dim(OTU_tab)[2]
n = dim(OTU_tab)[1]
k = 50
nd_mat <- matrix(rep(1, k*m), k, m)
l <- rep(1,k)

n_taxa = m
tree_height = k

#generate the node set
for(i in 1:n_taxa){
  print(i)
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


d1_mat <- matrix(0, nrow = m, ncol = m)
taxa_vec <- match(1:m, edges[,2])
#generate the distance matrix
for(i in 1:m){
  for(j in 1:m){
    int_sc <- intersect(nd_mat[,i], nd_mat[,j])
    leni <- sum(!is.na(int_sc))
    len1 <- sum(!is.na(nd_mat[,i]))
    len2 <- sum(!is.na(nd_mat[,j]))
    d1_mat[i, j] = len1 - leni + 1 + len2 - leni + 1
  }
}
diag(d1_mat) = 0

#process the otu_table
otu_tab_processed <- OTU_tab / rowSums(OTU_tab) * 10^6
otu_tab_processed <- log10(otu_tab_processed+1.01)

#construct x
head(chosen_meta)
x <- cbind(rep(1, dim(chosen_meta)[1]), chosen_meta)

Yu_otable <- otu_tab_processed
Yu_meta_data <- x
Yu_D_mat <- d1_mat

Yu_condition = Yu_meta_data$study_condition

zero_taxa <- c()
for(j in 1:dim(Yu_otable)[2]){
  for(cond in unique(Yu_condition)){
    cond_samp <- Yu_otable[Yu_condition == cond,j]
    if( sum(cond_samp != log10(1.01)) / length(cond_samp) < 0.1 ){
      zero_taxa <- c(zero_taxa, j)
    }
  }
}
zero_taxa <- unique(zero_taxa)

otu_tab = Yu_otable[,-zero_taxa]
D = Yu_D_mat[-zero_taxa , -zero_taxa ]


y_sim = otu_tab@.Data
Yu_raw = y_sim
Yu_otu_tab = otu_tab
Yu_D_mat = D
Yu_hierarchy.info = Yu_hierarchy.info[-zero_taxa , ]


##read in Zeller data
Zeller_cd <- curatedMetagenomicData("ZellerG_2014.relative_abundance", dryrun = FALSE, counts = TRUE)
tree_sum_exp <- Zeller_cd[["2021-03-31.ZellerG_2014.relative_abundance"]]
assay_data <- assays(tree_sum_exp)
relative_abundance_matrix <- assay_data[["relative_abundance"]]
Zeller_OTU <- relative_abundance_matrix

col_data <- colData(tree_sum_exp)
str(col_data)
all_list_data <- col_data@listData
Zeller_meta <- data.frame(all_list_data)

phylo_list <- rowTree(tree_sum_exp)
edges <- phylo_list$edge
row_ranges <- rowRanges(tree_sum_exp)
element_metadata <- row_ranges@elementMetadata
list_data <- element_metadata@listData
species <- list_data$species
genus <- list_data$genus
family <- list_data$family
Zeller_hierarchy.info <- data.frame(
  plant_id = 1:length(species),
  species = species,
  genus = genus,
  family = family
)


#
Zeller_OTU <- t(Zeller_OTU)


#get a reasonable subset of meta features
head(Zeller_meta)
chosen_meta <- Zeller_meta[,c("study_condition", "age_category", "gender", "BMI", "number_reads", "country")]
# chosen_meta$study_condition[chosen_meta$study_condition == "adenoma"]=1
# chosen_meta$study_condition[chosen_meta$study_condition == "control"]=0
# chosen_meta$study_condition[chosen_meta$study_condition == "CRC"]=2
chosen_meta$age_category[chosen_meta$age_category == "adult"] = 1
chosen_meta$age_category[chosen_meta$age_category == "senior"] = 0
chosen_meta$gender[chosen_meta$gender == "male"] = 1
chosen_meta$gender[chosen_meta$gender == "female"] = 0
chosen_meta$BMI[is.na(chosen_meta$BMI)] = mean(chosen_meta$BMI[!is.na(chosen_meta$BMI)])
chosen_meta$BMI <- (chosen_meta$BMI-min(chosen_meta$BMI))/sqrt(var(chosen_meta$BMI))
chosen_meta$number_reads <- (chosen_meta$number_reads-min(chosen_meta$number_reads))/sqrt(var(chosen_meta$number_reads))
chosen_meta$country[chosen_meta$country == "DEU"] = 1
chosen_meta$country[chosen_meta$country == "FRA"] = 0


OTU_tab <- Zeller_OTU

#get the distance matrix D base on phylogenetic tree
m = dim(OTU_tab)[2]
n = dim(OTU_tab)[1]
k = 50
nd_mat <- matrix(rep(1, k*m), k, m)
l <- rep(1,k)

n_taxa = m
tree_height = k

#generate the node set
for(i in 1:n_taxa){
  print(i)
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


d1_mat <- matrix(0, nrow = m, ncol = m)
#records the position of 1:m in the edges set.
taxa_vec <- match(1:m, edges[,2])
#generate the distance matrix
for(i in 1:m){
  for(j in 1:m){
    int_sc <- intersect(nd_mat[,i], nd_mat[,j])
    leni <- sum(!is.na(int_sc))
    len1 <- sum(!is.na(nd_mat[,i]))
    len2 <- sum(!is.na(nd_mat[,j]))
    d1_mat[i, j] = len1 - leni + 1 + len2 - leni + 1
  }
}
diag(d1_mat) = 0

#process the otu_table
otu_tab_processed <- OTU_tab / rowSums(OTU_tab) * 10^6
otu_tab_processed <- log10(otu_tab_processed+1.01)

#construct x
head(chosen_meta)
x <- cbind(rep(1, dim(chosen_meta)[1]), chosen_meta)

Zeller_otable <- otu_tab_processed
Zeller_meta_data <- x
Zeller_D_mat <- d1_mat

Zeller_condition = Zeller_meta_data$study_condition
Zeller_otable <- Zeller_otable[Zeller_condition != "adenoma",]
Zeller_meta_data <- Zeller_meta_data[Zeller_condition != "adenoma",]
Zeller_condition = Zeller_condition [Zeller_condition != "adenoma"]

zero_taxa <- c()
for(j in 1:dim(Zeller_otable)[2]){
  for(cond in unique(Zeller_condition)){
    cond_samp <- Zeller_otable[Zeller_condition == cond,j]
    if( sum(cond_samp != log10(1.01)) / length(cond_samp) < 0.1 ){
      zero_taxa <- c(zero_taxa, j)
    }
  }
}
zero_taxa <- unique(zero_taxa)

otu_tab = Zeller_otable[,-zero_taxa]
D = Zeller_D_mat[-zero_taxa , -zero_taxa ]


y_sim = otu_tab@.Data
Zeller_raw = y_sim
Zeller_otu_tab = otu_tab
Zeller_D_mat = D
Zeller_hierarchy.info = Zeller_hierarchy.info[-zero_taxa , ]




###Imputation:TphPMF
library(BHPMF)
#Qin_imp
dist_obj <- as.dist(Qin_D_mat)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 13)
genus_clusters <- cutree(hc, h = 30)
family_clusters <- cutree(hc, h = 45)
microbe_data <- data.frame(
  plant_id = 1:179,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)

microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
Qin_hierarchy.info<-microbe_data


tolerance <- .Machine$double.eps^0.5
Qin_raw[abs(Qin_raw - log10(1.01)) < tolerance] <- NA
y_sim <- t(Qin_raw)
trait.info<- y_sim
hierarchy.info <- Qin_hierarchy.info
setwd("/Users/hanxinyu/Desktop/real_data_analysis/Qin_BHPMF/")
tmp.dir<-dirname("/Users/hanxinyu/Desktop/real_data_analysis/Qin_BHPMF/tmp/")
GapFilling(trait.info,hierarchy.info,num.samples=26, burn=3, gaps=1,
           num.latent.feats=6,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)


data <- read.table("/Users/hanxinyu/Desktop/real_data_analysis/Qin_BHPMF/mean_gap_filled.txt", header = TRUE, sep = "\t")
bhmpf_pre <- as.matrix(data)
Qin_imp<-t(bhmpf_pre)
colnames(Qin_imp) <- colnames(Qin_raw)


#Karlsson_imp
dist_obj <- as.dist(Karlsson_D_mat)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 6)
genus_clusters <- cutree(hc, h = 30)
family_clusters <- cutree(hc, h = 45)
microbe_data <- data.frame(
  plant_id = 1:181,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)

microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
Karlsson_hierarchy.info<-microbe_data


tolerance <- .Machine$double.eps^0.5
Karlsson_raw[abs(Karlsson_raw - log10(1.01)) < tolerance] <- NA

y_sim <- t(Karlsson_raw)
trait.info<- y_sim
hierarchy.info <- Karlsson_hierarchy.info
setwd("/Users/hanxinyu/Desktop/real_data_analysis/Karlsson_BHPMF/")
tmp.dir<-dirname("/Users/hanxinyu/Desktop/real_data_analysis/Karlsson_BHPMF/tmp/")
GapFilling(trait.info,hierarchy.info,num.samples=26, burn=3, gaps=1,
           num.latent.feats=6,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)


data <- read.table("/Users/hanxinyu/Desktop/real_data_analysis/Karlsson_BHPMF/mean_gap_filled.txt", header = TRUE, sep = "\t")
bhmpf_pre <- as.matrix(data)
Karlsson_imp<-t(bhmpf_pre)
colnames(Karlsson_imp) <- colnames(Karlsson_raw)

#Feng_imp
dist_obj <- as.dist(Feng_D_mat)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 8)
genus_clusters <- cutree(hc, h = 30)
family_clusters <- cutree(hc, h = 46)
microbe_data <- data.frame(
  plant_id = 1:216,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)
microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
Feng_hierarchy.info<-microbe_data


tolerance <- .Machine$double.eps^0.5
Feng_raw[abs(Feng_raw - log10(1.01)) < tolerance] <- NA

y_sim <- t(Feng_raw)
trait.info<- y_sim
hierarchy.info <- Feng_hierarchy.info
setwd("/Users/hanxinyu/Desktop/real_data_analysis/Feng_BHPMF/")
tmp.dir<-dirname("/Users/hanxinyu/Desktop/real_data_analysis/Feng_BHPMF/tmp/")
GapFilling(trait.info,hierarchy.info,num.samples=26, burn=3, gaps=1,
           num.latent.feats=6,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)


data <- read.table("/Users/hanxinyu/Desktop/real_data_analysis/Feng_BHPMF/mean_gap_filled.txt", header = TRUE, sep = "\t")
bhmpf_pre <- as.matrix(data)
Feng_imp<-t(bhmpf_pre)
colnames(Feng_imp) <- colnames(Feng_raw)


#Vogtmann_imp
dist_obj <- as.dist(Vogtmann_D_mat)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 10)
genus_clusters <- cutree(hc, h = 21)
family_clusters <- cutree(hc, h = 46)
microbe_data <- data.frame(
  plant_id = 1:210,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)

microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
Vogtmann_hierarchy.info<-microbe_data


tolerance <- .Machine$double.eps^0.5
Vogtmann_raw[abs(Vogtmann_raw - log10(1.01)) < tolerance] <- NA

y_sim <- t(Vogtmann_raw)
trait.info<- y_sim
hierarchy.info <- Vogtmann_hierarchy.info
setwd("/Users/hanxinyu/Desktop/real_data_analysis/Vogtmann_BHPMF/")
tmp.dir<-dirname("/Users/hanxinyu/Desktop/real_data_analysis/Vogtmann_BHPMF/tmp/")
GapFilling(trait.info,hierarchy.info,num.samples=26, burn=3, gaps=1,
           num.latent.feats=6,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)


data <- read.table("/Users/hanxinyu/Desktop/real_data_analysis/Vogtmann_BHPMF/mean_gap_filled.txt", header = TRUE, sep = "\t")
bhmpf_pre <- as.matrix(data)
Vogtmann_imp<-t(bhmpf_pre)
colnames(Vogtmann_imp) <- colnames(Vogtmann_raw)


#Yu_imp
dist_obj <- as.dist(Yu_D_mat)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 13)
genus_clusters <- cutree(hc, h = 30)
family_clusters <- cutree(hc, h = 45)
microbe_data <- data.frame(
  plant_id = 1:199,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)

microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
Yu_hierarchy.info<-microbe_data


tolerance <- .Machine$double.eps^0.5
Yu_raw[abs(Yu_raw - log10(1.01)) < tolerance] <- NA

y_sim <- t(Yu_raw)
trait.info<- y_sim
hierarchy.info <- Yu_hierarchy.info
setwd("/Users/hanxinyu/Desktop/real_data_analysis/Yu_BHPMF/")
tmp.dir<-dirname("/Users/hanxinyu/Desktop/real_data_analysis/Yu_BHPMF/tmp/")
GapFilling(trait.info,hierarchy.info,num.samples=26, burn=3, gaps=1,
           num.latent.feats=6,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)


data <- read.table("/Users/hanxinyu/Desktop/real_data_analysis/Yu_BHPMF/mean_gap_filled.txt", header = TRUE, sep = "\t")
bhmpf_pre <- as.matrix(data)
Yu_imp<-t(bhmpf_pre)
colnames(Yu_imp) <- colnames(Yu_raw)


#Zeller_imp
dist_obj <- as.dist(Zeller_D_mat)
hc <- hclust(dist_obj, method = "complete")
species_clusters <- cutree(hc, h = 9)
genus_clusters <- cutree(hc, h = 30)
family_clusters <- cutree(hc, h = 46)
microbe_data <- data.frame(
  plant_id = 1:237,
  species = as.factor(species_clusters),
  genus = as.factor(genus_clusters),
  family = as.factor(family_clusters)
)
microbe_data$species <- paste("S", microbe_data$species, sep = "_")
microbe_data$genus <- paste("G", microbe_data$genus, sep = "_")
microbe_data$family <- paste("F", microbe_data$family, sep = "_")
Zeller_hierarchy.info<-microbe_data


tolerance <- .Machine$double.eps^0.5
Zeller_raw[abs(Zeller_raw - log10(1.01)) < tolerance] <- NA
y_sim <- t(Zeller_raw)
trait.info<- y_sim
hierarchy.info <- Zeller_hierarchy.info
setwd("/Users/hanxinyu/Desktop/real_data_analysis/Zeller_BHPMF/")
tmp.dir<-dirname("/Users/hanxinyu/Desktop/real_data_analysis/Zeller_BHPMF/tmp/")
GapFilling(trait.info,hierarchy.info,num.samples=26, burn=3, gaps=1,
           num.latent.feats=6,num.folds.tuning=10,mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled.txt"),std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled.txt"),tmp.dir=tmp.dir,rmse.plot.test.data=TRUE)


data <- read.table("/Users/hanxinyu/Desktop/real_data_analysis/Zeller_BHPMF/mean_gap_filled.txt", header = TRUE, sep = "\t")
bhmpf_pre <- as.matrix(data)
Zeller_imp<-t(bhmpf_pre)
colnames(Zeller_imp) <- colnames(Zeller_raw)



##Required R Packages for classification
if (!require("glmnet")) {
  install.packages("glmnet", dependencies = TRUE)
  library(glmnet)
}
if (!require("randomForest")) {
  install.packages("randomForest", dependencies = TRUE)
  library(randomForest)
}
if (!require("e1071")) {
  install.packages("e1071", dependencies = TRUE)
  library(e1071)
}
if (!require("neuralnet")) {
  install.packages("neuralnet", dependencies = TRUE)
  library(neuralnet)
}
if (!require("xgboost")) {
  install.packages("xgboost", dependencies = TRUE)
  library(xgboost)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("lattice")) {
  install.packages("lattice", dependencies = TRUE)
  library(lattice)
}
if (!require("factoextra")) {
  install.packages("factoextra", dependencies = TRUE)
  library(factoextra)
}
if (!require("PRROC")) {
  install.packages("PRROC", dependencies = TRUE)
  library(PRROC)
}
if (!require("nproc")) {
  install.packages("nproc", dependencies = TRUE)
  library(nproc)
}
if (!require("igraph")) {
  install.packages("igraph", dependencies = TRUE)
  library(igraph)
}
if (!require("bnviewer")) {
  install.packages("bnviewer", dependencies = TRUE)
  library(bnviewer)
}
library(graph)
require("parallel")
library(bnlearn)




##Functions
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
  otu_colnames <- setdiff(colnames(otu_data), "Sample.ID")
  var_colnames <- setdiff(colnames(var_data), "Sample.ID")
  data_comp <- merge(otu_data, var_data, by="Sample.ID", all.y=TRUE)
  data_comp_colnames <- c("Sample.ID", otu_colnames, var_colnames)
  data_comp <- setNames(data_comp, data_comp_colnames)
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
  
  #########################################
  ### Code to extract surrogate p-value
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
  
  
  ### Bubble plot
  
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

real_data_meta_analysis <- function(OTU_mat, condition, meta_data, p_thresh = 0.05){
  OTU_count <- floor(10^OTU_mat + 1.01)
  taxa_species <- colnames(OTU_mat) 
  # Wilcoxon test
  less_pval_wilcoxon <- apply(OTU_count, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
  greater_pval_wilcoxon <- apply(OTU_count, 2, FUN = function(x) pairwise.wilcox.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
  
  less_identified <- which(p.adjust(less_pval_wilcoxon, method = "fdr") < p_thresh)
  greater_identified <- which(p.adjust(greater_pval_wilcoxon, method = "fdr") < p_thresh)
  
  Wilcox_result <- list("T2D_leq_cond" = taxa_species[less_identified], "T2D_grt_cond" = taxa_species[greater_identified])
  
  # t test
  less_pval_t_test <- apply(OTU_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "less", p.adjust.method = "none")$p.value)
  greater_pval_t_test <- apply(OTU_mat, 2, FUN = function(x) pairwise.t.test(x, condition, alternative = "greater", p.adjust.method = "none")$p.value)
  
  less_identified <- which(p.adjust(less_pval_t_test, method = "fdr") < p_thresh)
  greater_identified <- which(p.adjust(greater_pval_t_test, method = "fdr") < p_thresh)
  t_test_result <- list("T2D_leq_cond" = taxa_species[less_identified], "T2D_grt_cond" = taxa_species[greater_identified])
  
  # ANCOM
  Vardat <- cbind(meta_data, condition)
  if(is.null(dim(meta_data))){
    colnames(Vardat) <- c("feature1", "condition")
  }else{
    colnames(Vardat) <- c(colnames(meta_data), "condition")
  }
  Vardat <- cbind(1:length(condition), Vardat)
  colnames(Vardat)[1] <- "Sample.ID"
  OTUdat <- OTU_count
  OTUdat <- cbind(1:length(condition), OTUdat)
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
                             sig=p_thresh,
                             prev.cut=0.90,
                             longitudinal = F)
  
  ANCOM_detected <- comparison_test$W.taxa$otu.names[comparison_test$W.taxa$detected_0.6]
  ANCOM_result <- list("T2D_leq_cond" = ANCOM_detected, "T2D_grt_cond" = ANCOM_detected)
  
  # ZINB
  Vardat <- cbind(meta_data, condition)
  if(is.null(dim(meta_data))){
    colnames(Vardat) <- c("feature1", "condition")
  }else{
    colnames(Vardat) <- c(colnames(meta_data), "condition")
  }
  
  OTUdat <- OTU_count
  
  Vardat <- as.data.frame(Vardat)
  OTUdat <- as.data.frame(OTUdat)
  
  library(pscl)
  
  df_taxa <- c()
  p_value_zinb <- c()
  for(i in 1:dim(OTUdat)[2]){
    print(i)
   
    out <- tryCatch({
      sr <- summary(ml <- zeroinfl(OTUdat[,i] ~ condition, dist = "negbin"))
    }, error = function(cond) {
      message("error")
      return(TRUE)
    },
    warning = function(cond) {
      message("warning")
      return(TRUE)
    },
    finally =  {
      p_val = NA
    }
    )
    if(!(isTRUE(out))){
      if(sr$coefficients$count[2,1] > 0){
        p_val <- sr$coefficients$count[2,4]
      }else{
        p_val <- 1
      }
      
    }
    p_value_zinb <- c(p_value_zinb, p_val)
  }
  p_value_zinb[is.na(p_value_zinb)] <- 1
  p_value_zinb_greater <- p_value_zinb
  ZINB_ident_greater <- which(p.adjust(p_value_zinb, method = "fdr") < p_thresh)
  
  df_taxa <- c()
  p_value_zinb <- c()
  for(i in 1:dim(OTUdat)[2]){
    print(i)
    
    out <- tryCatch({
      sr <- summary(ml <- zeroinfl(OTUdat[,i] ~ condition, dist = "negbin"))
    }, error = function(cond) {
      message("error")
      return(TRUE)
    },
    warning = function(cond) {
      message("warning")
      return(TRUE)
    },
    finally =  {
      p_val = NA
    }
    )
    if(!(isTRUE(out))){
      if(sr$coefficients$count[2,1] < 0){
        p_val <- sr$coefficients$count[2,4]
      }else{
        p_val <- 1
      }
      
    }
    p_value_zinb <- c(p_value_zinb, p_val)
  }
  p_value_zinb_less <- p_value_zinb
  p_value_zinb_less[is.na(p_value_zinb_less)] <- 1
  ZINB_ident_less <- which(p.adjust(p_value_zinb, method = "fdr") < p_thresh)
  ZINB_result <- list("T2D_leq_cond" = taxa_species[ZINB_ident_less], "T2D_grt_cond" = taxa_species[ZINB_ident_greater])
  
  # NB
  library(MASS)
  Vardat <- cbind(meta_data, condition)
  if(is.null(dim(meta_data))){
    colnames(Vardat) <- c("feature1", "condition")
  }else{
    colnames(Vardat) <- c(colnames(meta_data), "condition")
  }
  
  OTUdat <- OTU_count
  
  Vardat <- as.data.frame(Vardat)
  OTUdat <- as.data.frame(OTUdat)
  df_taxa <- c()
  p_value_nb <- c()
  for(i in 1:dim(OTUdat)[2]){
    print(i)
    out <- tryCatch({
      snb <- summary(glm.nb(OTUdat[,i] ~ condition))
    }, error = function(cond) {
      message("error")
      return(TRUE)
    },
    warning = function(cond) {
      message("warning")
      return(TRUE)
    },
    finally =  {
      p_val = NA
    })
  
    if(!(isTRUE(out))){
      if(snb[12]$coefficients[2,1] > 0){
        p_val <- snb[12]$coefficients[2,4]
      }else{
        p_val <- 1
      }
    }
    p_value_nb <- c(p_value_nb, p_val)
  }
  p_value_nb[is.na(p_value_nb)] <- 1
  p_value_nb_greater <- p_value_nb
  nb_ident_greater <- which(p.adjust(p_value_nb, method = "fdr") < p_thresh)
  
  df_taxa <- c()
  p_value_nb <- c()
  for(i in 1:dim(OTUdat)[2]){
    print(i)
    out <- tryCatch({
      snb <- summary(glm.nb(OTUdat[,i] ~ condition))
    }, error = function(cond) {
      message("error")
      return(TRUE)
    },
    warning = function(cond) {
      message("warning")
      return(FALSE)
    },
    finally =  {
      p_val = NA
    })
    
    
    if(!(isTRUE(out))){
      if(snb[12]$coefficients[2,1] < 0){
        p_val <- snb[12]$coefficients[2,4]
      }else{
        p_val <- 1
      }
    }
    p_value_nb <- c(p_value_nb, p_val)
  }
  p_value_nb_less <- p_value_nb
  p_value_nb_less[is.na(p_value_nb_less)] <- 1
  nb_ident_less <- which(p.adjust(p_value_nb, method = "fdr") < p_thresh)
  NB_result <- list("T2D_leq_cond" = taxa_species[nb_ident_less], "T2D_grt_cond" = taxa_species[nb_ident_greater])
  
  # metagenomic seq
  otumat = t(OTU_count)
  OTU = otu_table(otumat, taxa_are_rows = TRUE)
  physeq = phyloseq(OTU)
  rownames(Vardat) <- rownames(OTU_count)
  sampledata = sample_data(Vardat)
  physeq2 = phyloseq(OTU,sampledata)
  
  mseq_obj <- phyloseq_to_metagenomeSeq(physeq2)
  pd <- pData(mseq_obj)
  mod <- model.matrix(~1 + condition, data = pd)
  ran_seq = fitFeatureModel(mseq_obj, mod)
  less_identified <- which( ran_seq@fitZeroLogNormal$logFC < 0 & p.adjust(ran_seq@pvalues, method = "fdr") < p_thresh)
  greater_identified <- which( ran_seq@fitZeroLogNormal$logFC > 0 & p.adjust(ran_seq@pvalues, method = "fdr") < p_thresh)
  Metagenomicseq_result <- list("T2D_leq_cond" = taxa_species[less_identified], "T2D_grt_cond" = taxa_species[greater_identified])
  
  # DEseq2
  Deseq2_obj <- phyloseq_to_deseq2(physeq2, ~ condition)
  results = DESeq(Deseq2_obj, test="Wald", fitType="parametric")
  res = results(results, cooksCutoff = FALSE)
  less_identified <- which( res$log2FoldChange < 0 & p.adjust(res$pvalue, method = "fdr") < p_thresh)
  greater_identified <- which( res$log2FoldChange > 0 & p.adjust(res$pvalue, method = "fdr") < p_thresh)
  DEseq2_result <- list("T2D_leq_cond" = taxa_species[less_identified], "T2D_grt_cond" = taxa_species[greater_identified])
  
  # Omnibus test
  #library(mbzinb)
  Vardat_zinb <- Vardat
  OTUdat_zinb <- as.matrix(t(OTUdat))
  colnames(OTUdat_zinb) <- rownames(Vardat_zinb)
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
  mbzinb_data <- mbzinb.dataset(OTUdat_zinb, Vardat_zinb)
  mbzinb_test_result <- mbzinb.test(mbzinb_data, group = "condition")
  less_identified <- which( res$log2FoldChange < 0 & p.adjust(mbzinb_test_result$results$PValue, method = "fdr") < p_thresh)
  greater_identified <- which( res$log2FoldChange > 0 & p.adjust(mbzinb_test_result$results$PValue, method = "fdr") < p_thresh)
  mbzinb_result <- list("T2D_leq_cond" = taxa_species[less_identified], "T2D_grt_cond" = taxa_species[greater_identified])
  
  return(result = list("Wilcox" = Wilcox_result, "t_test" = t_test_result, "ANCOM" = ANCOM_result, "ZINB" = ZINB_result, "NB" = NB_result, "metagenomicseq" = Metagenomicseq_result, "DEseq2" = DEseq2_result, "Omnibus" = mbzinb_result, "DEseq2_pval" = res$pvalue, "Omnibus_pval" = mbzinb_test_result$results$PValue))
}


library(exactRankTests)
library(pscl)
Qin_condition <- Qin_meta_data$study_condition
Qin_raw_result <- real_data_meta_analysis(OTU_mat = Qin_raw, condition = Qin_condition, meta_data = Qin_meta_data[,-c(1,2)], p_thresh = 0.05)
Qin_imp_result <- real_data_meta_analysis(OTU_mat = Qin_imp, condition = Qin_condition, meta_data = Qin_meta_data[,-c(1,2)], p_thresh = 0.05)


Karlsson_condition <- Karlsson_meta_data$study_condition
Karlsson_raw_result <- real_data_meta_analysis(Karlsson_raw, Karlsson_condition, meta_data = Karlsson_meta_data[,-c(1,2)], p_thresh = 0.05)
Karlsson_imp_result <- real_data_meta_analysis(Karlsson_imp, Karlsson_condition, meta_data = Karlsson_meta_data[,-c(1,2)], p_thresh = 0.05)

Feng_condition <- Feng_meta_data$study_condition
Feng_raw_result <- real_data_meta_analysis(Feng_raw, Feng_condition, meta_data = Feng_meta_data[,-c(1,2)], p_thresh = 0.05)
Feng_imp_result <- real_data_meta_analysis(Feng_imp, Feng_condition, meta_data = Feng_meta_data[,-c(1,2)], p_thresh = 0.05)

Vogtmann_condition <- Vogtmann_meta_data$study_condition
Vogtmann_raw_result <- real_data_meta_analysis(Vogtmann_raw, Vogtmann_condition, meta_data = Vogtmann_meta_data[,-c(1,2)], p_thresh = 0.05)
Vogtmann_imp_result <- real_data_meta_analysis(Vogtmann_imp, Vogtmann_condition, meta_data = Vogtmann_meta_data[,-c(1,2)], p_thresh = 0.05)

Yu_condition <- Yu_meta_data$study_condition
Yu_raw_result <- real_data_meta_analysis(Yu_raw, Yu_condition, meta_data = Yu_meta_data[,-c(1,2)], p_thresh = 0.05)
Yu_imp_result <- real_data_meta_analysis(Yu_imp, Yu_condition, meta_data = Yu_meta_data[,-c(1,2)], p_thresh = 0.05)

Zeller_condition <- Zeller_meta_data$study_condition
Zeller_raw_result <- real_data_meta_analysis(Zeller_raw, Zeller_condition, meta_data = Zeller_meta_data[,-c(1,2)], p_thresh = 0.05)
Zeller_imp_result <- real_data_meta_analysis(Zeller_imp, Zeller_condition, meta_data = Zeller_meta_data[,-c(1,2)], p_thresh = 0.05)



DA_result <- list()
for(i in 1:7){
  raw_leq <- length(intersect( Qin_raw_result[[i]]$T2D_leq_cond, Karlsson_raw_result[[i]]$T2D_leq_cond ))
  raw_grt <- length(intersect( Qin_raw_result[[i]]$T2D_grt_cond, Karlsson_raw_result[[i]]$T2D_grt_cond ))
  
  imp_leq <- length(intersect( Qin_imp_result[[i]]$T2D_leq_cond, Karlsson_imp_result[[i]]$T2D_leq_cond ))
  imp_grt <-  length(intersect( Qin_imp_result[[i]]$T2D_grt_cond, Karlsson_imp_result[[i]]$T2D_grt_cond ))
  DA_result[[i]] <- c("raw_leq" =  raw_leq, "raw_grt" = raw_grt, "imp_leq" = imp_leq, "imp_grt" = imp_grt)
}

Qin_Wilcox <- length(Qin_raw_result$Wilcox$T2D_leq_cond)+length(Qin_raw_result$Wilcox$T2D_grt_cond)
Karlsson_Wilcox <-length(Karlsson_raw_result$Wilcox$T2D_leq_cond)+length(Karlsson_raw_result$Wilcox$T2D_grt_cond)
Feng_Wilcox <-length(Feng_raw_result$Wilcox$T2D_leq_cond)+length(Feng_raw_result$Wilcox$T2D_grt_cond)
Vogtmann_Wilcox <-length(Vogtmann_raw_result$Wilcox$T2D_leq_cond)+length(Vogtmann_raw_result$Wilcox$T2D_grt_cond)
Yu_Wilcox <-length(Yu_raw_result$Wilcox$T2D_leq_cond)+length(Yu_raw_result$Wilcox$T2D_grt_cond)
Zeller_Wilcox <-length(Zeller_raw_result$Wilcox$T2D_leq_cond)+length(Zeller_raw_result$Wilcox$T2D_grt_cond)

Qin_ANCOM <-length(Qin_raw_result$ANCOM$T2D_leq_cond)+length(Qin_raw_result$ANCOM$T2D_grt_cond)
Karlsson_ANCOM <-length(Karlsson_raw_result$ANCOM$T2D_leq_cond)+length(Karlsson_raw_result$ANCOM$T2D_grt_cond)
Feng_ANCOM <-length(Feng_raw_result$ANCOM$T2D_leq_cond)+length(Feng_raw_result$ANCOM$T2D_grt_cond)
Vogtmann_ANCOM <-length(Vogtmann_raw_result$ANCOM$T2D_leq_cond)+length(Vogtmann_raw_result$ANCOM$T2D_grt_cond)
Yu_ANCOM <-length(Yu_raw_result$ANCOM$T2D_leq_cond)+length(Yu_raw_result$ANCOM$T2D_grt_cond)
Zeller_ANCOM <-length(Zeller_raw_result$ANCOM$T2D_leq_cond)+length(Zeller_raw_result$ANCOM$T2D_grt_cond)

Qin_metagenomicseq <-length(Qin_raw_result$metagenomicseq$T2D_leq_cond)+length(Qin_raw_result$metagenomicseq$T2D_grt_cond)
Karlsson_metagenomicseq <-length(Karlsson_raw_result$metagenomicseq$T2D_leq_cond)+length(Karlsson_raw_result$metagenomicseq$T2D_grt_cond)
Feng_metagenomicseq <-length(Feng_raw_result$metagenomicseq$T2D_leq_cond)+length(Feng_raw_result$metagenomicseq$T2D_grt_cond)
Vogtmann_metagenomicseq <-length(Vogtmann_raw_result$metagenomicseq$T2D_leq_cond)+length(Vogtmann_raw_result$metagenomicseq$T2D_grt_cond)
Yu_metagenomicseq <-length(Yu_raw_result$metagenomicseq$T2D_leq_cond)+length(Yu_raw_result$metagenomicseq$T2D_grt_cond)
Zeller_metagenomicseq <-length(Zeller_raw_result$metagenomicseq$T2D_leq_cond)+length(Zeller_raw_result$metagenomicseq$T2D_grt_cond)


Qin_DEseq2 <-length(Qin_raw_result$DEseq2$T2D_leq_cond)+length(Qin_raw_result$DEseq2$T2D_grt_cond)
Karlsson_DEseq2 <-length(Karlsson_raw_result$DEseq2$T2D_leq_cond)+length(Karlsson_raw_result$DEseq2$T2D_grt_cond)
Feng_DEseq2 <-length(Feng_raw_result$DEseq2$T2D_leq_cond)+length(Feng_raw_result$DEseq2$T2D_grt_cond)
Vogtmann_DEseq2 <-length(Vogtmann_raw_result$DEseq2$T2D_leq_cond)+length(Vogtmann_raw_result$DEseq2$T2D_grt_cond)
Yu_DEseq2 <-length(Yu_raw_result$DEseq2$T2D_leq_cond)+length(Yu_raw_result$DEseq2$T2D_grt_cond)
Zeller_DEseq2 <-length(Zeller_raw_result$DEseq2$T2D_leq_cond)+length(Zeller_raw_result$DEseq2$T2D_grt_cond)

Qin_Omnibus <-length(Qin_raw_result$Omnibus$T2D_leq_cond)+length(Qin_raw_result$Omnibus$T2D_grt_cond)
Karlsson_Omnibus <-length(Karlsson_raw_result$Omnibus$T2D_leq_cond)+length(Karlsson_raw_result$Omnibus$T2D_grt_cond)
Feng_Omnibus <-length(Feng_raw_result$Omnibus$T2D_leq_cond)+length(Feng_raw_result$Omnibus$T2D_grt_cond)
Vogtmann_Omnibus <-length(Vogtmann_raw_result$Omnibus$T2D_leq_cond)+length(Vogtmann_raw_result$Omnibus$T2D_grt_cond)
Yu_Omnibus <-length(Yu_raw_result$Omnibus$T2D_leq_cond)+length(Yu_raw_result$Omnibus$T2D_grt_cond)
Zeller_Omnibus <-length(Zeller_raw_result$Omnibus$T2D_leq_cond)+length(Zeller_raw_result$Omnibus$T2D_grt_cond)

results_table <- data.frame(
  Wilcoxon = c(Qin_Wilcox,Karlsson_Wilcox,Feng_Wilcox,Vogtmann_Wilcox,Yu_Wilcox,Zeller_Wilcox), 
  ANCOM = c(Qin_ANCOM,Karlsson_ANCOM,Feng_ANCOM,Vogtmann_ANCOM,Yu_ANCOM,Zeller_ANCOM),
  MetagenomeSeq = c(Qin_metagenomicseq,Karlsson_metagenomicseq,Feng_metagenomicseq,Vogtmann_metagenomicseq,Yu_metagenomicseq,Zeller_metagenomicseq),
  DESeq2_phyloseq = c(Qin_DEseq2,Karlsson_DEseq2,Feng_DEseq2,Vogtmann_DEseq2,Yu_DEseq2,Zeller_DEseq2),
  Omnibust_test = c(Qin_Omnibus,Karlsson_Omnibus,Feng_Omnibus,Vogtmann_Omnibus,Yu_Omnibus,Zeller_Omnibus)
  
)

rownames(results_table) <- c("Qin et al.", "Karlsson et al.", "Feng et al.", "Vogtmann et al.", "Yu et al.", "Zeller et al.")



Qin_Wilcox <- length(Qin_imp_result$Wilcox$T2D_leq_cond)+length(Qin_imp_result$Wilcox$T2D_grt_cond)
Karlsson_Wilcox <-length(Karlsson_imp_result$Wilcox$T2D_leq_cond)+length(Karlsson_imp_result$Wilcox$T2D_grt_cond)
Feng_Wilcox <-length(Feng_imp_result$Wilcox$T2D_leq_cond)+length(Feng_imp_result$Wilcox$T2D_grt_cond)
Vogtmann_Wilcox <-length(Vogtmann_imp_result$Wilcox$T2D_leq_cond)+length(Vogtmann_imp_result$Wilcox$T2D_grt_cond)
Yu_Wilcox <-length(Yu_imp_result$Wilcox$T2D_leq_cond)+length(Yu_imp_result$Wilcox$T2D_grt_cond)
Zeller_Wilcox <-length(Zeller_imp_result$Wilcox$T2D_leq_cond)+length(Zeller_imp_result$Wilcox$T2D_grt_cond)

Qin_ANCOM <-length(Qin_imp_result$ANCOM$T2D_leq_cond)+length(Qin_imp_result$ANCOM$T2D_grt_cond)
Karlsson_ANCOM <-length(Karlsson_imp_result$ANCOM$T2D_leq_cond)+length(Karlsson_imp_result$ANCOM$T2D_grt_cond)
Feng_ANCOM <-length(Feng_imp_result$ANCOM$T2D_leq_cond)+length(Feng_imp_result$ANCOM$T2D_grt_cond)
Vogtmann_ANCOM <-length(Vogtmann_imp_result$ANCOM$T2D_leq_cond)+length(Vogtmann_imp_result$ANCOM$T2D_grt_cond)
Yu_ANCOM <-length(Yu_imp_result$ANCOM$T2D_leq_cond)+length(Yu_imp_result$ANCOM$T2D_grt_cond)
Zeller_ANCOM <-length(Zeller_imp_result$ANCOM$T2D_leq_cond)+length(Zeller_imp_result$ANCOM$T2D_grt_cond)

Qin_metagenomicseq <-length(Qin_imp_result$metagenomicseq$T2D_leq_cond)+length(Qin_imp_result$metagenomicseq$T2D_grt_cond)
Karlsson_metagenomicseq <-length(Karlsson_imp_result$metagenomicseq$T2D_leq_cond)+length(Karlsson_imp_result$metagenomicseq$T2D_grt_cond)
Feng_metagenomicseq <-length(Feng_imp_result$metagenomicseq$T2D_leq_cond)+length(Feng_imp_result$metagenomicseq$T2D_grt_cond)
Vogtmann_metagenomicseq <-length(Vogtmann_imp_result$metagenomicseq$T2D_leq_cond)+length(Vogtmann_imp_result$metagenomicseq$T2D_grt_cond)
Yu_metagenomicseq <-length(Yu_imp_result$metagenomicseq$T2D_leq_cond)+length(Yu_imp_result$metagenomicseq$T2D_grt_cond)
Zeller_metagenomicseq <-length(Zeller_imp_result$metagenomicseq$T2D_leq_cond)+length(Zeller_imp_result$metagenomicseq$T2D_grt_cond)


Qin_DEseq2 <-length(Qin_imp_result$DEseq2$T2D_leq_cond)+length(Qin_imp_result$DEseq2$T2D_grt_cond)
Karlsson_DEseq2 <-length(Karlsson_imp_result$DEseq2$T2D_leq_cond)+length(Karlsson_imp_result$DEseq2$T2D_grt_cond)
Feng_DEseq2 <-length(Feng_imp_result$DEseq2$T2D_leq_cond)+length(Feng_imp_result$DEseq2$T2D_grt_cond)
Vogtmann_DEseq2 <-length(Vogtmann_imp_result$DEseq2$T2D_leq_cond)+length(Vogtmann_imp_result$DEseq2$T2D_grt_cond)
Yu_DEseq2 <-length(Yu_imp_result$DEseq2$T2D_leq_cond)+length(Yu_imp_result$DEseq2$T2D_grt_cond)
Zeller_DEseq2 <-length(Zeller_imp_result$DEseq2$T2D_leq_cond)+length(Zeller_imp_result$DEseq2$T2D_grt_cond)

Qin_Omnibus <-length(Qin_imp_result$Omnibus$T2D_leq_cond)+length(Qin_imp_result$Omnibus$T2D_grt_cond)
Karlsson_Omnibus <-length(Karlsson_imp_result$Omnibus$T2D_leq_cond)+length(Karlsson_imp_result$Omnibus$T2D_grt_cond)
Feng_Omnibus <-length(Feng_imp_result$Omnibus$T2D_leq_cond)+length(Feng_imp_result$Omnibus$T2D_grt_cond)
Vogtmann_Omnibus <-length(Vogtmann_imp_result$Omnibus$T2D_leq_cond)+length(Vogtmann_imp_result$Omnibus$T2D_grt_cond)
Yu_Omnibus <-length(Yu_imp_result$Omnibus$T2D_leq_cond)+length(Yu_imp_result$Omnibus$T2D_grt_cond)
Zeller_Omnibus <-length(Zeller_imp_result$Omnibus$T2D_leq_cond)+length(Zeller_imp_result$Omnibus$T2D_grt_cond)

results_table_imp <- data.frame(
  Wilcoxon = c(Qin_Wilcox,Karlsson_Wilcox,Feng_Wilcox,Vogtmann_Wilcox,Yu_Wilcox,Zeller_Wilcox), 
  ANCOM = c(Qin_ANCOM,Karlsson_ANCOM,Feng_ANCOM,Vogtmann_ANCOM,Yu_ANCOM,Zeller_ANCOM),
  MetagenomeSeq = c(Qin_metagenomicseq,Karlsson_metagenomicseq,Feng_metagenomicseq,Vogtmann_metagenomicseq,Yu_metagenomicseq,Zeller_metagenomicseq),
  DESeq2_phyloseq = c(Qin_DEseq2,Karlsson_DEseq2,Feng_DEseq2,Vogtmann_DEseq2,Yu_DEseq2,Zeller_DEseq2),
  Omnibust_test = c(Qin_Omnibus,Karlsson_Omnibus,Feng_Omnibus,Vogtmann_Omnibus,Yu_Omnibus,Zeller_Omnibus)
  
)

rownames(results_table_imp) <- c("Qin et al.", "Karlsson et al.", "Feng et al.", "Vogtmann et al.", "Yu et al.", "Zeller et al.")
print(results_table_imp)


set.seed(2020)

Zeller_data <- Zeller_raw
Zeller_data <- apply(Zeller_data, 2, FUN = function(x){
  return(as.numeric(as.character(x)))
})

input_data <- as.matrix(Zeller_data) 

Zeller_condition[Zeller_condition == "control"]=0
Zeller_condition[Zeller_condition == "CRC"]=1
Zeller_condition <- as.numeric(Zeller_condition)
condition <- Zeller_condition
features_raw <- Zeller_raw_result
features_imputed <- Zeller_imp_result

classification_results <- function(input_data, condition){
  # Use 10xCV to generate CV_testing error.
  # method used are logistic regression, SVM, elastic net, randomforest and XGboost
  
  # single result for all the taxa species.
  ## Define the K folds
  K <- 5
  CV_split_idx <- sample(1:K, size=length(input_data[,1]), replace=T)
  
  input_data <- data.frame(input_data)
  cv_lr <- sapply(1:K, FUN=function(k) {
    test_idx <- which(CV_split_idx == k)
    x_train <- input_data[-test_idx,]
    x_test <- input_data[test_idx,]
    y_train <- condition[-test_idx]
    y_test <- condition[test_idx]
    data_train <- cbind(y_train, x_train)
    colnames(data_train)[1] <- "y"
    
    lr_model.orig <- glm(y~., family=binomial(), data=data_train)
    lr_pred.orig <- predict(lr_model.orig, newdata=x_test, type="response")
    lr_prauc.orig <- pr.curve(scores.class0=lr_pred.orig, weights.class0=as.numeric(as.factor(y_test)) -1)$auc.integral
    lr_rocauc.orig <- roc.curve(scores.class0=lr_pred.orig, weights.class0=as.numeric(as.factor(y_test)) -1)$auc
    
    return(list("prauc" = lr_prauc.orig, "rocauc" = lr_rocauc.orig))
  })
  
  ## logistic regression with elastic net penalization
  ### choose the lambda
  input_data <- as.matrix(input_data)
  lambda.enet <- cv.glmnet(x=input_data, y= as.factor(condition), alpha=0.5, family="binomial")$lambda.min
  cv_lr.enet <- sapply(1:K, FUN=function(k) {
    test_idx <- which(CV_split_idx == k)
    x_train <- input_data[-test_idx,]
    x_test <- input_data[test_idx,]
    y_train <- condition[-test_idx]
    y_test <- condition[test_idx]
    ## original sample ratio
    lr_model.orig <- glmnet(x=x_train, y=y_train, alpha=0.5, family="binomial")
    lr_pred.orig <- predict(lr_model.orig, newx=x_test, s=lambda.enet, type="response")
    lr_prauc.orig <- pr.curve(scores.class0=lr_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc.integral
    # enet_taxa <- coef.glmnet(lr_model.orig, s = lambda.enet)
    # enet_taxa <- enet_taxa@Dimnames[[1]][which(as.numeric(enet_taxa) > 0)]
    lr_rocauc.orig <- roc.curve(scores.class0=lr_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc
    
    return(list("prauc" = lr_prauc.orig, "rocauc" = lr_rocauc.orig))
  })
  
  condition <- as.factor(condition)
  cv_rf <- sapply(1:K, FUN=function(k) {
    test_idx <- which(CV_split_idx == k)
    x_train <- input_data[-test_idx,]
    x_test <- input_data[test_idx,]
    y_train <- condition[-test_idx]
    y_test <- condition[test_idx]
    ## original sample ratio
    rf_model.orig <- randomForest(x=x_train, y=y_train, xtest=x_test)
    rf_pred.orig <- rf_model.orig$test$votes[,2]
    rf_prauc.orig <- pr.curve(scores.class0=rf_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc.integral
    rf_rocauc.orig <- roc.curve(scores.class0=rf_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc
    
    return(list("prauc" = rf_prauc.orig, "rocauc" = rf_rocauc.orig))
  })
  
  cv_xgboost <- sapply(1:K, FUN=function(k) {
    test_idx <- which(CV_split_idx == k)
    x_train <- input_data[-test_idx,]
    x_test <- input_data[test_idx,]
    y_train <- condition[-test_idx]
    y_test <- condition[test_idx]
    ## original sample ratio
    model.orig <- xgboost(data=x_train, label=as.numeric(y_train)-1, max.depth=4, eta=1, nthread=detectCores()-1, nrounds=7, objective="binary:logistic")
    pred.orig <- predict(model.orig, newdata=x_test)
    prauc.orig <- pr.curve(scores.class0=pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc.integral
    rocauc.orig <- roc.curve(scores.class0=pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc
    
    return(list("prauc" = prauc.orig, "rocauc" = rocauc.orig))
  })
  
  ## Scale the feature matrix
  input_data <- scale(input_data)
  if(sum(is.na(colSums(input_data))) != 0){
    input_data <- input_data[,-which(is.na(colSums(input_data)))]
  }
  ## support vector machines with linear kernel
  cv_svm.linear <- sapply(1:K, FUN=function(k) {
    test_idx <- which(CV_split_idx == k)
    x_train <- input_data[-test_idx,]
    x_test <- input_data[test_idx,]
    y_train <- condition[-test_idx]
    y_test <- condition[test_idx]
    ## original sample ratio
    
    svm_model.orig <- svm(x=x_train, y=y_train, kernel="linear", probability=T)
    svm_pred.orig <- predict(svm_model.orig, newdata=x_test, probability=T)
    svm_pred.orig <- attr(svm_pred.orig, "probabilities")[,2]
    svm_prauc.orig <- pr.curve(scores.class0=svm_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc.integral
    svm_rocauc.orig <- roc.curve(scores.class0=svm_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc
    
    return(list("prauc" = svm_prauc.orig, "rocauc" = svm_rocauc.orig))
  })
  
  cv_svm.gauss <- sapply(1:K, FUN=function(k) {
    test_idx <- which(CV_split_idx == k)
    x_train <- input_data[-test_idx,]
    x_test <- input_data[test_idx,]
    y_train <- condition[-test_idx]
    y_test <- condition[test_idx]
    ## original sample ratio
    
    svm_model.orig <- svm(x=x_train, y=y_train, probability=T)
    svm_pred.orig <- predict(svm_model.orig, newdata=x_test, probability=T)
    svm_pred.orig <- attr(svm_pred.orig, "probabilities")[,2]
    svm_prauc.orig <- pr.curve(scores.class0=svm_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc.integral
    svm_rocauc.orig <- roc.curve(scores.class0=svm_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc
    
    return(list("prauc" = svm_prauc.orig, "rocauc" = svm_rocauc.orig))
  })
  return(list("logistic.prauc" = mean(unlist(cv_lr[1,])), "logistic.rocauc" = mean(unlist(cv_lr[2,])), "elastic_net.prauc" = mean(unlist(cv_lr.enet[1,])), "elastic_net.rocauc" = mean(unlist(cv_lr.enet[2,])), "random_forest.prauc" = mean(unlist(cv_rf[1,])), "random_forest.rocauc" = mean(unlist(cv_rf[2,])), "xgboost.prauc" = mean(unlist(cv_xgboost[1,])), "xgboost.rocauc" = mean(unlist(cv_xgboost[2,])), "svm.lr.prauc" = mean(unlist(cv_svm.linear[1,])), "svm.lr.rocauc" = mean(unlist(cv_svm.linear[2,])), "svm.gauss.prauc" =  mean(unlist(cv_svm.gauss[1,])), "svm.gauss.rocauc" = mean(unlist(cv_svm.gauss[2,]))))
}


matrix_ori <- as.matrix(Zeller_raw)
features_imputed <- Zeller_imp_result
imputed_results <- matrix(nrow = 8, ncol = 12)
for(i in 1:8){
  print(i)
  combined_conds <- union(features_imputed[[i]]$T2D_leq_cond, features_imputed[[i]]$T2D_grt_cond)
  print(length(combined_conds))
  if(length(combined_conds) > 1){
    overlap_taxa <- intersect(combined_conds, colnames(matrix_ori))
    print(length(overlap_taxa))
    if(length(overlap_taxa) != 0){
      input_data <- as.matrix(matrix_ori[,overlap_taxa])
      imputed_results[i,] <- unlist(classification_results(input_data = as.matrix(Zeller_raw), condition = Zeller_condition))}}}



full_data_result <- classification_results(input_data = as.matrix(Zeller_data), condition = condition)
classification_archive <- function(input_data, condition, features_raw, features_imputed){
  set.seed(2020)
  matrix_ori <- as.matrix(input_data)
  full_data_result <- classification_results(input_data = matrix_ori, condition = condition)
  raw_results <- matrix(nrow = 7, ncol = 12)
  for(i in 1:length(features_raw)){
    print(i)
    if(length(union(features_raw[[i]]$T2D_leq_cond, features_raw[[i]]$T2D_grt_cond)) != 0){
      input_data <- as.matrix(matrix_ori[,union(features_raw[[i]]$T2D_leq_cond, features_raw[[i]]$T2D_grt_cond)])
      raw_results[i,] <- unlist(classification_results(input_data = input_data, condition = condition))
    }else{
    }
  }
  rownames(raw_results) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2")
  colnames(raw_results) <- c("logistic_reg.prauc", "logistic_reg.rocauc", "elastic_net.prauc", "elastic_net.rocauc", "random_forest.prauc", "random_forest.rocauc", "xgboost.prauc", "xgboost.rocauc", "svm.lr.prauc", "svm.lr.rocauc", "svm.guass.prauc", "svm.guass.rocauc")
  
  imputed_results <- matrix(nrow = 7, ncol = 12)
  for(i in 1:length(features_imputed)){
    print(i)
    if(length(union(features_imputed[[i]]$T2D_leq_cond, features_imputed[[i]]$T2D_grt_cond)) != 0){
      input_data <- as.matrix(matrix_ori[,union(features_imputed[[i]]$T2D_leq_cond, features_imputed[[i]]$T2D_grt_cond)])
      imputed_results[i,] <- unlist(classification_results(input_data = input_data, condition = condition))
    }else{
    }
  }
  rownames(imputed_results) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2")
  colnames(imputed_results) <- c("logistic_reg.prauc", "logistic_reg.rocauc", "elastic_net.prauc", "elastic_net.rocauc", "random_forest.prauc", "random_forest.rocauc", "xgboost.prauc", "xgboost.rocauc", "svm.lr.prauc", "svm.lr.rocauc", "svm.guass.prauc", "svm.guass.rocauc")
  
  # rbind(unlist(full_data_result), raw_results, imputed_results)
  return(list("full_data_result" = unlist(full_data_result), "raw_result" = raw_results, "imputed_result" = imputed_results))
}

classification <- function(input_data, condition, features_raw, features_imputed){
  set.seed(2020)
  matrix_ori <- as.matrix(input_data)
  full_data_result <- classification_results(input_data = matrix_ori, condition = condition)
  raw_results <- matrix(nrow = 8, ncol = 12)
  for(i in 1:8){
    print(i)
    if(length(union(features_raw[[i]]$T2D_leq_cond, features_raw[[i]]$T2D_grt_cond)) > 1){
      overlap_taxa <- intersect(union(features_raw[[i]]$T2D_leq_cond, features_raw[[i]]$T2D_grt_cond), colnames(matrix_ori))
      if(length(overlap_taxa) != 0){
        input_data <- as.matrix(matrix_ori[,overlap_taxa])
        raw_results[i,] <- unlist(classification_results(input_data = input_data, condition = condition))
      }
    }else{
    }
  }
  rownames(raw_results) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2", "Omnibus")
  colnames(raw_results) <- c("logistic_reg.prauc", "logistic_reg.rocauc", "elastic_net.prauc", "elastic_net.rocauc", "random_forest.prauc", "random_forest.rocauc", "xgboost.prauc", "xgboost.rocauc", "svm.lr.prauc", "svm.lr.rocauc", "svm.guass.prauc", "svm.guass.rocauc")
  
  imputed_results <- matrix(nrow = 8, ncol = 12)
  for(i in 1:8){
    print(i)
    if(length(union(features_imputed[[i]]$T2D_leq_cond, features_imputed[[i]]$T2D_grt_cond)) > 1){
      overlap_taxa <- intersect(union(features_imputed[[i]]$T2D_leq_cond, features_imputed[[i]]$T2D_grt_cond), colnames(matrix_ori))
      if(length(overlap_taxa) != 0){
        input_data <- as.matrix(matrix_ori[,overlap_taxa])
        results <- classification_results(input_data = input_data, condition = condition)
        if (any(is.na(unlist(results)))) {
          print(paste("在迭代", i, "中发现 NA 值"))
        } else {
          imputed_results[i,] <- unlist(results)
        }
        #imputed_results[i,] <- unlist(classification_results(input_data = input_data, condition = condition))
      }
    }else{
    }
  }
  rownames(imputed_results) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2", "Omnibus")
  colnames(imputed_results) <- c("logistic_reg.prauc", "logistic_reg.rocauc", "elastic_net.prauc", "elastic_net.rocauc", "random_forest.prauc", "random_forest.rocauc", "xgboost.prauc", "xgboost.rocauc", "svm.lr.prauc", "svm.lr.rocauc", "svm.guass.prauc", "svm.guass.rocauc")
  
  # rbind(unlist(full_data_result), raw_results, imputed_results)
  return(list("full_data_result" = unlist(full_data_result), "raw_result" = raw_results, "imputed_result" = imputed_results))
}


Zeller_condition[Zeller_condition == "control"]=0
Zeller_condition[Zeller_condition == "CRC"]=1
Zeller_condition <- as.numeric(Zeller_condition)
Zeller_result <- classification(input_data = as.matrix(Zeller_raw), condition = Zeller_condition, features_raw = Zeller_raw_result, features_imputed = Zeller_imp_result)

Yu_condition[Yu_condition == "control"]=0
Yu_condition[Yu_condition == "CRC"]=1
Yu_condition <- as.numeric(Yu_condition)
Yu_result <- classification(input_data = as.matrix(Yu_raw), condition = Yu_condition, features_raw = Yu_raw_result, features_imputed = Yu_imp_result)

Vogtmann_condition[Vogtmann_condition == "control"]=0
Vogtmann_condition[Vogtmann_condition == "CRC"]=1
Vogtmann_condition <- as.numeric(Vogtmann_condition)
Vogtmann_result <- classification(input_data = as.matrix(Vogtmann_raw), condition = Vogtmann_condition, features_raw = Vogtmann_raw_result, features_imputed = Vogtmann_imp_result)

Feng_condition[Feng_condition == "control"]=0
Feng_condition[Feng_condition == "CRC"]=1
Feng_condition <- as.numeric(Feng_condition)
Feng_result <- classification(input_data = as.matrix(Feng_raw), condition = Feng_condition, features_raw = Feng_raw_result, features_imputed = Feng_imp_result)

Karlsson_condition[Karlsson_condition == "control"]=0
Karlsson_condition[Karlsson_condition == "T2D"]=1
Karlsson_condition <- as.numeric(Karlsson_condition)
Karlsson_result <- classification(input_data = as.matrix(Karlsson_raw), condition = Karlsson_condition, features_raw = Karlsson_raw_result, features_imputed = Karlsson_imp_result)

Qin_condition[Qin_condition == "control"]=0
Qin_condition[Qin_condition == "T2D"]=1
Qin_condition <- as.numeric(Qin_condition)
Qin_result <- classification(input_data = as.matrix(Qin_raw), condition = Qin_condition, features_raw = Qin_raw_result, features_imputed = Qin_imp_result)



###DESeq2与Omnibus的p值分布图
png("Real data p value of DESeq2 distribution/Qin_raw.png", width = 480, height = 350)
hist(Qin_raw_result$DEseq2_pval, ylim = c(0, 90), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of DESeq2 distribution/Qin_imp.png", width = 480, height = 350)
hist(Qin_imp_result$DEseq2_pval, ylim = c(0, 90), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of DESeq2 distribution/Karlsson_raw.png", width = 480, height = 350)
hist(Karlsson_raw_result$DEseq2_pval, ylim = c(0, 60), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of DESeq2 distribution/Karlsson_imp.png", width = 480, height = 350)
hist(Karlsson_imp_result$DEseq2_pval, ylim = c(0, 60), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of DESeq2 distribution/Feng_raw.png", width = 480, height = 350)
hist(Feng_raw_result$DEseq2_pval, ylim = c(0, 90), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of DESeq2 distribution/Feng_imp.png", width = 480, height = 350)
hist(Feng_imp_result$DEseq2_pval, ylim = c(0, 90), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Yu_raw.png", width = 480, height = 350)
hist(Yu_raw_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Yu_imp.png", width = 480, height = 350)
hist(Yu_imp_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of DESeq2 distribution/Vogtmann_raw.png", width = 480, height = 350)
hist(Vogtmann_raw_result$DEseq2_pval, ylim = c(0, 90), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of DESeq2 distribution/Vogtmann_imp.png", width = 480, height = 350)
hist(Vogtmann_imp_result$DEseq2_pval, ylim = c(0, 90), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of DESeq2 distribution/Zeller_raw.png", width = 480, height = 350)
hist(Zeller_raw_result$DEseq2_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of DESeq2 distribution/Zeller_imp.png", width = 480, height = 350)
hist(Zeller_imp_result$DEseq2_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()


# Open a PDF device with the correct filename extension and size in inches
pdf("/Users/hanxinyu/Desktop/real_data_analysis/Real_data_p_value_DESeq2_distribution_combined_1.pdf", width = 20, height = 14)

# Set up the graphic parameters for a 2x3 plot layout
par(mfrow = c(2, 3), mar = c(4, 4.5, 2, 1))

# Create the histograms for the Qin dataset
hist(Qin_raw_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis = 1.5, main = "Qin_Raw", xlab = "P-value")
hist(Karlsson_raw_result$DEseq2_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis = 1.5, main = "Karlsson_Raw", xlab = "P-value")
hist(Feng_raw_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis = 1.5, main = "Feng_Raw", xlab = "P-value")
hist(Qin_imp_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis = 1.5, main = "Qin_Imp", xlab = "P-value")
hist(Karlsson_imp_result$DEseq2_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis = 1.5, main = "Karlsson_Imp", xlab = "P-value")
hist(Feng_imp_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis = 1.5, main = "Feng_Imp", xlab = "P-value")

dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Real_data_p_value_DESeq2_distribution_combined_2.pdf", width = 20, height = 14)

# Set up the graphic parameters for a 2x3 plot layout
par(mfrow = c(2, 3), mar = c(4, 4.5, 2, 1))

# Create the histograms for the Qin dataset
hist(Yu_raw_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis=1.5, main = "Yu_Raw", xlab = "P-value")
hist(Vogtmann_raw_result$DEseq2_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "Vogtmann_Raw", xlab = "P-value")
hist(Zeller_raw_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis=1.5, main = "Zeller_Raw", xlab = "P-value")
hist(Yu_imp_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis=1.5, main = "Yu_Imp", xlab = "P-value")
hist(Vogtmann_imp_result$DEseq2_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "Vogtmann_Imp", xlab = "P-value")
hist(Zeller_imp_result$DEseq2_pval, ylim = c(0, 120), cex.lab = 1.5, cex.axis=1.5, main = "Zeller_Imp", xlab = "P-value")

dev.off()


# Open a PDF device with the correct filename extension and size in inches
pdf("/Users/hanxinyu/Desktop/real_data_analysis/Real_data_p_value_Omnibus_distribution_combined_1.pdf", width = 20, height = 14)

# Set up the graphic parameters for a 2x3 plot layout
par(mfrow = c(2, 3), mar = c(4, 4.5, 2, 1))

# Create the histograms for the Qin dataset
hist(Qin_raw_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis = 1.5, main = "Qin_Raw", xlab = "P-value")
hist(Karlsson_raw_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis = 1.5, main = "Karlsson_Raw", xlab = "P-value")
hist(Feng_raw_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis = 1.5, main = "Feng_Raw", xlab = "P-value")
hist(Qin_imp_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis = 1.5, main = "Qin_Imp", xlab = "P-value")
hist(Karlsson_imp_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis = 1.5, main = "Karlsson_Imp", xlab = "P-value")
hist(Feng_imp_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis = 1.5, main = "Feng_Imp", xlab = "P-value")

dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Real_data_p_value_Omnibus_distribution_combined_2.pdf", width = 20, height = 14)

# Set up the graphic parameters for a 2x3 plot layout
par(mfrow = c(2, 3), mar = c(4, 4.5, 2, 1))

# Create the histograms for the Qin dataset
hist(Yu_raw_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "Yu_Raw", xlab = "P-value")
hist(Vogtmann_raw_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "Vogtmann_Raw", xlab = "P-value")
hist(Zeller_raw_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "Zeller_Raw", xlab = "P-value")
hist(Yu_imp_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "Yu_Imp", xlab = "P-value")
hist(Vogtmann_imp_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "Vogtmann_Imp", xlab = "P-value")
hist(Zeller_imp_result$Omnibus_pval, ylim = c(0, 100), cex.lab = 1.5, cex.axis=1.5, main = "Zeller_Imp", xlab = "P-value")

dev.off()





png("Real data p value of Omnibus distribution/Qin_raw.png", width = 480, height = 350)
hist(Qin_raw_result$Omnibus_pval, ylim = c(0, 60), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Qin_imp.png", width = 480, height = 350)
hist(Qin_imp_result$Omnibus_pval, ylim = c(0, 60), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Karlsson_raw.png", width = 480, height = 350)
hist(Karlsson_raw_result$Omnibus_pval, ylim = c(0, 40), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Karlsson_imp.png", width = 480, height = 350)
hist(Karlsson_imp_result$Omnibus_pval, ylim = c(0, 40), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Feng_raw.png", width = 480, height = 350)
hist(Feng_raw_result$Omnibus_pval, ylim = c(0, 70), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Feng_imp.png", width = 480, height = 350)
hist(Feng_imp_result$Omnibus_pval, ylim = c(0, 70), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Yu_raw.png", width = 480, height = 350)
hist(Yu_raw_result$Omnibus_pval, ylim = c(0, 70), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Yu_imp.png", width = 480, height = 350)
hist(Yu_imp_result$Omnibus_pval, ylim = c(0, 70), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Vogtmann_raw.png", width = 480, height = 350)
hist(Vogtmann_raw_result$Omnibus_pval, ylim = c(0, 40), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Vogtmann_imp.png", width = 480, height = 350)
hist(Vogtmann_imp_result$Omnibus_pval, ylim = c(0, 40), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Zeller_raw.png", width = 480, height = 350)
hist(Zeller_raw_result$Omnibus_pval, ylim = c(0, 70), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()
png("Real data p value of Omnibus distribution/Zeller_imp.png", width = 480, height = 350)
hist(Zeller_imp_result$Omnibus_pval, ylim = c(0, 70), cex.lab = 1.5, cex.axis=1.5, main = "", xlab = "")
dev.off()




raw_imp_DEseq2_result <- c(Qin_result$raw_result[7,9], Qin_result$imputed_result[7,9], Karlsson_result$raw_result[7,9], Karlsson_result$imputed_result[7,9],
                           Feng_result$raw_result[7,9], Feng_result$imputed_result[7,9], Yu_result$raw_result[7,9], Yu_result$imputed_result[7,9],
                           Vogtmann_result$raw_result[7,9], Vogtmann_result$imputed_result[7,9], Zeller_result$raw_result[7,9], Zeller_result$imputed_result[7,9])
raw_imp_DEseq2_result <- cbind(raw_imp_DEseq2_result, c("Qin", "Qin", "Karlsson", "Karlsson", "Feng", "Feng", "Yu", "Yu", "Vogtmann", "Vogtmann", "Zeller", "Zeller"))
raw_imp_DEseq2_result <- cbind(raw_imp_DEseq2_result, c("raw", "imp", "raw", "imp", "raw", "imp", "raw", "imp", "raw", "imp", "raw", "imp"))
colnames(raw_imp_DEseq2_result) <- c("values", "data", "type")
raw_imp_DEseq2_result <- data.frame(raw_imp_DEseq2_result)
raw_imp_DEseq2_result$values <- as.numeric(as.character(raw_imp_DEseq2_result$values))
raw_imp_DEseq2_result$data <- factor(raw_imp_DEseq2_result$data, levels = c("Qin", "Karlsson", "Feng", "Yu", "Vogtmann", "Zeller"))
raw_imp_DEseq2_result$type <- factor(raw_imp_DEseq2_result$type, levels = c("raw", "imp"))


pdf("/Users/hanxinyu/Desktop/real_data_analysis/DEseq2_result11.pdf", width = 15, height = 5)
ggplot(raw_imp_DEseq2_result, aes(fill=type, y=values, x=data)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_x_discrete(labels=c("Qin et al.", "Karlsson et al.", "Feng et al.", 
                            "Yu et al.", "Vogtmann et al.", "Zeller et al.")) +
  scale_fill_manual(name="Method",
                    labels=c("DESeq2-phyloseq", "TphPMF+DESeq2-phyloseq"),
                    values=c("raw"=alpha("#636363", 0.5), "imp"="#636363")) +
  ylab("PR-AUC") +
  coord_cartesian(ylim=c(0.2, 0.8)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=12, angle=45, hjust=1),
        axis.title.y=element_text(size=14),
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        legend.position="right") +
  guides(fill=guide_legend(override.aes=list(alpha=c(0.5,1))))


dev.off()



raw_imp_DEseq2_result <- c(Qin_result$raw_result[7,10], Qin_result$imputed_result[7,10], Karlsson_result$raw_result[7,10], Karlsson_result$imputed_result[7,10],
                           Feng_result$raw_result[7,10], Feng_result$imputed_result[7,10], Yu_result$raw_result[7,10], Yu_result$imputed_result[7,10],
                           Vogtmann_result$raw_result[7,10], Vogtmann_result$imputed_result[7,10], Zeller_result$raw_result[7,10], Zeller_result$imputed_result[7,10])
raw_imp_DEseq2_result <- cbind(raw_imp_DEseq2_result, c("Qin", "Qin", "Karlsson", "Karlsson", "Feng", "Feng", "Yu", "Yu", "Vogtmann", "Vogtmann", "Zeller", "Zeller"))
raw_imp_DEseq2_result <- cbind(raw_imp_DEseq2_result, c("raw", "imp", "raw", "imp", "raw", "imp", "raw", "imp", "raw", "imp", "raw", "imp"))
colnames(raw_imp_DEseq2_result) <- c("values", "data", "type")
raw_imp_DEseq2_result <- data.frame(raw_imp_DEseq2_result)
raw_imp_DEseq2_result$values <- as.numeric(as.character(raw_imp_DEseq2_result$values))
raw_imp_DEseq2_result$data <- factor(raw_imp_DEseq2_result$data, levels = c("Qin", "Karlsson", "Feng", "Yu", "Vogtmann", "Zeller"))
raw_imp_DEseq2_result$type <- factor(raw_imp_DEseq2_result$type, levels = c("raw", "imp"))


pdf("/Users/hanxinyu/Desktop/real_data_analysis/DEseq2_result2.pdf", width = 15, height = 5)
ggplot(raw_imp_DEseq2_result, aes(fill=type, y=values, x=data)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_x_discrete(labels=c("Qin et al.", "Karlsson et al.", "Feng et al.", 
                            "Yu et al.", "Vogtmann et al.", "Zeller et al.")) +
  scale_fill_manual(name="Method",
                    labels=c("DESeq2-phyloseq", "BHPMF+DESeq2-phyloseq"),
                    values=c("raw"=alpha("#636363", 0.5), "imp"="#636363")) +
  ylab("ROC-AUC") +
  coord_cartesian(ylim=c(0, 0.8)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=12, angle=45, hjust=1),
        axis.title.y=element_text(size=14),
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        legend.position="right") +
  guides(fill=guide_legend(override.aes=list(alpha=c(0.5,1))))


dev.off()






raw_imp_Omnibus_result <- c(Qin_result$raw_result[8,12], Qin_result$imputed_result[8,12], Karlsson_result$raw_result[8,12], Karlsson_result$imputed_result[8,12],
                            Feng_result$raw_result[8,12], Feng_result$imputed_result[8,12], Yu_result$raw_result[8,12], Yu_result$imputed_result[8,12],
                            Vogtmann_result$raw_result[8,12], Vogtmann_result$imputed_result[8,12], Zeller_result$raw_result[8,12], Zeller_result$imputed_result[8,12])
raw_imp_Omnibus_result <- cbind(raw_imp_Omnibus_result, c("Qin", "Qin", "Karlsson", "Karlsson", "Feng", "Feng", "Yu", "Yu", "Vogtmann", "Vogtmann", "Zeller", "Zeller"))
raw_imp_Omnibus_result <- cbind(raw_imp_Omnibus_result, c("raw", "imp", "raw", "imp", "raw", "imp", "raw", "imp", "raw", "imp", "raw", "imp"))
colnames(raw_imp_Omnibus_result) <- c("values", "data", "type")
raw_imp_Omnibus_result <- data.frame(raw_imp_Omnibus_result)
raw_imp_Omnibus_result$values <- as.numeric(as.character(raw_imp_Omnibus_result$values))
raw_imp_Omnibus_result$data <- factor(raw_imp_Omnibus_result$data, levels = c("Qin", "Karlsson", "Feng", "Yu", "Vogtmann", "Zeller"))
raw_imp_Omnibus_result$type <- factor(raw_imp_Omnibus_result$type, levels = c("raw", "imp"))


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Omnibus_result_svm.guass_2.pdf", width = 15, height = 5)
ggplot(raw_imp_Omnibus_result, aes(fill=type, y=values, x=data)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_x_discrete(labels=c("Qin et al.", "Karlsson et al.", "Feng et al.", 
                            "Yu et al.", "Vogtmann et al.", "Zeller et al.")) +
  scale_fill_manual(name="Method",
                    labels=c("Omnibus", "BHPMF+Omnibus"),
                    values=c("raw"=alpha("#636363", 0.5), "imp"="#636363")) +
  ylab("ROC-AUC") +
  coord_cartesian(ylim=c(0, 0.8)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=12, angle=45, hjust=1),
        axis.title.y=element_text(size=14),
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        legend.position="right") +
  guides(fill=guide_legend(override.aes=list(alpha=c(0.5,1))))


dev.off()









T2D_grt_cond_overlap <- matrix(nrow = 8, ncol = 8)
for(i in 1:8){
  T2D_grt_cond_overlap[i,1] <- length(Qin_raw_result[[i]]$T2D_grt_cond)
  T2D_grt_cond_overlap[i,2] <- length(Qin_imp_result[[i]]$T2D_grt_cond)
  T2D_grt_cond_overlap[i,3] <- length(Karlsson_raw_result[[i]]$T2D_grt_cond)
  T2D_grt_cond_overlap[i,4] <- length(Karlsson_imp_result[[i]]$T2D_grt_cond)
  T2D_grt_cond_overlap[i,5] <- length(intersect(Qin_raw_result[[i]]$T2D_grt_cond, Karlsson_raw_result[[i]]$T2D_grt_cond))
  T2D_grt_cond_overlap[i,6] <- length(intersect(Qin_imp_result[[i]]$T2D_grt_cond, Karlsson_imp_result[[i]]$T2D_grt_cond))
  T2D_grt_cond_overlap[i,7] <- T2D_grt_cond_overlap[i,5] / length(union(Qin_raw_result[[i]]$T2D_grt_cond, Karlsson_raw_result[[i]]$T2D_grt_cond))
  T2D_grt_cond_overlap[i,8] <- T2D_grt_cond_overlap[i,6] / length(union(Qin_imp_result[[i]]$T2D_grt_cond, Karlsson_imp_result[[i]]$T2D_grt_cond))
}
colnames(T2D_grt_cond_overlap) <- c("Qin_raw", "Qin_imp", "Karlsson_raw", "Karlsson_imp", "intersect_raw", "intersect_imp", "DR_raw", "DR_imp")
rownames(T2D_grt_cond_overlap) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2","Omnibus")


T2D_leq_cond_overlap <- matrix(nrow = 8, ncol = 8)
for(i in 1:8){
  T2D_leq_cond_overlap[i,1] <- length(Qin_raw_result[[i]]$T2D_leq_cond)
  T2D_leq_cond_overlap[i,2] <- length(Qin_imp_result[[i]]$T2D_leq_cond)
  T2D_leq_cond_overlap[i,3] <- length(Karlsson_raw_result[[i]]$T2D_leq_cond)
  T2D_leq_cond_overlap[i,4] <- length(Karlsson_imp_result[[i]]$T2D_leq_cond)
  T2D_leq_cond_overlap[i,5] <- length(intersect(Qin_raw_result[[i]]$T2D_leq_cond, Karlsson_raw_result[[i]]$T2D_leq_cond))
  T2D_leq_cond_overlap[i,6] <- length(intersect(Qin_imp_result[[i]]$T2D_leq_cond, Karlsson_imp_result[[i]]$T2D_leq_cond))
  T2D_leq_cond_overlap[i,7] <- T2D_leq_cond_overlap[i,5] / length(union(Qin_raw_result[[i]]$T2D_leq_cond, Karlsson_raw_result[[i]]$T2D_leq_cond))
  T2D_leq_cond_overlap[i,8] <- T2D_leq_cond_overlap[i,6] / length(union(Qin_imp_result[[i]]$T2D_leq_cond, Karlsson_imp_result[[i]]$T2D_leq_cond))
}
colnames(T2D_leq_cond_overlap) <- c("Qin_raw", "Qin_imp", "Karlsson_raw", "Karlsson_imp", "intersect_raw", "intersect_imp", "DR_raw", "DR_imp")
rownames(T2D_leq_cond_overlap) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2","Omnibus")






CRC_grt_overlap_raw <- matrix(nrow = 8, ncol = 26)
for(i in 1:8){
  CRC_grt_overlap_raw[i,1] <- length(Feng_raw_result[[i]]$T2D_grt_cond)
  CRC_grt_overlap_raw[i,2] <- length(Vogtmann_raw_result[[i]]$T2D_grt_cond)
  CRC_grt_overlap_raw[i,3] <- length(Yu_raw_result[[i]]$T2D_grt_cond)
  CRC_grt_overlap_raw[i,4] <- length(Zeller_raw_result[[i]]$T2D_grt_cond)
  
  all_intersect <- intersect(intersect(intersect(Vogtmann_raw_result[[i]]$T2D_grt_cond, Yu_raw_result[[i]]$T2D_grt_cond),
                                       Zeller_raw_result[[i]]$T2D_grt_cond),Feng_raw_result[[i]]$T2D_grt_cond)
  all_union <- union(union(union(Feng_raw_result[[i]]$T2D_grt_cond, Vogtmann_raw_result[[i]]$T2D_grt_cond), 
                           Yu_raw_result[[i]]$T2D_grt_cond),  Zeller_raw_result[[i]]$T2D_grt_cond)
  FV <- intersect(Feng_raw_result[[i]]$T2D_grt_cond, Vogtmann_raw_result[[i]]$T2D_grt_cond)
  FY <- intersect(Feng_raw_result[[i]]$T2D_grt_cond, Yu_raw_result[[i]]$T2D_grt_cond)
  FZ <- intersect(Feng_raw_result[[i]]$T2D_grt_cond, Zeller_raw_result[[i]]$T2D_grt_cond)
  VY <- intersect(Vogtmann_raw_result[[i]]$T2D_grt_cond, Yu_raw_result[[i]]$T2D_grt_cond)
  VZ <- intersect(Vogtmann_raw_result[[i]]$T2D_grt_cond, Zeller_raw_result[[i]]$T2D_grt_cond)
  YZ <- intersect(Yu_raw_result[[i]]$T2D_grt_cond, Zeller_raw_result[[i]]$T2D_grt_cond)
  FVY <- intersect(intersect(Feng_raw_result[[i]]$T2D_grt_cond, Vogtmann_raw_result[[i]]$T2D_grt_cond),
                   Yu_raw_result[[i]]$T2D_grt_cond)
  FVZ <- intersect(intersect(Feng_raw_result[[i]]$T2D_grt_cond, Vogtmann_raw_result[[i]]$T2D_grt_cond),
                   Zeller_raw_result[[i]]$T2D_grt_cond)
  FYZ <- intersect(intersect(Feng_raw_result[[i]]$T2D_grt_cond, Yu_raw_result[[i]]$T2D_grt_cond),
                   Zeller_raw_result[[i]]$T2D_grt_cond)
  VYZ <- intersect(intersect(Vogtmann_raw_result[[i]]$T2D_grt_cond, Yu_raw_result[[i]]$T2D_grt_cond),
                   Zeller_raw_result[[i]]$T2D_grt_cond)
  
  CRC_grt_overlap_raw[i,15] <- length(all_intersect)
  CRC_grt_overlap_raw[i,11] <- length(FVY[!(FVY %in% all_intersect)])
  CRC_grt_overlap_raw[i,12] <- length(FVZ[!(FVZ %in% all_intersect)])
  CRC_grt_overlap_raw[i,13] <- length(FYZ[!(FYZ %in% all_intersect)])
  CRC_grt_overlap_raw[i,14] <- length(VYZ[!(VYZ %in% all_intersect)])
  
  CRC_grt_overlap_raw[i,5] <- length(FV[!(FV %in% union(FVY, FVZ))])
  CRC_grt_overlap_raw[i,6] <- length(FY[!(FY %in% union(FVY, FYZ))])
  CRC_grt_overlap_raw[i,7] <- length(FZ[!(FZ %in% union(FVZ, FYZ))])
  CRC_grt_overlap_raw[i,8] <- length(VY[!(VY %in% union(FVY, VYZ))])
  CRC_grt_overlap_raw[i,9] <- length(VZ[!(VZ %in% union(FVZ, VYZ))])
  CRC_grt_overlap_raw[i,10] <- length(YZ[!(YZ %in% union(FYZ, VYZ))])
  
  FV <- union(Feng_raw_result[[i]]$T2D_grt_cond, Vogtmann_raw_result[[i]]$T2D_grt_cond)
  FY <- union(Feng_raw_result[[i]]$T2D_grt_cond, Yu_raw_result[[i]]$T2D_grt_cond)
  FZ <- union(Feng_raw_result[[i]]$T2D_grt_cond, Zeller_raw_result[[i]]$T2D_grt_cond)
  VY <- union(Vogtmann_raw_result[[i]]$T2D_grt_cond, Yu_raw_result[[i]]$T2D_grt_cond)
  VZ <- union(Vogtmann_raw_result[[i]]$T2D_grt_cond, Zeller_raw_result[[i]]$T2D_grt_cond)
  YZ <- union(Yu_raw_result[[i]]$T2D_grt_cond, Zeller_raw_result[[i]]$T2D_grt_cond)
  FVY <- union(union(Feng_raw_result[[i]]$T2D_grt_cond, Vogtmann_raw_result[[i]]$T2D_grt_cond),
               Yu_raw_result[[i]]$T2D_grt_cond)
  FVZ <- union(union(Feng_raw_result[[i]]$T2D_grt_cond, Vogtmann_raw_result[[i]]$T2D_grt_cond),
               Zeller_raw_result[[i]]$T2D_grt_cond)
  FYZ <- union(union(Feng_raw_result[[i]]$T2D_grt_cond, Yu_raw_result[[i]]$T2D_grt_cond),
               Zeller_raw_result[[i]]$T2D_grt_cond)
  VYZ <- union(union(Vogtmann_raw_result[[i]]$T2D_grt_cond, Yu_raw_result[[i]]$T2D_grt_cond),
               Zeller_raw_result[[i]]$T2D_grt_cond)
  
  CRC_grt_overlap_raw[i,16] <- CRC_grt_overlap_raw[i,5] / length(FV)
  CRC_grt_overlap_raw[i,17] <- CRC_grt_overlap_raw[i,6] / length(FY)
  CRC_grt_overlap_raw[i,18] <- CRC_grt_overlap_raw[i,7] / length(FZ)
  CRC_grt_overlap_raw[i,19] <- CRC_grt_overlap_raw[i,8] / length(VY)
  CRC_grt_overlap_raw[i,20] <- CRC_grt_overlap_raw[i,9] / length(VZ)
  CRC_grt_overlap_raw[i,21] <- CRC_grt_overlap_raw[i,10] / length(YZ)
  
  CRC_grt_overlap_raw[i,22] <- CRC_grt_overlap_raw[i,11] / length(FVY)
  CRC_grt_overlap_raw[i,23] <- CRC_grt_overlap_raw[i,12] / length(FVZ)
  CRC_grt_overlap_raw[i,24] <- CRC_grt_overlap_raw[i,13] / length(FYZ)
  CRC_grt_overlap_raw[i,25] <- CRC_grt_overlap_raw[i,14] / length(VYZ)
  CRC_grt_overlap_raw[i,26] <- CRC_grt_overlap_raw[i,15] / length(all_union)
}
rownames(CRC_grt_overlap_raw) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2","Omnibus")
colnames(CRC_grt_overlap_raw) <- c("Feng", "Vogtmann", "Yu", "Zeller", "FV", "FY", "FZ", "VY", "VZ", "YZ", "FVY", "FVZ", "FYZ", "VYZ", "FVYZ", "FV_rate", "FY_rate", "FZ_rate", "VY_rate", "VZ_rate", "YZ_rate", "FVY_rate", "FVZ_rate", "FYZ_rate", "VYZ_rate", "FVYZ_rate")

CRC_grt_overlap_imputed <- matrix(nrow = 8, ncol = 26)
for(i in 1:8){
  CRC_grt_overlap_imputed[i,1] <- length(Feng_imp_result[[i]]$T2D_grt_cond)
  CRC_grt_overlap_imputed[i,2] <- length(Vogtmann_imp_result[[i]]$T2D_grt_cond)
  CRC_grt_overlap_imputed[i,3] <- length(Yu_imp_result[[i]]$T2D_grt_cond)
  CRC_grt_overlap_imputed[i,4] <- length(Zeller_imp_result[[i]]$T2D_grt_cond)
  
  all_intersect <- intersect(intersect(intersect(Vogtmann_imp_result[[i]]$T2D_grt_cond, Yu_imp_result[[i]]$T2D_grt_cond),
                                       Zeller_imp_result[[i]]$T2D_grt_cond),Feng_imp_result[[i]]$T2D_grt_cond)
  all_union <- union(union(union(Feng_imp_result[[i]]$T2D_grt_cond, Vogtmann_imp_result[[i]]$T2D_grt_cond), 
                           Yu_imp_result[[i]]$T2D_grt_cond),  Zeller_imp_result[[i]]$T2D_grt_cond)
  FV <- intersect(Feng_imp_result[[i]]$T2D_grt_cond, Vogtmann_imp_result[[i]]$T2D_grt_cond)
  FY <- intersect(Feng_imp_result[[i]]$T2D_grt_cond, Yu_imp_result[[i]]$T2D_grt_cond)
  FZ <- intersect(Feng_imp_result[[i]]$T2D_grt_cond, Zeller_imp_result[[i]]$T2D_grt_cond)
  VY <- intersect(Vogtmann_imp_result[[i]]$T2D_grt_cond, Yu_imp_result[[i]]$T2D_grt_cond)
  VZ <- intersect(Vogtmann_imp_result[[i]]$T2D_grt_cond, Zeller_imp_result[[i]]$T2D_grt_cond)
  YZ <- intersect(Yu_imp_result[[i]]$T2D_grt_cond, Zeller_imp_result[[i]]$T2D_grt_cond)
  FVY <- intersect(intersect(Feng_imp_result[[i]]$T2D_grt_cond, Vogtmann_imp_result[[i]]$T2D_grt_cond),
                   Yu_imp_result[[i]]$T2D_grt_cond)
  FVZ <- intersect(intersect(Feng_imp_result[[i]]$T2D_grt_cond, Vogtmann_imp_result[[i]]$T2D_grt_cond),
                   Zeller_imp_result[[i]]$T2D_grt_cond)
  FYZ <- intersect(intersect(Feng_imp_result[[i]]$T2D_grt_cond, Yu_imp_result[[i]]$T2D_grt_cond),
                   Zeller_imp_result[[i]]$T2D_grt_cond)
  VYZ <- intersect(intersect(Vogtmann_imp_result[[i]]$T2D_grt_cond, Yu_imp_result[[i]]$T2D_grt_cond),
                   Zeller_imp_result[[i]]$T2D_grt_cond)
  
  CRC_grt_overlap_imputed[i,15] <- length(all_intersect)
  CRC_grt_overlap_imputed[i,11] <- length(FVY[!(FVY %in% all_intersect)])
  CRC_grt_overlap_imputed[i,12] <- length(FVZ[!(FVZ %in% all_intersect)])
  CRC_grt_overlap_imputed[i,13] <- length(FYZ[!(FYZ %in% all_intersect)])
  CRC_grt_overlap_imputed[i,14] <- length(VYZ[!(VYZ %in% all_intersect)])
  
  CRC_grt_overlap_imputed[i,5] <- length(FV[!(FV %in% union(FVY, FVZ))])
  CRC_grt_overlap_imputed[i,6] <- length(FY[!(FY %in% union(FVY, FYZ))])
  CRC_grt_overlap_imputed[i,7] <- length(FZ[!(FZ %in% union(FVZ, FYZ))])
  CRC_grt_overlap_imputed[i,8] <- length(VY[!(VY %in% union(FVY, VYZ))])
  CRC_grt_overlap_imputed[i,9] <- length(VZ[!(VZ %in% union(FVZ, VYZ))])
  CRC_grt_overlap_imputed[i,10] <- length(YZ[!(YZ %in% union(FYZ, VYZ))])
  
  FV <- union(Feng_imp_result[[i]]$T2D_grt_cond, Vogtmann_imp_result[[i]]$T2D_grt_cond)
  FY <- union(Feng_imp_result[[i]]$T2D_grt_cond, Yu_imp_result[[i]]$T2D_grt_cond)
  FZ <- union(Feng_imp_result[[i]]$T2D_grt_cond, Zeller_imp_result[[i]]$T2D_grt_cond)
  VY <- union(Vogtmann_imp_result[[i]]$T2D_grt_cond, Yu_imp_result[[i]]$T2D_grt_cond)
  VZ <- union(Vogtmann_imp_result[[i]]$T2D_grt_cond, Zeller_imp_result[[i]]$T2D_grt_cond)
  YZ <- union(Yu_imp_result[[i]]$T2D_grt_cond, Zeller_imp_result[[i]]$T2D_grt_cond)
  FVY <- union(union(Feng_imp_result[[i]]$T2D_grt_cond, Vogtmann_imp_result[[i]]$T2D_grt_cond),
               Yu_imp_result[[i]]$T2D_grt_cond)
  FVZ <- union(union(Feng_imp_result[[i]]$T2D_grt_cond, Vogtmann_imp_result[[i]]$T2D_grt_cond),
               Zeller_imp_result[[i]]$T2D_grt_cond)
  FYZ <- union(union(Feng_imp_result[[i]]$T2D_grt_cond, Yu_imp_result[[i]]$T2D_grt_cond),
               Zeller_imp_result[[i]]$T2D_grt_cond)
  VYZ <- union(union(Vogtmann_imp_result[[i]]$T2D_grt_cond, Yu_imp_result[[i]]$T2D_grt_cond),
               Zeller_imp_result[[i]]$T2D_grt_cond)
  
  CRC_grt_overlap_imputed[i,16] <- CRC_grt_overlap_imputed[i,5] / length(FV)
  CRC_grt_overlap_imputed[i,17] <- CRC_grt_overlap_imputed[i,6] / length(FY)
  CRC_grt_overlap_imputed[i,18] <- CRC_grt_overlap_imputed[i,7] / length(FZ)
  CRC_grt_overlap_imputed[i,19] <- CRC_grt_overlap_imputed[i,8] / length(VY)
  CRC_grt_overlap_imputed[i,20] <- CRC_grt_overlap_imputed[i,9] / length(VZ)
  CRC_grt_overlap_imputed[i,21] <- CRC_grt_overlap_imputed[i,10] / length(YZ)
  
  CRC_grt_overlap_imputed[i,22] <- CRC_grt_overlap_imputed[i,11] / length(FVY)
  CRC_grt_overlap_imputed[i,23] <- CRC_grt_overlap_imputed[i,12] / length(FVZ)
  CRC_grt_overlap_imputed[i,24] <- CRC_grt_overlap_imputed[i,13] / length(FYZ)
  CRC_grt_overlap_imputed[i,25] <- CRC_grt_overlap_imputed[i,14] / length(VYZ)
  CRC_grt_overlap_imputed[i,26] <- CRC_grt_overlap_imputed[i,15] / length(all_union)
}
rownames(CRC_grt_overlap_imputed) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2","Omnibus")
colnames(CRC_grt_overlap_imputed) <- c("Feng", "Vogtmann", "Yu", "Zeller", "FV", "FY", "FZ", "VY", "VZ", "YZ", "FVY", "FVZ", "FYZ", "VYZ", "FVYZ", "FV_rate", "FY_rate", "FZ_rate", "VY_rate", "VZ_rate", "YZ_rate", "FVY_rate", "FVZ_rate", "FYZ_rate", "VYZ_rate", "FVYZ_rate")





CRC_leq_overlap_raw <- matrix(nrow = 8, ncol = 26)
for(i in 1:8){
  CRC_leq_overlap_raw[i,1] <- length(Feng_raw_result[[i]]$T2D_leq_cond)
  CRC_leq_overlap_raw[i,2] <- length(Vogtmann_raw_result[[i]]$T2D_leq_cond)
  CRC_leq_overlap_raw[i,3] <- length(Yu_raw_result[[i]]$T2D_leq_cond)
  CRC_leq_overlap_raw[i,4] <- length(Zeller_raw_result[[i]]$T2D_leq_cond)
  
  all_intersect <- intersect(intersect(intersect(Vogtmann_raw_result[[i]]$T2D_leq_cond, Yu_raw_result[[i]]$T2D_leq_cond),
                                       Zeller_raw_result[[i]]$T2D_leq_cond),Feng_raw_result[[i]]$T2D_leq_cond)
  all_union <- union(union(union(Feng_raw_result[[i]]$T2D_leq_cond, Vogtmann_raw_result[[i]]$T2D_leq_cond), 
                           Yu_raw_result[[i]]$T2D_leq_cond),  Zeller_raw_result[[i]]$T2D_leq_cond)
  FV <- intersect(Feng_raw_result[[i]]$T2D_leq_cond, Vogtmann_raw_result[[i]]$T2D_leq_cond)
  FY <- intersect(Feng_raw_result[[i]]$T2D_leq_cond, Yu_raw_result[[i]]$T2D_leq_cond)
  FZ <- intersect(Feng_raw_result[[i]]$T2D_leq_cond, Zeller_raw_result[[i]]$T2D_leq_cond)
  VY <- intersect(Vogtmann_raw_result[[i]]$T2D_leq_cond, Yu_raw_result[[i]]$T2D_leq_cond)
  VZ <- intersect(Vogtmann_raw_result[[i]]$T2D_leq_cond, Zeller_raw_result[[i]]$T2D_leq_cond)
  YZ <- intersect(Yu_raw_result[[i]]$T2D_leq_cond, Zeller_raw_result[[i]]$T2D_leq_cond)
  FVY <- intersect(intersect(Feng_raw_result[[i]]$T2D_leq_cond, Vogtmann_raw_result[[i]]$T2D_leq_cond),
                   Yu_raw_result[[i]]$T2D_leq_cond)
  FVZ <- intersect(intersect(Feng_raw_result[[i]]$T2D_leq_cond, Vogtmann_raw_result[[i]]$T2D_leq_cond),
                   Zeller_raw_result[[i]]$T2D_leq_cond)
  FYZ <- intersect(intersect(Feng_raw_result[[i]]$T2D_leq_cond, Yu_raw_result[[i]]$T2D_leq_cond),
                   Zeller_raw_result[[i]]$T2D_leq_cond)
  VYZ <- intersect(intersect(Vogtmann_raw_result[[i]]$T2D_leq_cond, Yu_raw_result[[i]]$T2D_leq_cond),
                   Zeller_raw_result[[i]]$T2D_leq_cond)
  
  CRC_leq_overlap_raw[i,15] <- length(all_intersect)
  CRC_leq_overlap_raw[i,11] <- length(FVY[!(FVY %in% all_intersect)])
  CRC_leq_overlap_raw[i,12] <- length(FVZ[!(FVZ %in% all_intersect)])
  CRC_leq_overlap_raw[i,13] <- length(FYZ[!(FYZ %in% all_intersect)])
  CRC_leq_overlap_raw[i,14] <- length(VYZ[!(VYZ %in% all_intersect)])
  
  CRC_leq_overlap_raw[i,5] <- length(FV[!(FV %in% union(FVY, FVZ))])
  CRC_leq_overlap_raw[i,6] <- length(FY[!(FY %in% union(FVY, FYZ))])
  CRC_leq_overlap_raw[i,7] <- length(FZ[!(FZ %in% union(FVZ, FYZ))])
  CRC_leq_overlap_raw[i,8] <- length(VY[!(VY %in% union(FVY, VYZ))])
  CRC_leq_overlap_raw[i,9] <- length(VZ[!(VZ %in% union(FVZ, VYZ))])
  CRC_leq_overlap_raw[i,10] <- length(YZ[!(YZ %in% union(FYZ, VYZ))])
  
  FV <- union(Feng_raw_result[[i]]$T2D_leq_cond, Vogtmann_raw_result[[i]]$T2D_leq_cond)
  FY <- union(Feng_raw_result[[i]]$T2D_leq_cond, Yu_raw_result[[i]]$T2D_leq_cond)
  FZ <- union(Feng_raw_result[[i]]$T2D_leq_cond, Zeller_raw_result[[i]]$T2D_leq_cond)
  VY <- union(Vogtmann_raw_result[[i]]$T2D_leq_cond, Yu_raw_result[[i]]$T2D_leq_cond)
  VZ <- union(Vogtmann_raw_result[[i]]$T2D_leq_cond, Zeller_raw_result[[i]]$T2D_leq_cond)
  YZ <- union(Yu_raw_result[[i]]$T2D_leq_cond, Zeller_raw_result[[i]]$T2D_leq_cond)
  FVY <- union(union(Feng_raw_result[[i]]$T2D_leq_cond, Vogtmann_raw_result[[i]]$T2D_leq_cond),
               Yu_raw_result[[i]]$T2D_leq_cond)
  FVZ <- union(union(Feng_raw_result[[i]]$T2D_leq_cond, Vogtmann_raw_result[[i]]$T2D_leq_cond),
               Zeller_raw_result[[i]]$T2D_leq_cond)
  FYZ <- union(union(Feng_raw_result[[i]]$T2D_leq_cond, Yu_raw_result[[i]]$T2D_leq_cond),
               Zeller_raw_result[[i]]$T2D_leq_cond)
  VYZ <- union(union(Vogtmann_raw_result[[i]]$T2D_leq_cond, Yu_raw_result[[i]]$T2D_leq_cond),
               Zeller_raw_result[[i]]$T2D_leq_cond)
  
  CRC_leq_overlap_raw[i,16] <- CRC_leq_overlap_raw[i,5] / length(FV)
  CRC_leq_overlap_raw[i,17] <- CRC_leq_overlap_raw[i,6] / length(FY)
  CRC_leq_overlap_raw[i,18] <- CRC_leq_overlap_raw[i,7] / length(FZ)
  CRC_leq_overlap_raw[i,19] <- CRC_leq_overlap_raw[i,8] / length(VY)
  CRC_leq_overlap_raw[i,20] <- CRC_leq_overlap_raw[i,9] / length(VZ)
  CRC_leq_overlap_raw[i,21] <- CRC_leq_overlap_raw[i,10] / length(YZ)
  
  CRC_leq_overlap_raw[i,22] <- CRC_leq_overlap_raw[i,11] / length(FVY)
  CRC_leq_overlap_raw[i,23] <- CRC_leq_overlap_raw[i,12] / length(FVZ)
  CRC_leq_overlap_raw[i,24] <- CRC_leq_overlap_raw[i,13] / length(FYZ)
  CRC_leq_overlap_raw[i,25] <- CRC_leq_overlap_raw[i,14] / length(VYZ)
  CRC_leq_overlap_raw[i,26] <- CRC_leq_overlap_raw[i,15] / length(all_union)
}
rownames(CRC_leq_overlap_raw) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2","Omnibus")
colnames(CRC_leq_overlap_raw) <- c("Feng", "Vogtmann", "Yu", "Zeller", "FV", "FY", "FZ", "VY", "VZ", "YZ", "FVY", "FVZ", "FYZ", "VYZ", "FVYZ", "FV_rate", "FY_rate", "FZ_rate", "VY_rate", "VZ_rate", "YZ_rate", "FVY_rate", "FVZ_rate", "FYZ_rate", "VYZ_rate", "FVYZ_rate")

CRC_leq_overlap_imputed <- matrix(nrow = 8, ncol = 26)
for(i in 1:8){
  CRC_leq_overlap_imputed[i,1] <- length(Feng_imp_result[[i]]$T2D_leq_cond)
  CRC_leq_overlap_imputed[i,2] <- length(Vogtmann_imp_result[[i]]$T2D_leq_cond)
  CRC_leq_overlap_imputed[i,3] <- length(Yu_imp_result[[i]]$T2D_leq_cond)
  CRC_leq_overlap_imputed[i,4] <- length(Zeller_imp_result[[i]]$T2D_leq_cond)
  
  all_intersect <- intersect(intersect(intersect(Vogtmann_imp_result[[i]]$T2D_leq_cond, Yu_imp_result[[i]]$T2D_leq_cond),
                                       Zeller_imp_result[[i]]$T2D_leq_cond),Feng_imp_result[[i]]$T2D_leq_cond)
  all_union <- union(union(union(Feng_imp_result[[i]]$T2D_leq_cond, Vogtmann_imp_result[[i]]$T2D_leq_cond), 
                           Yu_imp_result[[i]]$T2D_leq_cond),  Zeller_imp_result[[i]]$T2D_leq_cond)
  FV <- intersect(Feng_imp_result[[i]]$T2D_leq_cond, Vogtmann_imp_result[[i]]$T2D_leq_cond)
  FY <- intersect(Feng_imp_result[[i]]$T2D_leq_cond, Yu_imp_result[[i]]$T2D_leq_cond)
  FZ <- intersect(Feng_imp_result[[i]]$T2D_leq_cond, Zeller_imp_result[[i]]$T2D_leq_cond)
  VY <- intersect(Vogtmann_imp_result[[i]]$T2D_leq_cond, Yu_imp_result[[i]]$T2D_leq_cond)
  VZ <- intersect(Vogtmann_imp_result[[i]]$T2D_leq_cond, Zeller_imp_result[[i]]$T2D_leq_cond)
  YZ <- intersect(Yu_imp_result[[i]]$T2D_leq_cond, Zeller_imp_result[[i]]$T2D_leq_cond)
  FVY <- intersect(intersect(Feng_imp_result[[i]]$T2D_leq_cond, Vogtmann_imp_result[[i]]$T2D_leq_cond),
                   Yu_imp_result[[i]]$T2D_leq_cond)
  FVZ <- intersect(intersect(Feng_imp_result[[i]]$T2D_leq_cond, Vogtmann_imp_result[[i]]$T2D_leq_cond),
                   Zeller_imp_result[[i]]$T2D_leq_cond)
  FYZ <- intersect(intersect(Feng_imp_result[[i]]$T2D_leq_cond, Yu_imp_result[[i]]$T2D_leq_cond),
                   Zeller_imp_result[[i]]$T2D_leq_cond)
  VYZ <- intersect(intersect(Vogtmann_imp_result[[i]]$T2D_leq_cond, Yu_imp_result[[i]]$T2D_leq_cond),
                   Zeller_imp_result[[i]]$T2D_leq_cond)
  
  CRC_leq_overlap_imputed[i,15] <- length(all_intersect)
  CRC_leq_overlap_imputed[i,11] <- length(FVY[!(FVY %in% all_intersect)])
  CRC_leq_overlap_imputed[i,12] <- length(FVZ[!(FVZ %in% all_intersect)])
  CRC_leq_overlap_imputed[i,13] <- length(FYZ[!(FYZ %in% all_intersect)])
  CRC_leq_overlap_imputed[i,14] <- length(VYZ[!(VYZ %in% all_intersect)])
  
  CRC_leq_overlap_imputed[i,5] <- length(FV[!(FV %in% union(FVY, FVZ))])
  CRC_leq_overlap_imputed[i,6] <- length(FY[!(FY %in% union(FVY, FYZ))])
  CRC_leq_overlap_imputed[i,7] <- length(FZ[!(FZ %in% union(FVZ, FYZ))])
  CRC_leq_overlap_imputed[i,8] <- length(VY[!(VY %in% union(FVY, VYZ))])
  CRC_leq_overlap_imputed[i,9] <- length(VZ[!(VZ %in% union(FVZ, VYZ))])
  CRC_leq_overlap_imputed[i,10] <- length(YZ[!(YZ %in% union(FYZ, VYZ))])
  
  FV <- union(Feng_imp_result[[i]]$T2D_leq_cond, Vogtmann_imp_result[[i]]$T2D_leq_cond)
  FY <- union(Feng_imp_result[[i]]$T2D_leq_cond, Yu_imp_result[[i]]$T2D_leq_cond)
  FZ <- union(Feng_imp_result[[i]]$T2D_leq_cond, Zeller_imp_result[[i]]$T2D_leq_cond)
  VY <- union(Vogtmann_imp_result[[i]]$T2D_leq_cond, Yu_imp_result[[i]]$T2D_leq_cond)
  VZ <- union(Vogtmann_imp_result[[i]]$T2D_leq_cond, Zeller_imp_result[[i]]$T2D_leq_cond)
  YZ <- union(Yu_imp_result[[i]]$T2D_leq_cond, Zeller_imp_result[[i]]$T2D_leq_cond)
  FVY <- union(union(Feng_imp_result[[i]]$T2D_leq_cond, Vogtmann_imp_result[[i]]$T2D_leq_cond),
               Yu_imp_result[[i]]$T2D_leq_cond)
  FVZ <- union(union(Feng_imp_result[[i]]$T2D_leq_cond, Vogtmann_imp_result[[i]]$T2D_leq_cond),
               Zeller_imp_result[[i]]$T2D_leq_cond)
  FYZ <- union(union(Feng_imp_result[[i]]$T2D_leq_cond, Yu_imp_result[[i]]$T2D_leq_cond),
               Zeller_imp_result[[i]]$T2D_leq_cond)
  VYZ <- union(union(Vogtmann_imp_result[[i]]$T2D_leq_cond, Yu_imp_result[[i]]$T2D_leq_cond),
               Zeller_imp_result[[i]]$T2D_leq_cond)
  
  CRC_leq_overlap_imputed[i,16] <- CRC_leq_overlap_imputed[i,5] / length(FV)
  CRC_leq_overlap_imputed[i,17] <- CRC_leq_overlap_imputed[i,6] / length(FY)
  CRC_leq_overlap_imputed[i,18] <- CRC_leq_overlap_imputed[i,7] / length(FZ)
  CRC_leq_overlap_imputed[i,19] <- CRC_leq_overlap_imputed[i,8] / length(VY)
  CRC_leq_overlap_imputed[i,20] <- CRC_leq_overlap_imputed[i,9] / length(VZ)
  CRC_leq_overlap_imputed[i,21] <- CRC_leq_overlap_imputed[i,10] / length(YZ)
  
  CRC_leq_overlap_imputed[i,22] <- CRC_leq_overlap_imputed[i,11] / length(FVY)
  CRC_leq_overlap_imputed[i,23] <- CRC_leq_overlap_imputed[i,12] / length(FVZ)
  CRC_leq_overlap_imputed[i,24] <- CRC_leq_overlap_imputed[i,13] / length(FYZ)
  CRC_leq_overlap_imputed[i,25] <- CRC_leq_overlap_imputed[i,14] / length(VYZ)
  CRC_leq_overlap_imputed[i,26] <- CRC_leq_overlap_imputed[i,15] / length(all_union)
}
rownames(CRC_leq_overlap_imputed) <- c("Wilcox", "t_test", "ANCOM", "ZINB", "NB", "Metaseq", "DEseq2","Omnibus")
colnames(CRC_leq_overlap_imputed) <- c("Feng", "Vogtmann", "Yu", "Zeller", "FV", "FY", "FZ", "VY", "VZ", "YZ", "FVY", "FVZ", "FYZ", "VYZ", "FVYZ", "FV_rate", "FY_rate", "FZ_rate", "VY_rate", "VZ_rate", "YZ_rate", "FVY_rate", "FVZ_rate", "FYZ_rate", "VYZ_rate", "FVYZ_rate")




Qin_Karlsson_result <- classification(input_data = Qin_raw, condition = Qin_condition, features_raw = Karlsson_raw_result, features_imputed = Karlsson_imp_result)

Karlsson_Qin_result <- classification(input_data = Karlsson_raw, condition = Karlsson_condition, features_raw = Qin_raw_result, features_imputed = Qin_imp_result)

Feng_Vogtmann_result <- classification(input_data = Feng_raw, condition = Feng_condition, features_raw = Vogtmann_raw_result, features_imputed = Vogtmann_imp_result)

Vogtmann_Feng_result <- classification(input_data = Vogtmann_raw, condition = Vogtmann_condition, features_raw = Feng_raw_result, features_imputed = Feng_imp_result)

Yu_Zeller_result <- classification(input_data = Yu_raw, condition = Yu_condition, features_raw = Zeller_raw_result, features_imputed = Zeller_imp_result)

Zeller_Yu_result <- classification(input_data = Zeller_raw, condition = Zeller_condition, features_raw = Yu_raw_result, features_imputed = Yu_imp_result)




classification_plot_pr <- function(result){
  hz_value <-  result$full_data_result[11]
  df <- rbind(cbind(c(result$raw_result[,11], result$imputed_result[,11]), c(paste0("raw_",  rownames(result$raw_result)), paste0("imp_",  rownames(result$raw_result))), rep( rownames(result$raw_result), 2)))
  df <- df[-c(5, 12),]
  df[4,3] <- "ZINB"
  df[11,3] <- "NB"
  colnames(df) <- c("values", "type", "test")
  df <- as.data.frame(df)
  df$values <- as.numeric(as.character(df$values))
  df$values[is.na(df$values)] <- 0.5
  df$type <- factor(as.character(df$type), levels = c(paste0("raw_",  rownames(result$raw_result)[-5]), paste0("imp_",  rownames(result$raw_result)[-4])))
  new_levels <- rownames(result$raw_result)[-5]
  new_levels[4] <- "ZINB/NB"
  
  df$test <- factor(as.character(df$test), levels = new_levels)
  df$type <- ifelse(grepl("raw_", df$type), "raw_result", "imputed_result")
  df$type <- factor(df$type, levels = c("raw_result", "imputed_result"))
  ggplot(df, aes(fill=type, y=values, x=test)) + 
    geom_bar(position = "dodge", stat="identity") + 
    scale_x_discrete(labels=c("Wilcox", "t_test", "ANCOM", 
                              "ZINB/NB", "Metaseq", "DEseq2", "Omnibus")) +
    scale_fill_manual(
      name="",
      values=c("raw_result"=alpha("#636363", 0.5), "imputed_result"="#636363"),
      labels=c("raw_result"="raw_result", "imputed_result"="imputed_result")
    ) + 
    #scale_fill_manual("method", values = c("ZI" = "#7e7e7e", "mbImpute" = "#1D9E78", "softImpute" = "#346FA0", "scImpute" = "#3383C2", "SAVER" = "#53AFDA", "MAGIC" = "#AAD0E6", "ALRA" = "#deebf7")) +
    ylab("PR-AUC") +
    coord_cartesian(ylim=c(0.2, 1)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(size=12, angle=45, hjust=1),
                       #axis.text.x=element_blank(),
                       axis.title.y=element_text(size=14),
                       legend.text=element_text(size=12),
                       panel.border = element_blank(),
                       axis.line = element_line(colour = "black"),
                       legend.position="right",
                       legend.title = element_blank()) +geom_hline(yintercept = hz_value, linetype="dashed", color = "#7e7e7e", size = 2)
}


pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Qin_Karlsson_result_random_forest.pdf")
classification_plot_pr(result = Qin_Karlsson_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Karlsson_Qin_result_random_forest.pdf")
classification_plot_pr(result = Karlsson_Qin_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Qin_Karlsson_result_xgboost.pdf")
classification_plot_pr(result = Qin_Karlsson_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Karlsson_Qin_result_xgboost.pdf")
classification_plot_pr(result = Karlsson_Qin_result)
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Qin_Karlsson_result_svm_lr.pdf")
classification_plot_pr(result = Qin_Karlsson_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Karlsson_Qin_result_svm_lr.pdf")
classification_plot_pr(result = Karlsson_Qin_result)
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Qin_Karlsson_result_svm_guass.pdf")
classification_plot_pr(result = Qin_Karlsson_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Karlsson_Qin_result_svm_guass.pdf")
classification_plot_pr(result = Karlsson_Qin_result)
dev.off()



pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Feng_Vogtmann_result_xgboost.pdf")
classification_plot_pr(result = Feng_Vogtmann_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Vogtmann_Feng_result_xgboost.pdf")
classification_plot_pr(result = Vogtmann_Feng_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Feng_Vogtmann_result_random_forest.pdf")
classification_plot_pr(result = Feng_Vogtmann_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Vogtmann_Feng_result_random_forest.pdf")
classification_plot_pr(result = Vogtmann_Feng_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Yu_Zeller_result_random_forest.pdf")
classification_plot_pr(result = Yu_Zeller_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Zeller_Yu_result_random_forest.pdf")
classification_plot_pr(result = Zeller_Yu_result)
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Yu_Zeller_result_xgboost.pdf")
classification_plot_pr(result = Yu_Zeller_result)
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/change/Zeller_Yu_result_xgboost.pdf")
classification_plot_pr(result = Zeller_Yu_result)
dev.off()



DE_comparison_Karlsson <- function(taxon_name, xlim = c(-0.3,6), ylim = c(0, 1.4)){
  par(mfrow = c(2,1))
  # hist(Qin_raw[,taxon_name])
  
  # hist(Karlsson_raw[,taxon_name])
  hist(Karlsson_raw[Karlsson_condition == "T2D",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("red", 0.3), main = taxon_name, xlab ="" )
  hist(Karlsson_raw[Karlsson_condition == "control",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("#636363", 0.3), xlab = "", add=T)
  
  # hist(Karlsson_imp[,taxon_name])
  hist(Karlsson_imp[Karlsson_condition == "T2D",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("red", 0.3), main = "", xlab = "")
  hist(Karlsson_imp[Karlsson_condition == "control",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("#636363", 0.3), add=T, xlab = "")
}
# p value for DESeq2 before and after imputation:
DE_comparison_Karlsson("s__Ruminococcus_sp_5_1_39BFAA")
DE_comparison_Karlsson("s__Ruminococcus_callidus")
DE_comparison_Karlsson("s__Ruminococcus_albus")
DE_comparison_Karlsson("s__Clostridium_hathewayi")
DE_comparison_Karlsson("s__Blautia_hydrogenotrophica")
DE_comparison_Karlsson("s__Bacteroides_finegoldii")
DE_comparison_Karlsson("s__Prevotella_copri")
DE_comparison_Karlsson("s__Enterococcus_faecalis")
DE_comparison_Karlsson("s__Clostridium_bolteae")

par(mfrow = c(1,1))
taxon_name = "s__Ruminococcus_sp_5_1_39BFAA"
xlim = c(-0.3,6)
ylim = c(0, 1)
pdf(paste0("Real data distribution/", taxon_name, "_raw.pdf"), width = 7, height = 5)
hist(Karlsson_raw[Karlsson_condition == "T2D",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("red", 0.3), main = "", xlab ="" )
hist(Karlsson_raw[Karlsson_condition == "control",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("#636363", 0.3), xlab = "", add=T)
dev.off()
pdf(paste0("Real data distribution/", taxon_name, "_imp.pdf"), width = 7, height = 5)  
# hist(Karlsson_imp[,taxon_name])
hist(Karlsson_imp[Karlsson_condition == "T2D",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("red", 0.3), main = "", xlab = "")
hist(Karlsson_imp[Karlsson_condition == "control",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("#636363", 0.3), add=T, xlab = "")
dev.off()


par(mfrow = c(1,1))
taxon_name = "s__Ruminococcus_callidus"
xlim = c(-0.3,6)
ylim = c(0, 1)
pdf(paste0("Real data distribution/", taxon_name, "_raw.pdf"), width = 7, height = 5)
hist(Karlsson_raw[Karlsson_condition == "T2D",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("red", 0.3), main = "", xlab ="" )
hist(Karlsson_raw[Karlsson_condition == "control",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("#636363", 0.3), xlab = "", add=T)
dev.off()
pdf(paste0("Real data distribution/", taxon_name, "_imp.pdf"), width = 7, height = 5)  
# hist(Karlsson_imp[,taxon_name])
hist(Karlsson_imp[Karlsson_condition == "T2D",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("red", 0.3), main = "", xlab = "")
hist(Karlsson_imp[Karlsson_condition == "control",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("#636363", 0.3), add=T, xlab = "")
dev.off()

par(mfrow = c(1,1))
taxon_name = "s__Ruminococcus_albus"
xlim = c(-0.3,6)
ylim = c(0, 1)
pdf(paste0("Real data distribution/", taxon_name, "_raw.pdf"), width = 7, height = 5)
hist(Karlsson_raw[Karlsson_condition == "T2D",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("red", 0.3), main = "", xlab ="" )
hist(Karlsson_raw[Karlsson_condition == "control",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("#636363", 0.3), xlab = "", add=T)
dev.off()
pdf(paste0("Real data distribution/", taxon_name, "_imp.pdf"), width = 7, height = 5)  
# hist(Karlsson_imp[,taxon_name])
hist(Karlsson_imp[Karlsson_condition == "T2D",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("red", 0.3), main = "", xlab = "")
hist(Karlsson_imp[Karlsson_condition == "control",taxon_name], freq = FALSE, breaks = seq(0, 6, 0.5), xlim = xlim, ylim = ylim, col = alpha("#636363", 0.3), add=T, xlab = "")
dev.off()






par(mfrow = c(1,1))
taxon_name = "s__Ruminococcus_sp_5_1_39BFAA"
xlim = c(-0.3,6)
pdf(paste0("Real data distribution/", taxon_name, "_raw.pdf"), width = 7, height = 5)
hist(Karlsson_raw[Karlsson_condition == "T2D",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim, col = alpha("red", 0.3), main = "", xlab ="" )
hist(Karlsson_raw[Karlsson_condition == "control",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim, col = alpha("#636363", 0.3), xlab = "", add=T)
dev.off()
pdf(paste0("Real data distribution/", taxon_name, "_imp.pdf"), width = 7, height = 5)  
# hist(Karlsson_imp[,taxon_name])
hist(Karlsson_imp[Karlsson_condition == "T2D",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim,  col = alpha("red", 0.3), main = "", xlab = "")
hist(Karlsson_imp[Karlsson_condition == "control",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim,  col = alpha("#636363", 0.3), add=T, xlab = "")
dev.off()

par(mfrow = c(1,1))
taxon_name = "s__Ruminococcus_callidus"
xlim = c(-0.3,6)
ylim = c(0, 1)
pdf(paste0("Real data distribution/", taxon_name, "_raw.pdf"), width = 7, height = 5)
hist(Karlsson_raw[Karlsson_condition == "T2D",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim,  col = alpha("red", 0.3), main = "", xlab ="" )
hist(Karlsson_raw[Karlsson_condition == "control",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim,  col = alpha("#636363", 0.3), xlab = "", add=T)
dev.off()
pdf(paste0("Real data distribution/", taxon_name, "_imp.pdf"), width = 7, height = 5)  
# hist(Karlsson_imp[,taxon_name])
hist(Karlsson_imp[Karlsson_condition == "T2D",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim,  col = alpha("red", 0.3), main = "", xlab = "")
hist(Karlsson_imp[Karlsson_condition == "control",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim,  col = alpha("#636363", 0.3), add=T, xlab = "")
dev.off()

par(mfrow = c(1,1))
taxon_name = "s__Ruminococcus_albus"
xlim = c(-0.3,6)
ylim = c(0, 1)
pdf(paste0("Real data distribution/", taxon_name, "_raw.pdf"), width = 7, height = 5)
hist(Karlsson_raw[Karlsson_condition == "T2D",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim,  col = alpha("red", 0.3), main = "", xlab ="" )
hist(Karlsson_raw[Karlsson_condition == "control",taxon_name], freq = TRUE, breaks = seq(0, 6, 0.5), xlim = xlim,  col = alpha("#636363", 0.3), xlab = "", add=T)
dev.off()



DE_comparison_Yu <- function(taxon_name){
  par(mfrow = c(2,1), mar = c(5, 8, 4, 2) + 0.1) 
  taxon_title <- unlist(strsplit(taxon_name, "[.]"))  
  taxon_title <- taxon_title[length(taxon_title)]     
  taxon_title <- sub(".*s__", "", taxon_title)       

  hist(Yu_raw[Yu_condition == "CRC",taxon_name], freq = TRUE, breaks = "Sturges", col = alpha("red", 0.3), main = taxon_title, xlab = "Taxon abundances in log-scale", ylab = "Frequency", xlim=c(0, 6) )
  hist(Yu_raw[Yu_condition == "control",taxon_name], freq = TRUE, breaks = "Sturges", col = alpha("#636363", 0.3), add=T, xlim=c(0, 6))
  y_range1 <- par("usr")[3:4]
  y_mid1 <- mean(y_range1)
  mtext("Raw\nYu et al.", side = 2, line = 6, at = y_mid1, las = 1,adj = 0.5) 
  p_value1 <- Yu_raw_result$DEseq2_pval[which(colnames(Yu_raw) == taxon_name)]
  x_pos <- usr[2] * 0.98 
  y_pos <- usr[4] * 0.93 
  text(x_pos, y_pos, labels = paste("p.adj = ", format(p_value1, digits = 2)), 
       pos = 2, cex = 0.8, col = "black")
  hist(Yu_imp[Yu_condition == "CRC",taxon_name], freq = TRUE, breaks = "Sturges", col = alpha("red", 0.3), main = taxon_title,xlab = "Taxon abundances in log-scale", ylab = "Frequency", xlim=c(0, 6))
  hist(Yu_imp[Yu_condition == "control",taxon_name], freq = TRUE, breaks = "Sturges", col = alpha("#636363", 0.3), add=T, xlim=c(0, 6))
  y_range2 <- par("usr")[3:4]
  y_mid2 <- mean(y_range2)
  mtext("Imputed\nYu et al.", side = 2, line = 6, at = y_mid2, las = 1,adj = 0.5) 
  p_value2 <- Yu_imp_result$DEseq2_pval[which(colnames(Yu_imp) == taxon_name)]
  usr <- par("usr")
  
  x_pos <- usr[2] * 0.98 
  y_pos <- usr[4] * 0.93 
  text(x_pos, y_pos, labels = paste("p.adj = ", format(p_value2, digits = 2)), 
       pos = 2, cex = 0.8, col = "black")
  usr <- par("usr")
  x_pos <- usr[1] + diff(usr[1:2]) * 0.05 
  y_pos <- usr[4] - diff(usr[3:4]) * 0.1  
  legend(x_pos, y_pos, legend=c("CRC", "control"), fill=c(alpha("red", 0.3), alpha("#636363", 0.3)), cex=0.8, bty="n", xjust=0, yjust=1)
}                                           
pdf("/Users/hanxinyu/Desktop/real_data_analysis/DE_comparison_Yu_plots.pdf", width = 8, height = 6)

for(x in Yu_imp_result$DEseq2$T2D_grt_cond[which(!Yu_imp_result$DEseq2$T2D_grt_cond %in% Yu_intersect)]){
  DE_comparison_Yu(x)
}
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/DE_comparison_Yu_plots1.pdf", width = 8, height = 6)

for(x in Yu_raw_result$DEseq2$T2D_grt_cond[which(!Yu_raw_result$DEseq2$T2D_grt_cond %in% Yu_intersect)]){
  DE_comparison_Yu(x)
}
dev.off()


DE_comparison_Zeller <- function(taxon_name){
  par(mfrow = c(2,1), mar = c(5, 8, 3, 2) + 0.1) 
  taxon_title <- unlist(strsplit(taxon_name, "[.]"))  
  taxon_title <- taxon_title[length(taxon_title)]     
  taxon_title <- sub(".*s__", "", taxon_title)       
  
  
  hist(Zeller_raw[Zeller_condition == "CRC",taxon_name], freq = TRUE, breaks = "Sturges", col = alpha("red", 0.3), main = taxon_title, xlab = "Taxon abundances in log-scale", ylab = "Frequency", xlim=c(0, 6) )
  hist(Zeller_raw[Zeller_condition == "control",taxon_name], freq = TRUE, breaks = "Sturges", col = alpha("#636363", 0.3), add=T, xlim=c(0, 6))
  y_range1 <- par("usr")[3:4]
  y_mid1 <- mean(y_range1)
  
  mtext("Raw\nZeller et al.", side = 2, line = 6, at = y_mid1, las = 1,adj = 0.5) 
  p_value1 <- Zeller_raw_result$DEseq2_pval[which(colnames(Zeller_raw) == taxon_name)]
  usr <- par("usr")
  
  x_pos <- usr[2] * 0.98 
  y_pos <- usr[4] * 0.93 
  text(x_pos, y_pos, labels = paste("p.adj = ", format(p_value1, digits = 2)), 
       pos = 2, cex = 0.8, col = "black")
  hist(Zeller_imp[Zeller_condition == "CRC",taxon_name], freq = TRUE, breaks = "Sturges", col = alpha("red", 0.3), main = taxon_title,xlab = "Taxon abundances in log-scale", ylab = "Frequency", xlim=c(0, 6))
  hist(Zeller_imp[Zeller_condition == "control",taxon_name], freq = TRUE, breaks = "Sturges", col = alpha("#636363", 0.3), add=T, xlim=c(0, 6))
  y_range2 <- par("usr")[3:4]
  y_mid2 <- mean(y_range2)
  mtext("Imputed\nZeller et al.", side = 2, line = 6, at = y_mid2, las = 1,adj = 0.5) 
  p_value2 <- Zeller_imp_result$DEseq2_pval[which(colnames(Zeller_imp) == taxon_name)]
  usr <- par("usr")
  x_pos <- usr[2] * 0.98 
  y_pos <- usr[4] * 0.93 
  
  text(x_pos, y_pos, labels = paste("p.adj = ", format(p_value2, digits = 2)), 
       pos = 2, cex = 0.8, col = "black")
  usr <- par("usr")
  
  x_pos <- usr[1] + diff(usr[1:2]) * 0.05 
  y_pos <- usr[4] - diff(usr[3:4]) * 0.1 
  legend(x_pos, y_pos, legend=c("CRC", "control"), fill=c(alpha("red", 0.3), alpha("#636363", 0.3)), cex=0.8, bty="n", xjust=0, yjust=1)
}            

Zeller_raw_result$DEseq2$T2D_grt_cond
Zeller_imp_result$DEseq2$T2D_grt_cond
length(Zeller_raw_result$DEseq2$T2D_grt_cond)
length(Zeller_imp_result$DEseq2$T2D_grt_cond)

Zeller_intersect <- intersect(Zeller_raw_result$DEseq2$T2D_grt_cond,
                              Zeller_imp_result$DEseq2$T2D_grt_cond)
length(Zeller_intersect)
Zeller_raw_result$DEseq2$T2D_grt_cond[which(!Zeller_raw_result$DEseq2$T2D_grt_cond %in% Zeller_intersect)]
Zeller_imp_result$DEseq2$T2D_grt_cond[which(!Zeller_imp_result$DEseq2$T2D_grt_cond %in% Zeller_intersect)]


pdf("/Users/hanxinyu/Desktop/real_data_analysis/DE_comparison_Zeller_plots.pdf", width = 8, height = 7)

for(x in Zeller_imp_result$DEseq2$T2D_grt_cond[which(!Zeller_imp_result$DEseq2$T2D_grt_cond %in% Zeller_intersect)]){
  DE_comparison_Zeller(x)
}
dev.off()



pdf("/Users/hanxinyu/Desktop/real_data_analysis/DE_comparison_Zeller_plots1.pdf", width = 8, height = 7)

for(x in Zeller_raw_result$DEseq2$T2D_grt_cond[which(!Zeller_raw_result$DEseq2$T2D_grt_cond %in% Zeller_intersect)]){
  DE_comparison_Zeller(x)
}
dev.off()




tolerance <- .Machine$double.eps^0.5
Qin_raw[abs(Qin_raw - log10(1.01)) < tolerance] <- log10(1.01)

y_sim = Qin_raw
y_imp <- Qin_imp

Raw_cor <- cor(y_sim)
imp_cor <- cor(y_imp)

cor(as.vector(Raw_cor), as.vector(imp_cor))
cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( imp_cor[upper.tri(imp_cor)] ))


dim(Raw_cor)[1]

nz_cor <- matrix(nrow = 179, ncol = 179)
for(i in 2:dim(y_imp)[2]){
  for(j in 1:(i-1)){
    nz_idx <- intersect(which(y_sim[,i] != log10(1.01)), which(y_sim[,j] != log10(1.01)))
    nz_cor[j,i] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
    nz_cor[i,j] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
  }
}
nz_cor[is.na(nz_cor)] <- 0

diag(nz_cor) = 1
colnames(nz_cor) = colnames(Raw_cor)
rownames(nz_cor) = colnames(Raw_cor)

Qin_raw_cor = Raw_cor
Qin_imp_cor = imp_cor
Qin_nz_cor = nz_cor

cor_analysis <- function(y_sim, y_imp){
  Raw_cor <- cor(y_sim)
  imp_cor <- cor(y_imp)
  
  cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( imp_cor[upper.tri(imp_cor)] ))
  nz_cor <- matrix(nrow = dim(Raw_cor)[1], ncol = dim(Raw_cor)[1])
  for(i in 2:dim(y_imp)[2]){
    for(j in 1:(i-1)){
      nz_idx <- intersect(which(y_sim[,i] != log10(1.01)), which(y_sim[,j] != log10(1.01)))
      nz_cor[j,i] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
    }
  }
  nz_cor[is.na(nz_cor)] <- 0
  
  cor_result <- matrix(nrow = 2, ncol = 4)
  cor_result[1,1] = cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( nz_cor[upper.tri(nz_cor)] ))
  cor_result[2,1] = cor(as.vector( imp_cor[upper.tri(imp_cor)] ), as.vector( nz_cor[upper.tri(nz_cor)] ))
  cor_result[2,2] = cor(as.vector( Raw_cor[upper.tri(imp_cor)] ), as.vector( imp_cor[upper.tri(nz_cor)] ))
  
  
  cor_result[1,3] = cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( nz_cor[upper.tri(nz_cor)] ), method = "spearman")
  cor_result[2,3] = cor(as.vector( imp_cor[upper.tri(imp_cor)] ), as.vector( nz_cor[upper.tri(nz_cor)] ), method = "spearman")
  cor_result[2,4] = cor(as.vector( Raw_cor[upper.tri(imp_cor)] ), as.vector( imp_cor[upper.tri(nz_cor)] ), method = "spearman")
  return(cor_result)
}

df_gen_v2 <- function(df1, df2, df3){
  data_tp <- c(rep("all" , 2) , rep("T2D", 2) , rep("control", 2))
  comp <- rep( c("raw_nz", "imp_nz") , 3)
  value <- c(as.vector(df1[,1]), as.vector(df2[,1]), as.vector(df3[,1]))
  data <- data.frame(data_tp, comp, value)
  data$data_tp <- factor(data$data_tp, levels = c("all", "T2D", "control"))
  data$comp <- factor(data$comp, levels = c("raw_nz", "imp_nz"))
  return(data)
} 

df_gen_v3 <- function(df1, df2, df3){
  data_tp <- c(rep("all" , 2) , rep("T2D", 2) , rep("control", 2))
  comp <- rep( c("raw_nz", "imp_nz"), 3)
  value <- c(as.vector(df1[,3]), as.vector(df2[,3]), as.vector(df3[,3]))
  data <- data.frame(data_tp, comp, value)
  data$data_tp <- factor(data$data_tp, levels = c("all", "T2D", "control"))
  data$comp <- factor(data$comp, levels = c("raw_nz", "imp_nz"))
  return(data)
} 




Qin_condition[Qin_condition == "0"]="control"
Qin_condition[Qin_condition == "1"]="T2D"
df1 <- cor_analysis(y_sim[Qin_condition == "T2D",], y_imp[Qin_condition == "T2D",])

df2 <- cor_analysis(y_sim[Qin_condition == "control",], y_imp[Qin_condition == "control",])

df3 <- cor_analysis(y_sim, y_imp)

df <- df_gen_v2(df3, df1, df2)
df_spearman <- df_gen_v3(df3, df1, df2)




pdf("/Users/hanxinyu/Desktop/real_data_analysis/Qin/Qin_correlation.pdf", width = 10, height = 5)
ggplot(df, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) +  
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "T2D", "control"),
    labels = c("all", "T2D", "control")
  ) +
  ylab("correlation") +
  ggtitle("Qin et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),     
    axis.text.x = element_text(size = 20),   
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5),  
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Qin/Qin_correlation_spearman.pdf", width = 10, height = 5)
ggplot(df_spearman, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) +  
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "T2D", "control"),
    labels = c("all", "T2D", "control")
  ) +
  ylab("correlation") +
  ggtitle("Qin et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),       
    axis.text.x = element_text(size = 20),     
    axis.title = element_text(size = 20),      
    plot.title = element_text(size = 22, hjust = 0.5),  
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()


tolerance <- .Machine$double.eps^0.5
Karlsson_raw[abs(Karlsson_raw - log10(1.01)) < tolerance] <- log10(1.01)

y_sim = Karlsson_raw
y_imp <- Karlsson_imp

Raw_cor <- cor(y_sim)
imp_cor <- cor(y_imp)

cor(as.vector(Raw_cor), as.vector(imp_cor))
cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( imp_cor[upper.tri(imp_cor)] ))

dim(Raw_cor)[1]

nz_cor <- matrix(nrow = 181, ncol = 181)
for(i in 2:dim(y_imp)[2]){
  for(j in 1:(i-1)){
    nz_idx <- intersect(which(y_sim[,i] != log10(1.01)), which(y_sim[,j] != log10(1.01)))
    nz_cor[j,i] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
    nz_cor[i,j] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
  }
}
nz_cor[is.na(nz_cor)] <- 0

diag(nz_cor) = 1
colnames(nz_cor) = colnames(Raw_cor)
rownames(nz_cor) = colnames(Raw_cor)

Karlsson_raw_cor = Raw_cor
Karlsson_imp_cor = imp_cor
Karlsson_nz_cor = nz_cor

Karlsson_condition[Karlsson_condition == "0"]="control"
Karlsson_condition[Karlsson_condition == "1"]="T2D"

df1 <- cor_analysis(y_sim[Karlsson_condition == "T2D",], y_imp[Karlsson_condition == "T2D",])

df2 <- cor_analysis(y_sim[Karlsson_condition == "control",], y_imp[Karlsson_condition == "control",])

df3 <- cor_analysis(y_sim, y_imp)

df <- df_gen_v2(df3, df1, df2)
df_spearman <- df_gen_v3(df3, df1, df2)

pdf("/Users/hanxinyu/Desktop/real_data_analysis/Karlsson/Karlsson_correlation.pdf", width = 10, height = 5)
ggplot(df, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) +  
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "T2D", "control"),
    labels = c("all", "T2D", "control")
  ) +
  ylab("correlation") +
  ggtitle("Karlsson et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),       
    axis.text.x = element_text(size = 20),    
    axis.title = element_text(size = 20),      
    plot.title = element_text(size = 22, hjust = 0.5), 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Karlsson/Karlsson_correlation_spearman.pdf", width = 10, height = 5)
ggplot(df_spearman, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) +  
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "T2D", "control"),
    labels = c("all", "T2D", "control")
  ) +
  ylab("correlation") +
  ggtitle("Karlsson et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),       
    axis.text.x = element_text(size = 20),     
    axis.title = element_text(size = 20),      
    plot.title = element_text(size = 22, hjust = 0.5), 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()



df_gen_v2 <- function(df1, df2, df3){
  data_tp <- c(rep("all" , 2) , rep("CRC", 2) , rep("control", 2))
  comp <- rep( c("raw_nz", "imp_nz") , 3)
  value <- c(as.vector(df1[,1]), as.vector(df2[,1]), as.vector(df3[,1]))
  data <- data.frame(data_tp, comp, value)
  data$data_tp <- factor(data$data_tp, levels = c("all", "CRC", "control"))
  data$comp <- factor(data$comp, levels = c("raw_nz", "imp_nz"))
  return(data)
} 

df_gen_v3 <- function(df1, df2, df3){
  data_tp <- c(rep("all" , 2) , rep("CRC", 2) , rep("control", 2))
  comp <- rep( c("raw_nz", "imp_nz"), 3)
  value <- c(as.vector(df1[,3]), as.vector(df2[,3]), as.vector(df3[,3]))
  data <- data.frame(data_tp, comp, value)
  data$data_tp <- factor(data$data_tp, levels = c("all", "CRC", "control"))
  data$comp <- factor(data$comp, levels = c("raw_nz", "imp_nz"))
  return(data)
} 


tolerance <- .Machine$double.eps^0.5
Feng_raw[abs(Feng_raw - log10(1.01)) < tolerance] <- log10(1.01)

y_sim = Feng_raw
y_imp <- Feng_imp

Raw_cor <- cor(y_sim)
imp_cor <- cor(y_imp)

cor(as.vector(Raw_cor), as.vector(imp_cor))
cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( imp_cor[upper.tri(imp_cor)] ))


dim(Raw_cor)[1]

nz_cor <- matrix(nrow = 216, ncol = 216)
for(i in 2:dim(y_imp)[2]){
  for(j in 1:(i-1)){
    nz_idx <- intersect(which(y_sim[,i] != log10(1.01)), which(y_sim[,j] != log10(1.01)))
    nz_cor[j,i] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
    nz_cor[i,j] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
  }
}
nz_cor[is.na(nz_cor)] <- 0

diag(nz_cor) = 1
colnames(nz_cor) = colnames(Raw_cor)
rownames(nz_cor) = colnames(Raw_cor)

Feng_raw_cor = Raw_cor
Feng_imp_cor = imp_cor
Feng_nz_cor = nz_cor

Feng_condition[Feng_condition == "0"]="control"
Feng_condition[Feng_condition == "1"]="CRC"

df1 <- cor_analysis(y_sim[Feng_condition == "CRC",], y_imp[Feng_condition == "CRC",])

df2 <- cor_analysis(y_sim[Feng_condition == "control",], y_imp[Feng_condition == "control",])

df3 <- cor_analysis(y_sim, y_imp)

df <- df_gen_v2(df3, df1, df2)
df_spearman <- df_gen_v3(df3, df1, df2)

pdf("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Feng_correlation.pdf", width = 10, height = 5)
ggplot(df, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) +  
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "CRC", "control"),
    labels = c("all", "CRC", "control")
  ) +
  ylab("correlation") +
  ggtitle("Feng et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),       
    axis.text.x = element_text(size = 20),     
    axis.title = element_text(size = 20),     
    plot.title = element_text(size = 22, hjust = 0.5), 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Feng_correlation_spearman.pdf", width = 10, height = 5)
ggplot(df_spearman, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + 
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "CRC", "control"),
    labels = c("all", "CRC", "control")
  ) +
  ylab("correlation") +
  ggtitle("Feng et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),       
    axis.text.x = element_text(size = 20),     
    axis.title = element_text(size = 20),      
    plot.title = element_text(size = 22, hjust = 0.5),  
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()



tolerance <- .Machine$double.eps^0.5
Yu_raw[abs(Yu_raw - log10(1.01)) < tolerance] <- log10(1.01)

y_sim = Yu_raw
y_imp <- Yu_imp

Raw_cor <- cor(y_sim)
imp_cor <- cor(y_imp)

cor(as.vector(Raw_cor), as.vector(imp_cor))
cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( imp_cor[upper.tri(imp_cor)] ))

dim(Raw_cor)[1]

nz_cor <- matrix(nrow = 199, ncol = 199)
for(i in 2:dim(y_imp)[2]){
  for(j in 1:(i-1)){
    nz_idx <- intersect(which(y_sim[,i] != log10(1.01)), which(y_sim[,j] != log10(1.01)))
    nz_cor[j,i] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
    nz_cor[i,j] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
  }
}
nz_cor[is.na(nz_cor)] <- 0

diag(nz_cor) = 1
colnames(nz_cor) = colnames(Raw_cor)
rownames(nz_cor) = colnames(Raw_cor)

Yu_raw_cor = Raw_cor
Yu_imp_cor = imp_cor
Yu_nz_cor = nz_cor

Yu_condition[Yu_condition == "0"]="control"
Yu_condition[Yu_condition == "1"]="CRC"

df1 <- cor_analysis(y_sim[Yu_condition == "CRC",], y_imp[Yu_condition == "CRC",])

df2 <- cor_analysis(y_sim[Yu_condition == "control",], y_imp[Yu_condition == "control",])

df3 <- cor_analysis(y_sim, y_imp)

df <- df_gen_v2(df3, df1, df2)
df_spearman <- df_gen_v3(df3, df1, df2)

pdf("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Yu_correlation.pdf", width = 10, height = 5)
ggplot(df, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) +  
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "CRC", "control"),
    labels = c("all", "CRC", "control")
  ) +
  ylab("correlation") +
  ggtitle("Yu et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),       
    axis.text.x = element_text(size = 20),    
    axis.title = element_text(size = 20),     
    plot.title = element_text(size = 22, hjust = 0.5),  
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Yu_correlation_spearman.pdf", width = 10, height = 5)
ggplot(df_spearman, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) +  
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "CRC", "control"),
    labels = c("all", "CRC", "control")
  ) +
  ylab("correlation") +
  ggtitle("Yu et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),      
    axis.text.x = element_text(size = 20),     
    axis.title = element_text(size = 20),      
    plot.title = element_text(size = 22, hjust = 0.5),  
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()



tolerance <- .Machine$double.eps^0.5
Vogtmann_raw[abs(Vogtmann_raw - log10(1.01)) < tolerance] <- log10(1.01)

y_sim = Vogtmann_raw
y_imp <- Vogtmann_imp

Raw_cor <- cor(y_sim)
imp_cor <- cor(y_imp)

cor(as.vector(Raw_cor), as.vector(imp_cor))
cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( imp_cor[upper.tri(imp_cor)] ))


dim(Raw_cor)[1]

nz_cor <- matrix(nrow = 210, ncol = 210)
for(i in 2:dim(y_imp)[2]){
  for(j in 1:(i-1)){
    nz_idx <- intersect(which(y_sim[,i] != log10(1.01)), which(y_sim[,j] != log10(1.01)))
    nz_cor[j,i] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
    nz_cor[i,j] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
  }
}
nz_cor[is.na(nz_cor)] <- 0

diag(nz_cor) = 1
colnames(nz_cor) = colnames(Raw_cor)
rownames(nz_cor) = colnames(Raw_cor)

Vogtmann_raw_cor = Raw_cor
Vogtmann_imp_cor = imp_cor
Vogtmann_nz_cor = nz_cor

Vogtmann_condition[Vogtmann_condition == "0"]="control"
Vogtmann_condition[Vogtmann_condition == "1"]="CRC"

cor_analysis <- function(y_sim, y_imp){
  Raw_cor <- cor(y_sim)
  print(Raw_cor)
  imp_cor <- cor(y_imp)
  
  cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( imp_cor[upper.tri(imp_cor)] ))

  nz_cor <- matrix(nrow = dim(Raw_cor)[1], ncol = dim(Raw_cor)[1])
  for(i in 2:dim(y_imp)[2]){
    for(j in 1:(i-1)){
      nz_idx <- intersect(which(y_sim[,i] != log10(1.01)), which(y_sim[,j] != log10(1.01)))
      nz_cor[j,i] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
    }
  }
  nz_cor[is.na(nz_cor)] <- 0
  
  cor_result <- matrix(nrow = 2, ncol = 4)
  cor_result[1,1] = cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( nz_cor[upper.tri(nz_cor)] ))
  print(cor_result[1,1])
  cor_result[2,1] = cor(as.vector( imp_cor[upper.tri(imp_cor)] ), as.vector( nz_cor[upper.tri(nz_cor)] ))
  cor_result[2,2] = cor(as.vector( Raw_cor[upper.tri(imp_cor)] ), as.vector( imp_cor[upper.tri(nz_cor)] ))
  cor_result[1,3] = cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( nz_cor[upper.tri(nz_cor)] ), method = "spearman")
  cor_result[2,3] = cor(as.vector( imp_cor[upper.tri(imp_cor)] ), as.vector( nz_cor[upper.tri(nz_cor)] ), method = "spearman")
  cor_result[2,4] = cor(as.vector( Raw_cor[upper.tri(imp_cor)] ), as.vector( imp_cor[upper.tri(nz_cor)] ), method = "spearman")
  return(cor_result)
}

a <- y_sim[Vogtmann_condition == "CRC",]
a1 <- a[1:52,]
b <- y_sim[Vogtmann_condition == "control",]
b1 <- b[1:52,]
c <-y_imp[Vogtmann_condition == "CRC",]
c1 <- c[1:52,]
d <- y_imp[Vogtmann_condition == "control",]
d1 <- d[1:52,]
df1 <- cor_analysis(a1, c1)


df2 <- cor_analysis(b1, d1)

df3 <- cor_analysis(y_sim, y_imp)

df <- df_gen_v2(df3, df1, df2)
df_spearman <- df_gen_v3(df3, df1, df2)

pdf("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Vogtmann_correlation.pdf", width = 10, height = 5)
ggplot(df, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + 
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "CRC", "control"),
    labels = c("all", "CRC", "control")
  ) +
  ylab("correlation") +
  ggtitle("Vogtmann et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),   
    axis.text.x = element_text(size = 20),     
    axis.title = element_text(size = 20),     
    plot.title = element_text(size = 22, hjust = 0.5),  
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Vogtmann_correlation_spearman.pdf", width = 10, height = 5)
ggplot(df_spearman, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + 
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "CRC", "control"),
    labels = c("all", "CRC", "control")
  ) +
  ylab("correlation") +
  ggtitle("Vogtmann et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),       
    axis.text.x = element_text(size = 20),     
    axis.title = element_text(size = 20),      
    plot.title = element_text(size = 22, hjust = 0.5), 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()



tolerance <- .Machine$double.eps^0.5
Zeller_raw[abs(Zeller_raw - log10(1.01)) < tolerance] <- log10(1.01)

y_sim = Zeller_raw
y_imp <- Zeller_imp

Raw_cor <- cor(y_sim)
imp_cor <- cor(y_imp)

cor(as.vector(Raw_cor), as.vector(imp_cor))
cor(as.vector( Raw_cor[upper.tri(Raw_cor)] ), as.vector( imp_cor[upper.tri(imp_cor)] ))


dim(Raw_cor)[1]

nz_cor <- matrix(nrow = 237, ncol = 237)
for(i in 2:dim(y_imp)[2]){
  for(j in 1:(i-1)){
    nz_idx <- intersect(which(y_sim[,i] != log10(1.01)), which(y_sim[,j] != log10(1.01)))
    nz_cor[j,i] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
    nz_cor[i,j] <- cor(y_sim[nz_idx,i],y_sim[nz_idx,j])
  }
}
nz_cor[is.na(nz_cor)] <- 0

diag(nz_cor) = 1
colnames(nz_cor) = colnames(Raw_cor)
rownames(nz_cor) = colnames(Raw_cor)

Zeller_raw_cor = Raw_cor
Zeller_imp_cor = imp_cor
Zeller_nz_cor = nz_cor

Zeller_condition[Zeller_condition == "0"]="control"
Zeller_condition[Zeller_condition == "1"]="CRC"

df1 <- cor_analysis(y_sim[Zeller_condition == "CRC",], y_imp[Zeller_condition == "CRC",])

df2 <- cor_analysis(y_sim[Zeller_condition == "control",], y_imp[Zeller_condition == "control",])

df3 <- cor_analysis(y_sim, y_imp)

df <- df_gen_v2(df3, df1, df2)
df_spearman <- df_gen_v3(df3, df1, df2)

pdf("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Zeller_correlation.pdf", width = 10, height = 5)
ggplot(df, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) +  
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "CRC", "control"),
    labels = c("all", "CRC", "control")
  ) +
  ylab("correlation") +
  ggtitle("Zeller et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),       
    axis.text.x = element_text(size = 20),    
    axis.title = element_text(size = 20),      
    plot.title = element_text(size = 22, hjust = 0.5), 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Yu/Zeller_correlation_spearman.pdf", width = 10, height = 5)
ggplot(df_spearman, aes(fill= comp, y=value, x=data_tp)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) +  
  scale_fill_manual(values = c("raw_nz" = alpha("#636363", 0.5), "imp_nz" = "#636363")) +
  theme_bw() +
  scale_x_discrete(
    breaks = c("all", "CRC", "control"),
    labels = c("all", "CRC", "control")
  ) +
  ylab("correlation") +
  ggtitle("Zeller et al.") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),      
    axis.text.x = element_text(size = 20),    
    axis.title = element_text(size = 20),      
    plot.title = element_text(size = 22, hjust = 0.5),  
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()



Karlsson_control_ori <- Karlsson_raw[Karlsson_condition == "control",]
print(colnames(Karlsson_control_ori))
Karlsson_control_imputed <- Karlsson_imp[Karlsson_condition == "control",]
colnames(Karlsson_control_ori) <- gsub(".*s__", "s__", colnames(Karlsson_control_ori))
print(colnames(Karlsson_control_ori))
colnames(Karlsson_control_imputed) <- gsub(".*s__", "s__", colnames(Karlsson_control_imputed))


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Karlsson/Karlsson_CON_R.torques_E.eligens_raw.pdf", width = 7, height = 7)

plot(y = Karlsson_control_ori[,"s__Ruminococcus_torques"], x = Karlsson_control_ori[,"s__Eubacterium_eligens"],
     xlab = "Eubacterium_eligens", ylab = "Ruminococcus_troques", main = "Karlsson control data", 
     xlim = c(0, 6))
nz_idx <- intersect(which(Karlsson_control_ori[,"s__Ruminococcus_torques"] > log10(1.02)), which(Karlsson_control_ori[,"s__Eubacterium_eligens"] > log10(1.02)))
abline(coef = lm(Karlsson_control_ori[,"s__Ruminococcus_torques"] ~ Karlsson_control_ori[,"s__Eubacterium_eligens"])$coefficients, col = "black", lwd = 4)
abline(coef = lm(Karlsson_control_ori[nz_idx,"s__Ruminococcus_torques"] ~ Karlsson_control_ori[nz_idx,"s__Eubacterium_eligens"])$coefficients, col = alpha("blue", 0.5), lwd = 4)
legend(x='bottomleft', 
       legend= c(paste("Cor = ", round(cor(Karlsson_control_ori[,"s__Ruminococcus_torques"], 
                                           Karlsson_control_ori[,"s__Eubacterium_eligens"]), digits = 2), sep = ""),
                 paste("nz_Cor = ", round(cor(Karlsson_control_ori[nz_idx,"s__Ruminococcus_torques"], 
                                              Karlsson_control_ori[nz_idx,"s__Eubacterium_eligens"]), digits = 2), sep = "")),
       lty = c(1,1), 
       cex = 1, 
       col=c("black", "blue"))

dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Karlsson/Karlsson_CON_R.torques_E.eligens_imp1.pdf", width = 7, height = 7)

plot(y = Karlsson_control_imputed[,"s__Ruminococcus_torques"], x = Karlsson_control_imputed[,"s__Eubacterium_eligens"], 
     xlab = "Eubacterium_eligens", ylab = "Ruminococcus_troques", main = "TphPMF + Karlsson control data", 
     xlim = c(0, 6), ylim=c(0,5))
nz_idx <- intersect(which(Karlsson_control_imputed[,"s__Ruminococcus_torques"] > log10(1.02)), which(Karlsson_control_imputed[,"s__Eubacterium_eligens"] > log10(1.02)))
abline(coef = lm(Karlsson_control_imputed[,"s__Ruminococcus_torques"] ~ Karlsson_control_imputed[,"s__Eubacterium_eligens"])$coefficients, col = alpha("blue", 0.5), lwd = 4)
legend(x='bottomleft', legend= paste("Cor = ", round(cor(Karlsson_control_imputed[,"s__Ruminococcus_torques"], Karlsson_control_imputed[,"s__Eubacterium_eligens"]), digits = 2), sep = ""), lty = 1, cex = 0.8, col = alpha("blue", 0.5))
dev.off()



Karlsson_T2D_ori <- Karlsson_raw[Karlsson_condition == "T2D",]
Karlsson_T2D_imputed <- Karlsson_imp[Karlsson_condition == "T2D",]
colnames(Karlsson_T2D_ori) <- gsub(".*s__", "s__", colnames(Karlsson_T2D_ori))
colnames(Karlsson_T2D_imputed) <- gsub(".*s__", "s__", colnames(Karlsson_T2D_imputed))


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Karlsson/Karlsson_T2D_R.torques_E.eligens_raw.pdf", width = 7, height =7)
plot(y = Karlsson_T2D_ori[,"s__Ruminococcus_torques"], x = Karlsson_T2D_ori[,"s__Eubacterium_eligens"],
     xlab = "Eubacterium_eligens", ylab = "Ruminococcus_troques", main = "Karlsson T2D data", 
     xlim = c(0, 6))
nz_idx <- intersect(which(Karlsson_T2D_ori[,"s__Ruminococcus_torques"] > log10(1.02)), which(Karlsson_T2D_ori[,"s__Eubacterium_eligens"] > log10(1.02)))
abline(coef = lm(Karlsson_T2D_ori[,"s__Ruminococcus_torques"] ~ Karlsson_T2D_ori[,"s__Eubacterium_eligens"])$coefficients, col = "black", lwd = 4)
abline(coef = lm(Karlsson_T2D_ori[nz_idx,"s__Ruminococcus_torques"] ~ Karlsson_T2D_ori[nz_idx,"s__Eubacterium_eligens"])$coefficients, col = alpha("blue", 0.5), lwd = 4)
legend(x='bottomleft', legend= c( paste("Cor = ", round(cor(Karlsson_T2D_ori[,"s__Ruminococcus_torques"], Karlsson_T2D_ori[,"s__Eubacterium_eligens"]), digits = 2), sep = ""),                                  paste("nz_Cor = ", round(cor(Karlsson_T2D_ori[nz_idx,"s__Ruminococcus_torques"], Karlsson_T2D_ori[nz_idx,"s__Eubacterium_eligens"]), digits = 2), sep = "")), lty = c(1,1), cex = 0.8, col=c("black", "blue"))
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/Karlsson/Karlsson_T2D_R.torques_E.eligens_imp1.pdf", width = 7, height =7)

plot(y = Karlsson_T2D_imputed[,"s__Ruminococcus_torques"], x = Karlsson_T2D_imputed[,"s__Eubacterium_eligens"], 
     xlab = "Eubacterium_eligens", ylab = "Ruminococcus_troques", main = "TphPMF + Karlsson T2D data", 
     xlim = c(0, 6), ylim=c(0,5))
nz_idx <- intersect(which(Karlsson_T2D_imputed[,"s__Ruminococcus_torques"] > log10(1.02)), which(Karlsson_T2D_imputed[,"s__Eubacterium_eligens"] > log10(1.02)))
abline(coef = lm(Karlsson_T2D_imputed[,"s__Ruminococcus_torques"] ~ Karlsson_T2D_imputed[,"s__Eubacterium_eligens"])$coefficients, col = alpha("blue", 0.5), lwd = 4)
legend(x='bottomleft', legend= paste("Cor = ", round(cor(Karlsson_T2D_imputed[,"s__Ruminococcus_torques"], Karlsson_T2D_imputed[,"s__Eubacterium_eligens"]), digits = 2), sep = ""), lty = 1, cex = 0.8, col = alpha("blue", 0.5))

dev.off()


Qin_control_ori <- Qin_raw[Qin_condition == "control",]
Qin_control_imputed <- Qin_imp[Qin_condition == "control",]
colnames(Qin_control_ori) <- gsub(".*s__", "s__", colnames(Qin_control_ori))
colnames(Qin_control_imputed) <- gsub(".*s__", "s__", colnames(Qin_control_imputed))

print(Qin_control_ori[,"s__Dorea_formicigenerans"])

pdf("/Users/hanxinyu/Desktop/real_data_analysis/Qin/Qin_CON_B.obeum_A.putredinis_raw.pdf", width = 7, height =7)
plot(y = Qin_control_ori[,"s__Blautia_obeum"], x = Qin_control_ori[,"s__Alistipes_putredinis"], 
     xlab = "Alistipes_putredinise", ylab = "Blautia_obeum", main = "Qin control data", 
     xlim = c(0, 6), ylim=c(0,5))
nz_idx <- intersect(which(Qin_control_ori[,"s__Blautia_obeum"] > log10(1.02)), which(Qin_control_ori[,"s__Alistipes_putredinis"] > log10(1.02)))
abline(coef = lm(Qin_control_ori[,"s__Blautia_obeum"] ~ Qin_control_ori[,"s__Alistipes_putredinis"])$coefficients, lwd = 4)
abline(coef = lm(Qin_control_ori[nz_idx,"s__Blautia_obeum"] ~ Qin_control_ori[nz_idx,"s__Alistipes_putredinis"])$coefficients, lty = 1, col = alpha("blue", 0.5), lwd = 4)
legend(x='bottomright', legend= c( paste("Cor = ", round(cor(Qin_control_ori[,"s__Blautia_obeum"], Qin_control_ori[,"s__Alistipes_putredinis"]), digits = 2), sep = ""),                                  paste("nz_Cor = ", round(cor(Qin_control_ori[nz_idx,"s__Blautia_obeum"], Qin_control_ori[nz_idx,"s__Alistipes_putredinis"]), digits = 2), sep = "")), lty = c(1,1), cex = 0.8, col=c("black", "blue"))
dev.off()


pdf("/Users/hanxinyu/Desktop/real_data_analysis/Qin/Qin_CON_B.obeum_A.putredinis_imp1.pdf", width = 7, height =7)

plot(y = Qin_control_imputed[,"s__Blautia_obeum"], x = Qin_control_imputed[,"s__Alistipes_putredinis"],
     xlab = "Alistipes_putredinise", ylab = "Blautia_obeum", main = "TphPMF + Qin control data", 
     xlim = c(0, 6), ylim=c(0,5))
nz_idx <- intersect(which(Qin_control_imputed[,"s__Blautia_obeum"] > log10(1.02)), which(Qin_control_imputed[,"s__Alistipes_putredinis"] > log10(1.02)))
abline(coef = lm(Qin_control_imputed[,"s__Blautia_obeum"] ~ Qin_control_imputed[,"s__Alistipes_putredinis"])$coefficients, col = alpha("blue", 0.5), lwd = 4)
legend(x='bottomright', legend= paste("Cor = ", round(cor(Qin_control_imputed[,"s__Blautia_obeum"], Qin_control_imputed[,"s__Alistipes_putredinis"]), digits = 2), sep = ""), lty = 1, cex = 0.8, col = alpha("blue", 0.5))
dev.off()

Qin_T2D_ori <- Qin_raw[Qin_condition == "T2D",]
Qin_T2D_imputed <- Qin_imp[Qin_condition == "T2D",]
colnames(Qin_T2D_ori) <- gsub(".*s__", "s__", colnames(Qin_T2D_ori))
colnames(Qin_T2D_imputed) <- gsub(".*s__", "s__", colnames(Qin_T2D_imputed))

pdf("/Users/hanxinyu/Desktop/real_data_analysis/Qin/Qin_T2D_B.obeum_A.putredinis_raw.pdf", width = 7, height =7)
plot(y = Qin_T2D_ori[,"s__Blautia_obeum"], x = Qin_T2D_ori[,"s__Alistipes_putredinis"], 
     xlab = "Alistipes_putredinise", ylab = "Blautia_obeum", main = "Qin T2D data", 
     xlim = c(0, 6), ylim=c(0,5))
nz_idx <- intersect(which(Qin_T2D_ori[,"s__Blautia_obeum"] > log10(1.02)), which(Qin_T2D_ori[,"s__Alistipes_putredinis"] > log10(1.02)))
abline(coef = lm(Qin_T2D_ori[,"s__Blautia_obeum"] ~ Qin_T2D_ori[,"s__Alistipes_putredinis"])$coefficients, lwd = 4)
abline(coef = lm(Qin_T2D_ori[nz_idx,"s__Blautia_obeum"] ~ Qin_T2D_ori[nz_idx,"s__Alistipes_putredinis"])$coefficients, lty = 1, col = alpha("blue", 0.5), lwd = 4)
legend(x='bottomright', legend= c( paste("Cor = ", round(cor(Qin_T2D_ori[,"s__Blautia_obeum"], Qin_T2D_ori[,"s__Alistipes_putredinis"]), digits = 2), sep = ""),                                  paste("nz_Cor = ", round(cor(Qin_T2D_ori[nz_idx,"s__Blautia_obeum"], Qin_T2D_ori[nz_idx,"s__Alistipes_putredinis"]), digits = 2), sep = "")), lty = c(1,1), cex = 0.8, col=c("black", "blue"))
dev.off()

pdf("/Users/hanxinyu/Desktop/real_data_analysis/Qin/Qin_T2D_B.obeum_A.putredinis_imp1.pdf", width = 7, height =7)

plot(y = Qin_T2D_imputed[,"s__Blautia_obeum"], x = Qin_T2D_imputed[,"s__Alistipes_putredinis"], 
     xlab = "Alistipes_putredinise", ylab = "Blautia_obeum", main = "TphPMF + Qin T2D data", 
     xlim = c(0, 6), ylim=c(0,5))
nz_idx <- intersect(which(Qin_T2D_imputed[,"s__Blautia_obeum"] > log10(1.02)), which(Qin_T2D_imputed[,"s__Alistipes_putredinis"] > log10(1.02)))
abline(coef = lm(Qin_T2D_imputed[,"s__Blautia_obeum"] ~ Qin_T2D_imputed[,"s__Alistipes_putredinis"])$coefficients, col = alpha("blue", 0.5), lwd = 4)
legend(x='bottomright', legend= paste("Cor = ", round(cor(Qin_T2D_imputed[,"s__Blautia_obeum"], Qin_T2D_imputed[,"s__Alistipes_putredinis"]), digits = 2), sep = ""), lty = 1, cex = 0.8, col = alpha("blue", 0.5))
dev.off()












