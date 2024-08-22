## -------------------------------------------------------- ##
## ------------------- USER PARAMETERS -------------------- ##
## -------------------------------------------------------- ##

#Directory containing vcf file
Input_Dir <- "/Users/leonedmiidz/Documents/6970/Assignment_4"

vcf_files <- c("lct.vcf", "fut3.vcf", "mcm6.vcf", "fut2.vcf")

#Super-population codes of interest for analysis

pop_list <- c("YRI", "GBR", "FIN", "ESN")

## -------------------------------------------------------- ##
## --------------------- 1.LIBRARIES ---------------------- ##
## -------------------------------------------------------- ##

#Function to load or install packages
install_if_needed <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE) #Explicitly installing all dependencies
    library(package_name, character.only = TRUE)
  }
}

#Create a vector with the names of required packages
packages <- c("tidyverse", "ggplot2", "readxl", "broom", "data.table", "rentrez", "writexl", "vcfR", "tibble", "plotly", "class", "pROC", "randomForest","GGally", "cluster", "pheatmap", "factoextra","caret")

#Iterate over the vector to install/load each package
for (pkg in packages) {
  install_if_needed(pkg)
}

## -------------------------------------------------------- ##
## ------------------------- FILES ------------------------ ##
## -------------------------------------------------------- ##
setwd(Input_Dir)

#Read in 1000 Genomes samples population and gender labels
samples <- read.delim("samples.txt", sep = "\t")

#Filter for only super populations specified above, unless filter is set to FALSE
  samples <- samples %>%
    filter(pop %in% pop_list)

samp_ref_table <- samples %>% #Store a reference table for use later on
  select(sample, pop, super_pop, gender)

#Initilaist lists to store vcf info and genotype dataframes
vcf_data_list <- list()
genotypes_list <- list()

#Function to read in and concatonate vcf object information with vcfR
process_vcf_file <- function(vcf_filename) {
  #reading in a vcf file
  vcf <- read.vcfR(vcf_filename)
  
  #Extracting tidy objects
  tidy_vcf <- vcfR2tidy(vcf, single_frame = TRUE, info_types = TRUE, format_types = TRUE)
  VCF_long <- tidy_vcf$dat
  
  #Filter the VCF info based on filtered samples file
  VCF_long_filtered <- VCF_long %>%
    filter(Indiv %in% samples$sample)
  
  #Extract genotype information
  genotypes <- as.data.frame(extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = TRUE))
  genotypes <- rownames_to_column(genotypes, var = "Variant")
  
  #Modify 'Variant' column to only have numeric positions, handling cases with preceeding and trailing underscores
  genotypes$Variant <- gsub("[0-3000]+_", "", genotypes$Variant) 
  genotypes$Variant <- gsub("_+[0-3000]+", "", genotypes$Variant)
  
  #Return a list containing both pieces of data
  list(VCF_long_filtered = VCF_long_filtered, genotypes = genotypes)
}

#Apply the function to each VCF file
results_list <- lapply(vcf_files, process_vcf_file)

#Combine the lists now
vcf_data_list <- setNames(lapply(results_list, `[[`, "VCF_long_filtered"), vcf_files)
genotypes_list <- setNames(lapply(results_list, `[[`, "genotypes"), vcf_files)

#And combine lists into single dfs
genotypes <- bind_rows(genotypes_list)
VCF_long <-  bind_rows(vcf_data_list)

## ----------------- ##
##  QAULITY CONTROL  ##
## ----------------- ##

#Initial filter based on a few parameters
VCF_long <- VCF_long %>%
  filter(FILTER == "PASS", 
         VT == "SNP", #Only looking at SNPs
         MULTI_ALLELIC == "FALSE", #Remove multi-allelic variants
         nchar(REF) == 1, #Ensure only single REF and ALT Allele is listed
         nchar(ALT) == 1) 

#genotypes <- rownames_to_column(genotypes, var = "Variant") 
genotypes$Variant <- sub(".*_", "", genotypes$Variant)

#Filter genotype df to include only variants and samples in the VCF_long df after initial filtering
genotypes <- genotypes %>%
  filter(Variant %in% VCF_long$POS) %>%
  select(Variant, all_of(samples$sample)) #Include only sample cols present in filtered sample df

#Convert to numeric
convert_genotype <- function(gt) {
  ifelse(gt == "0|0", 0,
         ifelse(gt %in% c("0|1", "1|0"), 1,
                ifelse(gt == "1|1", 2, NA_integer_)))
}

genotypes[2:ncol(genotypes)] <- lapply(genotypes[2:ncol(genotypes)], convert_genotype)

genotypes <- genotypes %>%
  mutate(across(-1, as.numeric))

#Now compute HWE parameters using genotype df, to further filter variants

hwe <- data.frame(Variant = genotypes$Variant)

Total_individuals <- ncol(genotypes) - 1

#Calculate the sum for each genotype (AA, Aa, aa) based on counts 0, 1, or 2
hwe$Obs_AA <- rowSums(genotypes[, -1] == 0)
hwe$Obs_Aa <- rowSums(genotypes[, -1] == 1)
hwe$Obs_aa <- rowSums(genotypes[, -1] == 2)
hwe$AA_prop <- hwe$Obs_AA / Total_individuals
hwe$Aa_prop <- hwe$Obs_Aa / Total_individuals  
hwe$aa_prop <- hwe$Obs_aa / Total_individuals

#First, filter our variants with q < 0.001.

hwe$p_hat <- (2 * hwe$Obs_AA + hwe$Obs_Aa) / (2 * Total_individuals)
hwe$`1-p` <- (2 * hwe$Obs_aa + hwe$Obs_Aa) / (2 * Total_individuals)

hwe <- hwe %>%
  filter(`1-p` >= 0.001) #Filter for alternative allele freq less than 0.1%

#Expected genotype frequencies under HWE
hwe$Exp_AA <- Total_individuals * (hwe$p_hat^2)
hwe$Exp_Aa <- Total_individuals * 2 * hwe$p_hat * (1 - hwe$p_hat)
hwe$Exp_aa <- Total_individuals * ((1 - hwe$p_hat)^2)

#Calculate chi-square statistics for HWE
hwe$Chi_Square <- (hwe$Obs_AA - hwe$Exp_AA)^2 / hwe$Exp_AA +
  (hwe$Obs_Aa - hwe$Exp_Aa)^2 / hwe$Exp_Aa +
  (hwe$Obs_aa - hwe$Exp_aa)^2 / hwe$Exp_aa

#Calculate p-values from the chi-square statistics
hwe$P_Value <- pchisq(hwe$Chi_Square, df = 1, lower.tail = FALSE)

#Determine if the locus deviates significantly from HWE
alpha <- 0.05
hwe$HWE_Deviation <- hwe$P_Value < alpha

hwe <- hwe %>%
  filter(HWE_Deviation == "FALSE")

#We will not filter based on LD, since ll variants are on the same gene and we expect many of them will be in LD. 

genotypes <- genotypes %>%
  filter(Variant %in% hwe$Variant)

genotype_matrix <- as.matrix(genotypes) #Convert to matrix for further analysis

## -------------------------------------------------------- ##
## -------------------- CLUSTERS PREP --------------------- ##
## -------------------------------------------------------- ##

#Transpose the matrix to have individuals as rows and SNPs as columns for PCA
genotype_matrix <- t(genotype_matrix)
colnames(genotype_matrix) <- genotype_matrix[1, ]
genotype_matrix <- genotype_matrix[-1, ]

#Ensure matrix is numeric before proceeding
original_row_names <- rownames(genotype_matrix) #Preserve sample IDs
original_col_names <- colnames(genotype_matrix) #preserve variant POS

#Now convert to numeric
genotype_matrix <- apply(genotype_matrix, 2, as.numeric)

#Replace NA with the median 
genotype_matrix[is.na(genotype_matrix)] <- apply(genotype_matrix, 2, median, na.rm = TRUE)

# Replace NA with the median of each column
# This loop goes through each column individually
for(i in 1:ncol(genotype_matrix)) {
  na_indices <- is.na(genotype_matrix[, i])
  if (any(na_indices)) {  # Only calculate the median if there are NAs
    genotype_median <- median(genotype_matrix[, i], na.rm = TRUE)
    genotype_matrix[na_indices, i] <- genotype_median
  }
}

#Remove columns that have no variation (are constant)
non_constant_columns <- apply(genotype_matrix, 2, var) != 0
genotype_matrix_non_constant <- genotype_matrix[, non_constant_columns]

#Bring back the col and row names
rownames(genotype_matrix_non_constant) <- original_row_names
colnames(genotype_matrix_non_constant) <- original_col_names[non_constant_columns]

#Exploration of trends in the cleaned data
biplot(prcomp(genotype_matrix_non_constant, scale=T))
fviz_dist(as.dist(cor(genotype_matrix_non_constant)), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
fviz_dist(dist(genotype_matrix_non_constant), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

## -------------------------------------------------------- ##
## ------------------------- PCA -------------------------- ##
## -------------------------------------------------------- ##

pca_result <- prcomp(genotype_matrix_non_constant, center = TRUE, scale. = TRUE)

#Summary of PCA results
pca_summary <- summary(pca_result)

#Extract the eigenvalues for the scree plot.
eigenvalues <- pca_result$sdev^2

#Construct the scree plot.
plot(1:length(eigenvalues), eigenvalues, 
     type = "b", 
     xlab = "Principal Component", 
     ylab = "Eigenvalue",
     main = "Scree Plot")

#Compute variance explained for first 10 PCAs
var_explained <- pca_summary$importance[2, ] * 100

#Create a data frame from the PCA result for plotting
pca_data <- data.frame(PC1 = pca_result$x[, 1],
                       PC2 = pca_result$x[, 2])

#Add in population and superpopulation labels to the pca_data df
pca_data$sample <- rownames(pca_data) #Convert rownames to col
#Now merge with filtered sample df by sample col
pca_data <- merge(pca_data, samp_ref_table, by = "sample")


### PLOT PREPARATION ###

#Plot the PCA representing pop with colour and super population with shape
Super_pop <- ggplot(pca_data, aes(x = PC1, y = PC2, color = super_pop)) +
  geom_point(size = 0.35) +  
  theme_minimal() + 
  labs(title = "PCA by Super Population",
       x = paste("Principal Component 1 (", sprintf("%.2f%%", var_explained[1]), ")", sep=""),
       y = paste("Principal Component 2 (", sprintf("%.2f%%", var_explained[2]), ")", sep=""),
       color = "Population") +
  scale_color_brewer(palette = "Set1")

#Convert to plotly
Super_pop_PCA <- ggplotly(Super_pop)

Super_pop_PCA <- Super_pop_PCA %>%
  layout(
    plot_bgcolor = 'white',  #Background color
    legend = list(title = list(text = 'Super Population')),
    hovermode = 'closest'
  )

sub_pop_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(size = 0.35) + 
  theme_minimal() +
  labs(title = "PCA by Sub Population",
       x = paste("Principal Component 1 (", sprintf("%.2f%%", var_explained[1]), ")", sep=""),
       y = paste("Principal Component 2 (", sprintf("%.2f%%", var_explained[2]), ")", sep=""),
       color = "Population") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right") 

#Convert to plotly
sub_pop_PCA <- ggplotly(sub_pop_plot)

sub_pop_PCA <- sub_pop_PCA %>%
  layout(
    plot_bgcolor = 'white',  #background color
    legend = list(title = list(text = 'Super Population')),
    hovermode = 'closest'
  )

# ---VIEW PLOTS ----#

Super_pop_PCA
sub_pop_PCA

## -------------------------------------------------------- ##
## ---------------- HIERARCHICAL CLUSTERING --------------- ##
## -------------------------------------------------------- ##

#Compute Manhattan distance matrix for clustering
manhattan_dist <- dist(genotype_matrix_non_constant, method = "manhattan")

#Divisive cluster plot
plot(diana(manhattan_dist),which.plots=2, main="Divisive Clutering", cex = 0.5)

#Perform agglomerative clustering with Manhattan distance and single linkage.
plot(agnes(manhattan_dist, method = "single"), main="Single Linkage Clustering", which.plots=2, cex = 0.5)

#Repeat with average linkage
plot(agnes(manhattan_dist, method = "average"), main="Average Linkage Clustering", which.plots=2, cex = 0.5)

#Repeat with complete linkage
plot(agnes(manhattan_dist, method = "complete"), main="Complete Linkage Clustering", which.plots=2, cex = 0.5)
abline(h=120, col='red')

#Complete linkage appears to have clustered the best. Cut the tree to two clusters and analyze the split of the data.
cut_complete <- cutree(hclust(manhattan_dist, method = "complete"), h=120)

#Create a dataframe with population information to add the clustering grouping info
data_to_compare <- as.data.frame(genotype_matrix_non_constant)
data_to_compare$sample <- rownames(data_to_compare)
data_to_compare <- merge(data_to_compare, samp_ref_table, by = "sample")
rownames(data_to_compare) <- original_row_names

#Add the cluster grouping info to the data.
data_to_compare$cutcomp <- cut_complete

#Check to see how the data were clustered.
aggregate(cbind(pop, super_pop)~cutcomp,data_to_compare,table)

## -------------------------------------------------------- ##
## ------------------ KNN CLASSIFICATION ------------------ ##
## -------------------------------------------------------- ##

#Pre-process data: Scale the genotype matrix
genotypes_labeled <- merge(genotype_matrix_non_constant, samp_ref_table[, c("sample", "super_pop")], by.x = "row.names", by.y = "sample")

rownames(genotypes_labeled) <- genotypes_labeled$Row.names #Restore rownames as sample IDs
genotypes_labeled <- genotypes_labeled[,-1] #Remove redundant columns

genotypes_scaled <- scale(genotype_matrix_non_constant[, -ncol(genotype_matrix_non_constant)])

#Prepare data frame for KNN
data_for_knn <- as.data.frame(genotypes_scaled)
data_for_knn$super_pop <- samp_ref_table$super_pop[match(rownames(data_for_knn), samp_ref_table$sample)]

set.seed(346)  #Ensures reproducibility in your sample

#Splitting the scaled and labeled dataset
samp_size <- floor(0.75 * nrow(data_for_knn))
train_ind <- sample(seq_len(nrow(data_for_knn)), size = samp_size)
train_data <- data_for_knn[train_ind, ]
test_data <- data_for_knn[-train_ind, ]

#Preparing labels and data
train_labels <- train_data$super_pop
test_labels <- test_data$super_pop
train_data <- train_data[,-ncol(train_data)]
test_data <- test_data[,-ncol(test_data)]

#Applying KNN
knn_result <- knn(train = train_data, test = test_data, cl = train_labels, k = 5)

#Ensure test_labels and knn_result are factors with the same levels
test_labels <- factor(test_labels, levels = unique(genotypes_labeled$super_pop))
knn_result <- factor(knn_result, levels = levels(test_labels))

#Now add the loop to test different k values and calculate accuracies
k_values <- 1:20
accuracies <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  set.seed(346)
  knn_result_k <- knn(train = train_data, test = test_data, cl = train_labels, k = k_values[i])
  knn_result_k <- factor(knn_result_k, levels = levels(test_labels))
  
  #Compute accuracy
  conf_matrix_k <- table(Predicted = knn_result_k, Actual = test_labels)
  accuracies[i] <- sum(diag(conf_matrix_k)) / sum(conf_matrix_k)
}

#Now create the plot
accuracy_data <- data.frame(k = k_values, Accuracy = accuracies)
plot(accuracy_data$k, accuracy_data$Accuracy, type = 'b', col = 'blue',
     xlab = 'Number of Neighbors (k)', ylab = 'Accuracy',
     main = 'KNN Accuracy for Different k Values')

#Evaluation Metrics
conf_matrix <- table(Predicted = knn_result, Actual = test_labels)
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(conf_matrix)
print(paste("Accuracy:", accuracy))

#If using caret for extended metrics (optional)
confusionMatrix(knn_result, test_labels)

## -------------------------------------------------------- ##
## -------------- RANDOM FOREST CLASSIFICATION ------------ ##
## -------------------------------------------------------- ##

#Fit RF model with default parameters and 500 trees
model <- randomForest(x = train_data, y = as.factor(train_labels), ntree = 500)
#Predictions
pred_test_rf <- predict(model, newdata=test_data, type="response")

#Evaluation
rf_acc <- sum(pred_test_rf == test_labels) / length(test_labels)
print(paste("Random Forest Accuracy:", rf_acc))

#Confusion matrix
rf_conf_matrix <- table(Predicted = pred_test_rf, Actual = test_labels)
print(rf_conf_matrix)

#Look at error versus ntrees
plot(model)

#and variable importance
varImpPlot(model)

#Define training control
trcontrol <- trainControl(method="repeatedcv", 
                        number=10, 
                        repeats=1,
                        allowParallel = FALSE)

#Define the tuning grid to sample from
tunegrid <- expand.grid(.mtry=c(1:40))

#Perform the random search
rf_gridsearch_acc <- train(x = train_data, 
                           y = as.factor(train_labels),
                           method = "rf",
                           metric = "Accuracy",
                           trControl = trcontrol,
                           tuneGrid = tunegrid)

#Print and plot results
print(rf_gridsearch_acc)
plot(rf_gridsearch_acc)

#Training a final model with mtry =2
final_rf <- randomForest(x = train_data, y = as.factor(train_labels), mtry = 12, ntree = 500)

#Predictions
final_pred <- predict(final_rf, newdata=test_data, type="response")

rf_final_acc <- sum(final_pred == test_labels) / length(test_labels)
print(paste("Random Forest Accuracy:", rf_final_acc))

#Confusion matrix
final_conf_mtx <- table(Predicted = final_pred, Actual = test_labels)
print(final_conf_mtx)

#Extract error data for plotting
error_data <- as.data.frame(final_rf$err.rate)
error_data$Trees <- seq_along(final_rf$err.rate[,1])

ggplot(error_data, aes(x = Trees)) +
  geom_line(aes(y = OOB, color = "Overall Error")) +
  geom_line(aes(y = AFR, color = "Error AFR")) +
  geom_line(aes(y = EUR, color = "Error EUR")) +
  scale_color_manual(values = c("Overall Error" = "black", "Error AFR" = "purple", "Error EUR" = "orange")) +
  labs(x = "Number of Trees",y = "Error Rate", color = "Error Type",title = "OOB Error Rate by Number of Trees - Final Random Forest Model") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),  legend.text = element_text(size = 12))

#With randomForest we can get feature importance
importance_scores <- importance(final_rf)

#Filter for top 60 and plot
top_60_features <- head(importance_df, 60)
feature_names <- rownames(top_60_features)

par(mgp=c(3,0,0))

barplot(top_60_features$Importance,
        names.arg = feature_names, 
        las = 2,  
        main = "Top 60 SNPs in Final Random Forest Model", 
        col = 'orange', 
        ylab = "Importance (MeanDecreaseGini)", 
        xlab = "SNP Positions",
        cex.names = 0.7) 
