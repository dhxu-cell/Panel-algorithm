# Load necessary libraries
library(NMF)
library(Seurat)
library(dplyr)
library(harmony)
library(GSVA)
library(pheatmap)
library(AUCell)
library(msigdbr)
library(monocle3)
library(ggcor)
library(caret)
library(pROC)
library(CMScaller)

# Load Signature and Expression data
Signature <- readRDS('your/pathway/signature.rds')  # Signature for NMF classification
Exp <- readRDS('your/pathway/exp.rds')  # Expression matrix, rows are genes, columns are samples

# Perform GSVA analysis
GSVA_result <- gsva(expr = as.matrix(Exp), gset.idx.list = Signature, mx.diff = TRUE, method = 'gsva', kcdf = "Gaussian")

# Determine the optimal number of ranks for NMF
ranks <- 3:10
nmf_result <- nmf(GSVA_result, ranks, nrun = 100)
plot(estim)  # Select the optimal rank based on the plot

# Run NMF with the selected rank
optimal_rank <- 5  # Set this based on the plot
nmf_result_final <- nmf(GSVA_result, rank = optimal_rank, nrun = 100, method = "brunet")

# Extract features and group predictions
index <- extractFeatures(nmf_result_final, "max")
sig.order <- unlist(index)
group <- predict(nmf_result_final)
Label <- as.data.frame(group)  # Sample labels

# Split data into training and testing sets
set.seed(123)
trainIndex <- createDataPartition(Label$group, p = 0.7, list = FALSE)
Exp_train <- Exp[, trainIndex]
Exp_test <- Exp[, -trainIndex]
Label_train <- Label[trainIndex, , drop = FALSE]
Label_test <- Label[-trainIndex, , drop = FALSE]

# Set up RFE control
rfeControl <- rfeControl(functions = caretFuncs, method = "cv", number = 10)
threshold_auc <- 0.8
final_auc <- 0
best_features <- NULL

# Define the function for RFE and NTP evaluation
run_rfe_ntp <- function(Exp_train, Label_train, Exp_test, Label_test, threshold_auc) {
  auc <- 1
  final_auc <- 0
  best_features <- NULL
  
  while (auc > threshold_auc && nrow(Exp_train) > 5) {
    # Recursive Feature Elimination (RFE)
    rfeResult <- rfe(x = t(Exp_train), y = Label_train$group, sizes = c(nrow(Exp_train) * 0.9), rfeControl = rfeControl)
    
    # Get the current 90% of gene features
    selected_features <- predictors(rfeResult)
    
    # Update training and testing sets based on selected features
    Exp_train <- Exp_train[selected_features, ]
    Exp_test <- Exp_test[selected_features, ]
    
    # NTP validation
    ntp_result <- CMSclassifier(data = as.matrix(Exp_test), templates = as.matrix(Exp_train), annotation = Label_train$group)
    predictions <- ntp_result$prediction
    
    # Calculate AUC
    roc_obj <- roc(Label_test$group, predictions)
    auc <- auc(roc_obj)
    
    # Update best AUC and feature set
    if (auc > final_auc) {
      final_auc <- auc
      best_features <- selected_features
    }
    
    # Print current AUC and number of genes
    cat("Current AUC:", auc, "with", length(selected_features), "genes\n")
    
    # Check stop conditions
    if (auc < 0.8 || length(selected_features) < 5) {
      break
    }
  }
  
  return(list(auc = final_auc, features = best_features))
}

# Run the feature selection and validation
results <- run_rfe_ntp(Exp_train, Label_train, Exp_test, Label_test, threshold_auc)

# Output the final results
cat("Final AUC:", results$auc, "\n")
cat("Selected features:\n")
print(results$features)

# Save the best feature set
saveRDS(results$features, file = "best_features.rds")
