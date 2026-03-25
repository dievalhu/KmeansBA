# ==========================================
# K-MedoidsSC Optimized Version
# ==========================================

library(cluster)
library(proxy)
library(aricode)

# ------------------------------------------
# PREPARE DATA
# ------------------------------------------
prepare_data <- function(dataset) {
  
  dataset <- as.data.frame(dataset)
  
  if ("class" %in% colnames(dataset)) {
    y <- dataset$class
    X <- dataset[, setdiff(colnames(dataset), "class")]
  } else {
    y <- dataset[, ncol(dataset)]
    X <- dataset[, -ncol(dataset)]
  }
  
  list(X = as.matrix(X), y = y)
}

# ------------------------------------------
# SC_medoids Core Algorithm
# ------------------------------------------
SC_medoids <- function(D, k, E, C) {
  
  # Initial assignment to closest medoid
  cl <- max.col(-D[, C, drop = FALSE])
  
  # Order points by distance to nearest medoid
  sorted_points <- order(apply(D[, C, drop = FALSE], 1, min))
  
  # Enforce cardinality constraints
  for (i in 1:k) {
    cl[sorted_points[1:E[i]]] <- i
    sorted_points <- sorted_points[-(1:E[i])]
  }
  
  # Assign remaining points
  for (point in sorted_points) {
    cl[point] <- which.min(D[point, C])
  }
  
  cl
}

# ------------------------------------------
# MAIN FUNCTION (compatible with Testing2.R)
# ------------------------------------------
run_clustering <- function(dataset, target_cardinality, dataset_name) {
  
  data <- prepare_data(dataset)
  X <- data$X
  y <- data$y
  
  start_total <- proc.time()
  
  # Distance matrix (O(n^2 d))
  D <- as.matrix(proxy::dist(X, method = "euclidean"))
  
  k <- length(target_cardinality)
  
  # Initial medoids using PAM
  pam_result <- pam(X, k)
  C <- pam_result$id.med
  
  # Run SC algorithm
  label_pred <- SC_medoids(D, k, target_cardinality, C)
  
  # Silhouette (uses precomputed D)
  silhouette_values <- silhouette(label_pred, as.dist(D))
  mean_silhouette <- mean(silhouette_values[, "sil_width"])
  
  # External metrics
  ARI_value <- ARI(y, label_pred)
  AMI_value <- AMI(y, label_pred)
  NMI_value <- NMI(y, label_pred)
  
  total_time <- (proc.time() - start_total)[3]
  
  # Cardinalities
  cardinality_pred <- as.integer(table(label_pred))
  violations <- sum(abs(cardinality_pred - target_cardinality))
  
  results_df <- data.frame(
    Violations = violations,
    name = dataset_name,
    ARI = ARI_value,
    AMI = AMI_value,
    NMI = NMI_value,
    Mean_Silhouette = mean_silhouette,
    Clusters = length(unique(label_pred)),
    number_features = ncol(X),
    number_instances = nrow(X),
    Cardinality_pred = paste(cardinality_pred, collapse = " "),
    Cardinality_real = paste(target_cardinality, collapse = " "),
    Execution_Time = total_time,
    stringsAsFactors = FALSE
  )
  
  write.table(results_df,
              file = "results_KMedoidsSC.csv",
              sep = ",",
              row.names = FALSE,
              col.names = !file.exists("results_KMedoidsSC.csv"),
              append = TRUE)
  
  return(total_time)
}