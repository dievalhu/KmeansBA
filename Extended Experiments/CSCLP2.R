# ==========================================
# CSCLP Optimized Version
# ==========================================

library(cluster)
library(proxy)
library(aricode)
library(lpSolve)

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
# MAIN FUNCTION (compatible with Testing2.R)
# ------------------------------------------
run_clustering <- function(dataset, target_cardinality, dataset_name) {
  
  data <- prepare_data(dataset)
  X <- data$X
  y <- data$y
  
  start_total <- proc.time()
  
  # ------------------------------------------
  # Distance matrix (O(n^2 d))
  # ------------------------------------------
  D_full <- as.matrix(proxy::dist(X, method = "euclidean"))
  
  r <- nrow(D_full)
  k <- length(target_cardinality)
  
  # LP formulation uses n2 = cardinality - 1
  n2 <- matrix(target_cardinality - 1, ncol = 1)
  
  # ------------------------------------------
  # Medoid initialization (PAM)
  # ------------------------------------------
  pam_result <- pam(X, k)
  c_idx <- pam_result$id.med
  
  # Non-centroid indices
  v <- setdiff(1:r, c_idx)
  
  # ------------------------------------------
  # Reduced cost matrix
  # ------------------------------------------
  D <- D_full[-v, ]
  D <- D[, -c_idx]
  f <- as.vector(t(D))
  
  ncolD <- ncol(D)
  
  # ------------------------------------------
  # Build constraint matrix
  # ------------------------------------------
  Aeq1 <- matrix(0, nrow = k, ncol = length(f))
  for (i in 1:ncolD) {
    for (j in 1:k) {
      Aeq1[j, i + (j - 1) * ncolD] <- 1
    }
  }
  
  Aeq2 <- diag(ncolD)
  if (k > 1) {
    for (i in 1:(k - 1)) {
      Aeq2 <- cbind(Aeq2, diag(ncolD))
    }
  }
  
  Aeq <- rbind(Aeq1, Aeq2)
  beq <- c(n2, rep(1, ncolD))
  
  # ------------------------------------------
  # Solve Binary LP
  # ------------------------------------------
  result_lp <- lp(direction = "min",
                  objective.in = f,
                  const.mat = Aeq,
                  const.dir = rep("==", nrow(Aeq)),
                  const.rhs = beq,
                  all.bin = TRUE)
  
  x <- result_lp$solution
  
  # ------------------------------------------
  # Reconstruct clustering
  # ------------------------------------------
  s1 <- matrix(x, nrow = k, byrow = TRUE)
  s2 <- diag(k)
  st <- matrix(0, nrow = k, ncol = r)
  st[, c_idx] <- s2
  st[, v] <- s1
  
  label_pred <- numeric(r)
  for (i in 1:k) {
    ind <- which(st[i, ] == 1)
    label_pred[ind] <- i
  }
  
  # ------------------------------------------
  # Metrics
  # ------------------------------------------
  silhouette_values <- silhouette(label_pred, as.dist(D_full))
  mean_silhouette <- mean(silhouette_values[, "sil_width"])
  
  ARI_value <- ARI(y, label_pred)
  AMI_value <- AMI(y, label_pred)
  NMI_value <- NMI(y, label_pred)
  
  total_time <- (proc.time() - start_total)[3]
  
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
              file = "results_CSCLP.csv",
              sep = ",",
              row.names = FALSE,
              col.names = !file.exists("results_CSCLP.csv"),
              append = TRUE)
  
  return(total_time)
}