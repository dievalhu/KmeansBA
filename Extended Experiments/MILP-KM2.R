# ==========================================
# MILP-KM Optimized Version
# ==========================================

library(lpSolve)
library(aricode)
library(cluster)
library(proxy)
library(mclust)

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
# SOLVE MILP ASSIGNMENT
# ------------------------------------------
solve_milp_assignment <- function(data, centroids, size_constraints) {
  
  n <- nrow(data)
  k <- nrow(centroids)
  
  # Cost matrix O(n k d)
  cost_matrix <- as.matrix(proxy::dist(data, centroids, method = "cosine"))
  cost_vec <- as.vector(t(cost_matrix))
  
  # Constraint 1: each point assigned once
  constr1 <- matrix(0, n, n*k)
  for (i in 1:n) {
    constr1[i, ((i - 1) * k + 1):(i * k)] <- 1
  }
  
  # Constraint 2: cluster size constraints
  constr2 <- matrix(0, k, n*k)
  for (j in 1:k) {
    constr2[j, seq(j, n*k, by = k)] <- 1
  }
  
  f.con <- rbind(constr1, constr2)
  f.dir <- c(rep("=", n), rep("<=", k))
  f.rhs <- c(rep(1, n), size_constraints)
  
  result <- lp("min",
               cost_vec,
               f.con,
               f.dir,
               f.rhs,
               all.bin = TRUE)
  
  if (result$status != 0) {
    stop("MILP did not find a feasible solution")
  }
  
  x_opt <- matrix(result$solution, nrow = n, byrow = TRUE)
  p <- apply(x_opt, 1, which.max)
  
  return(p)
}

# ------------------------------------------
# MAIN ITERATIVE MILP-KM
# ------------------------------------------
clustering_with_size_constraints <- function(data,
                                             size_constraints,
                                             max_iter = 100,
                                             tol = 1e-6) {
  
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  k <- length(size_constraints)
  
  # Random centroid initialization
  set.seed(123)
  random_indices <- sample(1:n, k)
  centroids <- data[random_indices, , drop = FALSE]
  
  converged <- FALSE
  iteration <- 0
  
  while (!converged && iteration < max_iter) {
    
    iteration <- iteration + 1
    
    # Solve MILP
    p <- solve_milp_assignment(data, centroids, size_constraints)
    
    # Update centroids
    new_centroids <- matrix(0, nrow = k, ncol = d)
    for (j in 1:k) {
      cluster_points <- data[p == j, , drop = FALSE]
      if (nrow(cluster_points) > 0) {
        new_centroids[j, ] <- colMeans(cluster_points)
      } else {
        new_centroids[j, ] <- centroids[j, ]
      }
    }
    
    # Convergence check
    if (max(abs(centroids - new_centroids)) < tol) {
      converged <- TRUE
    }
    
    centroids <- new_centroids
  }
  
  list(p = p, centroids = centroids)
}

# ------------------------------------------
# RUN CLUSTERING (Compatible with Testing2.R)
# ------------------------------------------
run_clustering <- function(dataset, target_cardinality, dataset_name) {
  
  data <- prepare_data(dataset)
  X <- data$X
  y <- data$y
  
  start_total <- proc.time()
  
  resultado <- clustering_with_size_constraints(X, target_cardinality)
  
  pred_labels <- resultado$p
  
  # External metrics
  ari <- adjustedRandIndex(y, pred_labels)
  ami <- AMI(y, pred_labels)
  nmi <- NMI(y, pred_labels)
  
  # Silhouette using cosine distance
  D <- proxy::dist(X, method = "cosine")
  sil <- silhouette(pred_labels, D)
  sil_mean <- mean(sil[, 3])
  
  total_time <- (proc.time() - start_total)[3]
  
  num_variables <- ncol(X)
  num_instances <- nrow(X)
  num_clusters <- length(unique(pred_labels))
  
  pred_table <- as.numeric(table(pred_labels))
  true_table <- as.numeric(table(y))
  
  results_df <- data.frame(
    name = dataset_name,
    ARI = ari,
    AMI = ami,
    NMI = nmi,
    Mean_Silhouette = sil_mean,
    Clusters = num_clusters,
    number_features = num_variables,
    number_instances = num_instances,
    cardinality_MILP = paste(pred_table, collapse = " "),
    cardinality_REAL = paste(true_table, collapse = " "),
    Execution_Time = total_time,
    stringsAsFactors = FALSE
  )
  
  write.table(results_df,
              file = "results_MILP.csv",
              sep = ",",
              row.names = FALSE,
              col.names = !file.exists("results_MILP.csv"),
              append = TRUE)
  
  return(total_time)
}