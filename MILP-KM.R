# Install packages if needed
packages <- c("flexclust", "lpSolve", "mclust", "aricode", "cluster", "proxy")
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install) > 0) install.packages(to_install)

library(flexclust)
library(lpSolve)
library(mclust)     # ARI
library(aricode)    # AMI, NMI
library(cluster)    # Silhouette
library(proxy)

prepare_data <- function(dataset) {
  dataset <- as.data.frame(dataset)
  
  if ("class" %in% colnames(dataset)) {
    y <- dataset$class
    X <- dataset[, setdiff(colnames(dataset), "class"), drop = FALSE]
  } else if (is.factor(dataset[, ncol(dataset)]) || is.integer(dataset[, ncol(dataset)])) {
    X <- dataset[, -ncol(dataset), drop = FALSE]
    y <- dataset[, ncol(dataset)]
  } else {
    return(NULL)
  }
  
  list(X = X, y = y)
}

global_results_MILP <- data.frame(
  name = character(),
  ARI = numeric(),
  AMI = numeric(),
  NMI = numeric(),
  Mean_Silhouette = numeric(),
  Clusters = integer(),
  number_features = integer(),
  number_instances = integer(),
  cardinality_MILP = character(),
  cardinality_REAL = character(),
  Execution_Time = numeric(),
  stringsAsFactors = FALSE
)

#----------------------------------------------------------------------------------------------#
# MILP-KM
#----------------------------------------------------------------------------------------------#
run_MILP <- function(X, y, target_cardinality) {
  list(
    best_solution = best_cluster_assignment,
    mean_silhouette = mean_silhouette,
    AMI = AMI,
    ARI = ARI,
    NMI = NMI,
    Cluster = length(unique(best_cluster_assignment)),
    cardinality_real = target_cardinality,
    cardinality_pred = best_cluster_assignment
  )
}

# ---------- MILP Solver ----------
solve_milp_assignment <- function(data, centroids, size_constraints) {
  n <- nrow(data)
  k <- nrow(centroids)
  
  cost_matrix <- proxy::dist(data, centroids, method = "cosine")
  cost_matrix <- as.matrix(cost_matrix)
  cost_vec <- as.vector(t(cost_matrix))
  
  f.obj <- cost_vec
  
  constr1 <- matrix(0, n, n * k)
  for (i in 1:n) {
    constr1[i, ((i - 1) * k + 1):(i * k)] <- 1
  }
  
  constr2 <- matrix(0, k, n * k)
  for (j in 1:k) {
    constr2[j, seq(j, n * k, by = k)] <- 1
  }
  
  f.con <- rbind(constr1, constr2)
  f.dir <- c(rep("=", n), rep("<=", k))
  f.rhs <- c(rep(1, n), size_constraints)
  
  result <- lp("min", f.obj, f.con, f.dir, f.rhs, all.bin = TRUE)
  
  if (result$status != 0) {
    stop("MILP found no solution")
  }
  
  x_opt <- matrix(result$solution, nrow = n, byrow = TRUE)
  p <- apply(x_opt, 1, which.max)
  
  list(p = p, q = NULL)
}

# ---------- Main Function ----------
clustering_with_size_constraints <- function(data, size_constraints, max_iter = 100, tol = 1e-6) {
  data <- as.matrix(data)
  colnames(data) <- NULL
  d <- ncol(data)
  k <- length(size_constraints)
  
  random_indices <- sample(1:nrow(data), k)
  centroids <- data[random_indices, , drop = FALSE]
  converged <- FALSE
  iteration <- 0
  
  while (!converged && iteration < max_iter) {
    iteration <- iteration + 1
    
    milp_result <- solve_milp_assignment(data, centroids, size_constraints)
    p <- milp_result$p
    
    new_centroids <- matrix(0, nrow = k, ncol = d)
    for (j in 1:k) {
      cluster_points <- data[p == j, , drop = FALSE]
      if (nrow(cluster_points) > 0) {
        new_centroids[j, ] <- colMeans(cluster_points)
      } else {
        new_centroids[j, ] <- centroids[j, ]
      }
    }
    
    if (max(abs(centroids - new_centroids)) < tol) {
      converged <- TRUE
    }
    
    centroids <- new_centroids
  }
  
  list(p = p, centroids = centroids)
}

run_clustering <- function(dataset, target_cardinality, dataset_name) {
  data <- prepare_data(dataset)
  if (is.null(data)) {
    stop("Dataset is not in a valid format.")
  }
  
  X <- data$X
  y <- data$y
  
  start_algo <- Sys.time()
  result <- clustering_with_size_constraints(X, target_cardinality)
  end_algo <- Sys.time()
  
  pred_labels <- result$p
  ari <- adjustedRandIndex(y, pred_labels)
  ami <- AMI(y, pred_labels)
  nmi <- NMI(y, pred_labels)
  
  d <- proxy::dist(X, method = "cosine")
  sil <- silhouette(pred_labels, d)
  sil_mean <- mean(sil[, 3])
  
  algo_time <- as.numeric(difftime(end_algo, start_algo, units = "secs"))
  num_variables <- ncol(X)
  num_instances <- nrow(X)
  num_clusters <- length(unique(pred_labels))
  
  cat("Evaluation metrics:\n")
  cat(sprintf("ARI: %.4f\n", ari))
  cat(sprintf("AMI: %.4f\n", ami))
  cat(sprintf("NMI: %.4f\n", nmi))
  cat(sprintf("Silhouette Index: %.4f\n", sil_mean))
  
  pred_labels_table <- table(pred_labels)
  true_labels_table <- table(y)
  
  new_row <- data.frame(
    name = dataset_name,
    ARI = ari,
    AMI = ami,
    NMI = nmi,
    Mean_Silhouette = sil_mean,
    Clusters = num_clusters,
    number_features = num_variables,
    number_instances = num_instances,
    cardinality_MILP = paste(as.numeric(pred_labels_table), collapse = " "),
    cardinality_REAL = paste(as.numeric(true_labels_table), collapse = " "),
    Execution_Time = algo_time,
    stringsAsFactors = FALSE
  )
  
  global_results_MILP <<- rbind(global_results_MILP, new_row)
  
  out_file <- file.path(script_dir, "results_MILP.csv")
  write.table(
    new_row,
    file = out_file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(out_file),
    append = TRUE
  )
  
  return(algo_time)
}

algorithm_times <- numeric(nrow(odatasets_unique))

for (i in 1:nrow(odatasets_unique)) {
  cat("\n\n--- Executing for dataset at position:", i, "---\n")
  
  tryCatch({
    dataset <- odatasets_unique[i, ]$dataset[[1]]
    dataset_name <- odatasets_unique[i, ]$name
    prepared_data <- prepare_data(dataset)
    
    if (is.null(prepared_data)) {
      cat("Dataset at position", i, "is not in a valid format. Skipping.\n")
      algorithm_times[i] <- NA
      next
    }
    
    target_cardinality <- odatasets_unique[i, ]$class_distribution_vector[[1]]
    if (is.null(target_cardinality)) {
      cat("Target cardinality not available for position", i, ". Skipping.\n")
      algorithm_times[i] <- NA
      next
    }
    
    algo_time <- run_clustering(dataset, target_cardinality, dataset_name)
    algorithm_times[i] <- algo_time
    cat("Algorithm execution time for position", i, ":", algo_time, "seconds\n")
    
  }, error = function(e) {
    cat("Error processing dataset at position", i, ":", e$message, "\n")
    algorithm_times[i] <- NA
  })
}

if (nrow(global_results_MILP) > 0) {
  write.csv(global_results_MILP, file.path(script_dir, "results_MILP.csv"), row.names = FALSE)
}

cat("\nResults saved to results_MILP.csv\n")
