# ==========================================
# KMeansBA Optimized Version
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
# DISTANCE MATRIX (computed once)
# ------------------------------------------
calculate_distance_matrix <- function(X) {
  as.matrix(proxy::dist(X, method = "cosine"))
}

# ------------------------------------------
# ADJUST CARDINALITY
# ------------------------------------------
adjust_cardinality <- function(cluster_assignment, X, centroids, target_cardinality) {
  
  k <- length(target_cardinality)
  cluster_sizes <- tabulate(cluster_assignment, nbins = k)
  
  for (j in which(cluster_sizes > target_cardinality)) {
    while (cluster_sizes[j] > target_cardinality[j]) {
      
      idx <- which(cluster_assignment == j)
      element <- idx[1]
      
      distances <- colSums((t(centroids) - X[element, ])^2)
      available_clusters <- which(cluster_sizes < target_cardinality)
      valid_clusters <- setdiff(available_clusters, j)
      
      if (length(valid_clusters) == 0) break
      
      chosen <- valid_clusters[which.min(distances[valid_clusters])]
      cluster_assignment[element] <- chosen
      
      cluster_sizes[j] <- cluster_sizes[j] - 1
      cluster_sizes[chosen] <- cluster_sizes[chosen] + 1
    }
  }
  
  cluster_assignment
}

# ------------------------------------------
# CENTROIDS
# ------------------------------------------
calculate_centroids <- function(X, cluster_assignment, k) {
  centroids <- matrix(0, nrow = k, ncol = ncol(X))
  
  for (j in 1:k) {
    centroids[j, ] <- colMeans(X[cluster_assignment == j, , drop = FALSE])
  }
  
  centroids
}

# ------------------------------------------
# INITIAL SOLUTION
# ------------------------------------------
generate_initial_solution <- function(X, target_cardinality, seed) {
  set.seed(seed)
  km <- kmeans(X, centers = length(target_cardinality), nstart = 1)
  cluster_assignment <- km$cluster
  centroids <- calculate_centroids(X, cluster_assignment, length(target_cardinality))
  adjust_cardinality(cluster_assignment, X, centroids, target_cardinality)
}

# ------------------------------------------
# EVALUATION (uses precomputed D)
# ------------------------------------------
evaluate_solution <- function(cluster_assignment, D, target_cardinality, penalty_weight = 10) {
  
  ss <- silhouette(cluster_assignment, as.dist(D))
  mean_sil <- mean(ss[, "sil_width"])
  
  current_counts <- tabulate(cluster_assignment, nbins = length(target_cardinality))
  penalty <- penalty_weight * sum(abs(current_counts - target_cardinality))
  
  mean_sil - penalty
}

# ------------------------------------------
# BAT ALGORITHM (optimized)
# ------------------------------------------
run_bat_algorithm <- function(X, D, target_cardinality,
                              n_bats = 30,
                              max_iterations = 20,
                              f_min = 0,
                              f_max = 2,
                              loudness = 0.5,
                              pulse_rate = 0.5,
                              alpha = 0.9,
                              gamma = 0.9) {
  
  set.seed(1521)
  seeds <- sample(1:10000, n_bats)
  k <- length(target_cardinality)
  n <- nrow(X)
  
  bats <- lapply(1:n_bats, function(i) {
    list(
      position = generate_initial_solution(X, target_cardinality, seeds[i]),
      velocity = rep(0, n),
      seed = seeds[i],
      loudness = loudness,
      pulse_rate = pulse_rate
    )
  })
  
  best_solution <- bats[[1]]$position
  best_score <- evaluate_solution(best_solution, D, target_cardinality)
  best_seed <- bats[[1]]$seed
  
  for (iteration in 1:max_iterations) {
    for (i in 1:n_bats) {
      
      bat <- bats[[i]]
      freq <- runif(1, f_min, f_max)
      
      bat$velocity <- bat$velocity + (bat$position - best_solution) * freq
      new_position <- round(bat$position + bat$velocity)
      new_position <- pmin(pmax(new_position, 1), k)
      
      if (runif(1) > bat$pulse_rate) {
        new_position <- sample(1:k, n, replace = TRUE)
      }
      
      # Hard constraint
      counts <- tabulate(new_position, nbins = k)
      for (j in which(counts > target_cardinality)) {
        while (counts[j] > target_cardinality[j]) {
          idx <- which(new_position == j)
          chosen <- sample(idx, 1)
          new_cluster <- sample(setdiff(1:k, j), 1)
          new_position[chosen] <- new_cluster
          counts <- tabulate(new_position, nbins = k)
        }
      }
      
      new_score <- evaluate_solution(new_position, D, target_cardinality)
      
      if (new_score > best_score && runif(1) < bat$loudness) {
        bats[[i]]$position <- new_position
        best_solution <- new_position
        best_score <- new_score
        best_seed <- bat$seed
      }
    }
  }
  
  list(best_solution = best_solution,
       best_seed = best_seed)
}

# ------------------------------------------
# MAIN FUNCTION
# ------------------------------------------
run_clustering <- function(dataset, target_cardinality, dataset_name) {
  
  data <- prepare_data(dataset)
  X <- data$X
  y <- data$y
  
  start_total <- proc.time()
  
  # Distance matrix once
  D <- calculate_distance_matrix(X)
  
  results <- run_bat_algorithm(X, D, target_cardinality)
  
  best_solution <- results$best_solution
  best_seed <- results$best_seed
  
  ARI_value <- ARI(y, best_solution)
  AMI_value <- AMI(y, best_solution)
  NMI_value <- NMI(y, best_solution)
  
  silhouette_values <- silhouette(best_solution, as.dist(D))
  mean_silhouette <- mean(silhouette_values[, "sil_width"])
  
  total_time <- (proc.time() - start_total)[3]
  
  cardinality_pred <- as.integer(table(best_solution))
  violations <- sum(abs(cardinality_pred - target_cardinality))
  
  results_df <- data.frame(
    Violations = violations,
    name = dataset_name,
    Best_Seed = best_seed,
    ARI = ARI_value,
    AMI = AMI_value,
    NMI = NMI_value,
    Mean_Silhouette = mean_silhouette,
    Clusters = length(unique(best_solution)),
    number_features = ncol(X),
    number_instances = nrow(X),
    Cardinality_BAT = paste(cardinality_pred, collapse = " "),
    Cardinality_REAL = paste(target_cardinality, collapse = " "),
    Execution_Time = total_time,
    stringsAsFactors = FALSE
  )
  
  write.table(results_df,
              file = "results_KMeansBA.csv",
              sep = ",",
              row.names = FALSE,
              col.names = !file.exists("results_KMeansBA.csv"),
              append = TRUE)
  
  return(total_time)
}