library(cluster)
library(proxy)
library(mlr3oml)
library(mlr3)
library(pryr)
library(dplyr) 
library(aricode)
library(ggplot2)
library(corrplot)
library(clValid)
library(RColorBrewer)
library(factoextra) # For PCA visualization

#----------------------------------------------------------------------------------------------#
# KMeansBA
#----------------------------------------------------------------------------------------------#

# Function to prepare data with validation of the dependent variable
prepare_data <- function(dataset) {
  # Ensure that the dataset is a data.frame
  dataset <- as.data.frame(dataset)
  
  # Check if the "class" column exists
  if ("class" %in% colnames(dataset)) {
    # Separate the "class" variable from the rest of the columns
    y <- dataset$class
    X <- dataset[, setdiff(colnames(dataset), "class")]
  } 
  # Validate if the last column is categorical or integer
  else if (is.factor(dataset[, ncol(dataset)]) || is.integer(dataset[, ncol(dataset)])) {
    X <- dataset[, -ncol(dataset)]
    y <- dataset[, ncol(dataset)]
  } 
  # If no condition is met, return NULL to ignore this dataset
  else {
    return(NULL)
  }
  
  # Return a list with variables X and y
  list(X = X, y = y)
}

# Function to calculate the distance matrix
calculate_distance_matrix <- function(X) {
  D <- proxy::dist(as.matrix(X), method = "euclidean")
  as.matrix(D)
}

# Function to adjust cardinality
adjust_cardinality <- function(cluster_assignment, X, centroids, target_cardinality) {
  max_iterations <- 1000  # Avoid infinite loops
  iteration <- 0
  
  # Calculate initial cluster sizes (count vector)
  k <- length(target_cardinality)
  cluster_sizes <- sapply(1:k, function(c) sum(cluster_assignment == c, na.rm = TRUE))
  
  # While any cluster exceeds its cardinality and the iteration limit is not reached
  while(any(cluster_sizes > target_cardinality) && iteration < max_iterations) {
    iteration <- iteration + 1
    
    # Iterate over clusters that exceed the permitted size
    for (j in which(cluster_sizes > target_cardinality)) {
      # Get indices of elements in cluster j
      idx <- which(cluster_assignment == j)
      if (length(idx) == 0) next
      
      # Select the first element to readjust
      element <- idx[1]
      
      # Calculate distances from the element to all centroids (vectorized)
      # Convert X[element,] to numerical vector to avoid formatting issues
      distances <- colSums((t(centroids) - as.numeric(X[element, ]))^2)
      
      # Identify clusters with available space using pre-calculated counts
      available_clusters <- which(cluster_sizes < target_cardinality)
      # Exclude the current cluster (j)
      valid_clusters <- setdiff(available_clusters, j)
      
      if (length(valid_clusters) == 0) {
        # If no cluster has space, mark element as "unassigned"
        cluster_assignment[element] <- -1
        cluster_sizes[j] <- cluster_sizes[j] - 1
      } else {
        # Select the valid cluster with the minimum distance
        chosen <- valid_clusters[which.min(distances[valid_clusters])]
        # Update assignment and cluster sizes
        cluster_assignment[element] <- chosen
        cluster_sizes[j] <- cluster_sizes[j] - 1
        cluster_sizes[chosen] <- cluster_sizes[chosen] + 1
      }
    }
  }
  
  # For unassigned elements (-1), assign them to the cluster with the smallest size (and space)
  unassigned_idx <- which(cluster_assignment == -1)
  if (length(unassigned_idx) > 0) {
    for (element in unassigned_idx) {
      available_clusters <- which(cluster_sizes < target_cardinality)
      if (length(available_clusters) > 0) {
        chosen <- available_clusters[which.min(cluster_sizes[available_clusters])]
        cluster_assignment[element] <- chosen
        cluster_sizes[chosen] <- cluster_sizes[chosen] + 1
      } else {
        # Extreme case: assign arbitrarily to the first cluster
        cluster_assignment[element] <- 1
        cluster_sizes[1] <- cluster_sizes[1] + 1
      }
    }
  }
  
  if (iteration >= max_iterations) {
    warning("The adjustment process reached the maximum number of iterations.")
  }
  
  return(cluster_assignment)
}

# Function to generate an initial solution
generate_initial_solution <- function(X, target_cardinality, seed = 45) {
  set.seed(seed)
  km <- kmeans(X, centers = length(target_cardinality))
  cluster_assignment <- km$cluster
  centroids <- calculate_centroids(X, cluster_assignment, length(target_cardinality))
  adjust_cardinality(cluster_assignment, X, centroids, target_cardinality)
}

# Function to calculate centroids
calculate_centroids <- function(X, cluster_assignment, k) {
  centroids_df <- aggregate(X, by = list(cluster = cluster_assignment), FUN = mean)
  # Ensure the order of centroids
  centroids <- as.matrix(centroids_df[order(centroids_df$cluster), -1, drop = FALSE])
  return(centroids)
}

# Function to evaluate a solution
evaluate_solution <- function(cluster_assignment, X, target_cardinality, penalty_weight = 10) {
  ss <- silhouette(cluster_assignment, dist(X))
  current_counts <- tabulate(cluster_assignment, nbins = length(target_cardinality))
  penalty <- penalty_weight * sum(abs(current_counts - target_cardinality))
  return(mean(ss[, "sil_width"]) - penalty)
}

# Bat algorithm
run_bat_algorithm <- function(X, y, target_cardinality, n_bats = 30, max_iterations = 20,
                              f_min = 0, f_max = 2, loudness = 0.5, pulse_rate = 0.5,
                              alpha = 0.9, gamma = 0.9) {
  set.seed(1521) 
  seeds <- sample(1:10000, n_bats, replace = FALSE)
  bats <- lapply(1:n_bats, function(i) {
    list(
      position = generate_initial_solution(X, target_cardinality, seed = seeds[i]),
      velocity = rep(0, nrow(X)),
      frequency = runif(1, f_min, f_max),
      loudness = loudness,
      pulse_rate = pulse_rate,
      seed = seeds[i]
    )
  })
  best_solution <- bats[[1]]$position
  best_score <- evaluate_solution(best_solution, X, target_cardinality)
  best_seed <- bats[[1]]$seed
  for (iteration in 1:max_iterations) {
    for (i in 1:n_bats) {
      bat <- bats[[i]]
      bat$frequency <- runif(1, f_min, f_max)
      bat$velocity <- bat$velocity + (bat$position - best_solution) * bat$frequency
      new_position <- round(bat$position + bat$velocity)
      new_position <- pmin(pmax(new_position, 1), length(target_cardinality))
      if (runif(1) > bat$pulse_rate) {
        new_position <- sample(1:length(target_cardinality), nrow(X), replace = TRUE)
      }
      for (j in 1:length(target_cardinality)) {
        while (sum(new_position == j) > target_cardinality[j]) {
          idx <- which(new_position == j)
          new_position[sample(idx, 1)] <- sample(1:length(target_cardinality), 1)
        }
      }
      new_score <- evaluate_solution(new_position, X, target_cardinality)
      if (new_score > best_score && runif(1) < bat$loudness) {
        bats[[i]]$position <- new_position
        bats[[i]]$loudness <- alpha * bat$loudness
        bats[[i]]$pulse_rate <- pulse_rate * (1 - exp(-gamma * iteration))
        if (new_score > best_score) {
          best_solution <- new_position
          best_score <- new_score
          best_seed <- bats[[i]]$seed
        }
      }
    }
  }
  list(best_solution = best_solution, best_score = best_score, best_seed = best_seed, seeds = seeds)
}

print_results <- function(results, y, X, D, target_cardinality, dataset_name) {
  best_solution <- results$best_solution
  best_score <- results$best_score
  best_seed <- results$best_seed
  seeds <- results$seeds
  num_instances <- nrow(X)  
  num_variables <- ncol(X) + 1
  
  # Calculate ARI, AMI, and NMI
  ARI_value <- ARI(y, best_solution)
  AMI_value <- AMI(y, best_solution)
  NMI_value <- NMI(y, best_solution)
  
  # Calculate the silhouette coefficient for the final solution
  silhouette_values <- silhouette(x = best_solution, dist = as.dist(D))
  mean_silhouette <- mean(silhouette_values[, "sil_width"])
  
  # Convert silhouette values to a data frame
  silhouette_df <- as.data.frame(silhouette_values[, c("cluster", "sil_width")])
  
  # Save as CSV file
  write.csv(silhouette_df, "silhouette_results.csv", row.names = FALSE)
  
  cat("Silhouette coefficient results have been saved to 'silhouette_results.csv'.\n")
  
  # Count the number of clusters
  num_clusters <- length(unique(best_solution))
  class_dist = as.integer(table(best_solution))
  
  # Store results in the global results data frame
  global_results <<- rbind(global_results, data.frame(
    name = dataset_name,
    Best_Seed = best_seed,
    ARI = ARI_value,
    AMI = AMI_value,
    NMI = NMI_value,
    Mean_Silhouette = mean_silhouette,
    Clusters = num_clusters,
    number_features = num_variables,
    number_instances = num_instances,
    cardinality_BAT = I(list(class_dist)),
    cardinality_REAL = I(list(target_cardinality))
  ))
  
  # Show information in the console
  cat("Average silhouette coefficient:", mean_silhouette, "\n")
  cat("Seeds used for each bat:\n")
  print(seeds)
  cat("\nThe seed of the bat with the best solution was:", best_seed, "\n")
  
  cat("\nOptimal cluster assignment (Bat Algorithm):\n")
  print(table(best_solution))
  cat("Number of clusters:", num_clusters, "\n")
  
  cat("\nAdjusted Rand Index (ARI):", ARI_value, "\n")
  cat("Adjusted Mutual Information (AMI):", AMI_value, "\n")
  cat("Normalized Mutual Information (NMI):", NMI_value, "\n")
}


# Main function to run everything
run_clustering <- function(dataset, target_cardinality, dataset_name) {
  data <- prepare_data(dataset)
  X <- data$X
  y <- data$y
  D <- calculate_distance_matrix(X)
  
  # Measure only BAT algorithm execution:
  start_algo <- Sys.time()
  results <- run_bat_algorithm(X, y, target_cardinality)
  end_algo <- Sys.time()
  algo_time <- as.numeric(difftime(end_algo, start_algo, units = "secs"))
  
  # Print results (includes writing silhouette_results.csv)
  print_results(results, y, X, D, target_cardinality, dataset_name)
  
  # Return algorithm execution time
  return(algo_time)
}

#----------------------------------------------------------------------------------------------#
# Algorithm Execution
#----------------------------------------------------------------------------------------------#

# Create a global data frame to store results
global_results <- data.frame(
  name = character(),
  Best_Seed = integer(),
  ARI = numeric(),
  AMI = numeric(),
  NMI = numeric(),
  Mean_Silhouette = numeric(),
  Clusters = integer(),
  number_features = integer(),
  number_instances = integer(),
  cardinality_BAT = I(list()),
  cardinality_REAL = I(list()),
  stringsAsFactors = FALSE
)

# Vector to store the algorithm execution times (only BAT part)
algorithm_times <- numeric(nrow(odatasets_unique))

for (i in 1:nrow(odatasets_unique)) {
  cat("\n\n--- Executing for dataset at position:", i, "---\n")
  
  tryCatch({
    # Extract dataset, name and target cardinality
    dataset <- odatasets_unique[i]$dataset[[1]]
    dataset_name <- odatasets_unique[i]$name
    prepared_data <- prepare_data(dataset)
    if (is.null(prepared_data)) {
      cat("Dataset at position", i, "is not in a valid format. Skipping.\n")
      algorithm_times[i] <- NA
      next
    }
    
    target_cardinality <- odatasets_unique[i]$class_distribution_vector[[1]]
    if (is.null(target_cardinality)) {
      cat("Target cardinality not available for position", i, ". Skipping.\n")
      algorithm_times[i] <- NA
      next
    }
    
    # Run clustering and measure BAT algorithm execution time
    algo_time <- run_clustering(dataset, target_cardinality, dataset_name)
    algorithm_times[i] <- algo_time
    cat("Algorithm execution time for position", i, ":", algo_time, "seconds\n")
    
  }, error = function(e) {
    cat("Error processing dataset at position", i, ":", e$message, "\n")
    algorithm_times[i] <- NA
  })
}

# Convert lists into columns for visualization
global_results$Cardinality_BAT <- sapply(global_results$cardinality_BAT, paste, collapse = ", ")
global_results$Cardinality_REAL <- sapply(global_results$cardinality_REAL, paste, collapse = ", ")

# Prepare data for violations summary
violations_data <- data.frame(
  Violations = sapply(1:nrow(global_results), function(i) {
    real <- unlist(global_results$cardinality_REAL[i])
    bat <- unlist(global_results$cardinality_BAT[i])
    sum(abs(real - bat))
  })
)
algorithm_times <- algorithm_times[algorithm_times != 0 & !is.na(algorithm_times)]

# Incorporate the execution times into the final data frame
global_results_total <- cbind(violations_data, global_results)
global_results_total$Execution_Time <- algorithm_times

# Optional: Remove cardinality list columns before saving to CSV
global_results_total$cardinality_BAT <- NULL
global_results_total$cardinality_REAL <- NULL

write.csv(global_results_total, "results_KMeansBA.csv", row.names = FALSE)