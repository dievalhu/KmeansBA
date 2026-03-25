# ================================
# Testing2.R
# ================================

if (!requireNamespace("peakRAM", quietly = TRUE)) {
  install.packages("peakRAM")
}
library(peakRAM)

# --------------------------------
# Función para cargar dataset CSV
# --------------------------------
load_dataset <- function(file_path) {
  
  data <- read.csv(file_path, stringsAsFactors = FALSE)

  for(i in 1:(ncol(data)-1)){
    data[is.na(data[,i]), i] <- mean(data[,i], na.rm = TRUE)
  }
    
  # Convertir última columna en factor (ground truth)
  data[, ncol(data)] <- as.factor(data[, ncol(data)])
  
  X <- data[, -ncol(data)]
  y <- data[, ncol(data)]
  
  # Cardinalidad real
  target_cardinality <- as.numeric(table(y))
  
  list(
    dataset = data,
    X = X,
    y = y,
    target_cardinality = target_cardinality,
    dataset_name = tools::file_path_sans_ext(basename(file_path))
  )
}

# --------------------------------
# Ejecutar algoritmo
# --------------------------------
Execute_Test <- function(numalt, file_path) {
  
  cat("Loading dataset...\n")
  data_info <- load_dataset(file_path)
  
  dataset <- data_info$dataset
  target_cardinality <- data_info$target_cardinality
  dataset_name <- data_info$dataset_name
  
  cat("Dataset:", dataset_name, "\n")
  cat("Instances:", nrow(dataset), "\n")
  cat("Features:", ncol(dataset)-1, "\n")
  cat("Clusters:", length(target_cardinality), "\n\n")
  
  script_name <- switch(
    as.character(numalt),
    "1" = "KmedoidsSC2.R",
    "2" = "CSCLP2.R",
    "3" = "Bat2.R",
    "4" = "MILP-KM2.R",
    stop("Invalid option")
  )
  
  cat("Running algorithm:", script_name, "\n\n")
  
  # Cargar script del algoritmo
  source(script_name)
  
  # Ejecutar y medir
  result <- peakRAM({
    run_clustering(dataset, target_cardinality, dataset_name)
  })
  
  print(result)
  
  write.table(result,
              file = "resource_usage.csv",
              sep = ",",
              row.names = FALSE,
              col.names = !file.exists("resource_usage.csv"),
              append = TRUE)
}



# ================================
# EJEMPLO DE USO
# ================================
Execute_Test(1, "optdigits.csv")
Execute_Test(2, "optdigits.csv")
Execute_Test(3, "optdigits.csv")
Execute_Test(4, "optdigits.csv")
