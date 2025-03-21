# -----------------------------------------------------------------------------
# Cargar bibliotecas necesarias
# -----------------------------------------------------------------------------
if (!require("cluster")) install.packages("cluster")
library(cluster)
if (!require("proxy")) install.packages("proxy")
library(proxy)
if (!require("aricode")) install.packages("aricode")
library(aricode)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

# -----------------------------------------------------------------------------
# Función para preparar datos con validación de la variable dependiente
# -----------------------------------------------------------------------------
prepare_data <- function(dataset) {
  dataset <- as.data.frame(dataset)
  
  if ("class" %in% colnames(dataset)) {
    y <- dataset$class
    X <- dataset[, setdiff(colnames(dataset), "class")]
  } else if (is.factor(dataset[, ncol(dataset)]) || is.integer(dataset[, ncol(dataset)])) {
    X <- dataset[, -ncol(dataset)]
    y <- dataset[, ncol(dataset)]
  } else {
    return(NULL)
  }
  list(X = X, y = y)
}

# -----------------------------------------------------------------------------
# Función SC_medoids: clustering con KmedoidsSC
# -----------------------------------------------------------------------------
SC_medoids <- function(D, k, E, C = NULL) {
  if (is.null(C)) {
    C <- sample(1:nrow(D), k)
  }
  
  # Asigna inicialmente cada punto al medoid más cercano
  cl <- max.col(-D[, C, drop = FALSE])
  
  # Ordena los puntos por su distancia al medoid más cercano
  sorted_points <- order(apply(D[, C, drop = FALSE], 1, min))
  
  # Asigna los primeros E[i] puntos a cada grupo i
  for (i in 1:k) {
    cl[sorted_points[1:E[i]]] <- i
    sorted_points <- sorted_points[-(1:E[i])]
  }
  
  # Asigna los puntos restantes al medoid más cercano
  for (point in sorted_points) {
    cl[point] <- which.min(D[point, C])
  }
  
  # Genera etiquetas de cluster para cada punto
  labels <- numeric(nrow(D))
  for (i in 1:k) {
    ii <- which(cl == i)
    labels[ii] <- i
  }
  
  return(list(medoids = C, clustering = cl, labels = labels))
}

# -----------------------------------------------------------------------------
# Función para ejecutar KmedoidsSC en un dataset
# -----------------------------------------------------------------------------
run_KmedoidsSC <- function(dataset, target_cardinality, dataset_name) {
  # Preparar datos
  prepared <- prepare_data(dataset)
  if (is.null(prepared)) {
    stop("El dataset no está en un formato válido.")
  }
  X <- prepared$X
  y <- prepared$y
  
  # Calcular la matriz de distancias utilizando la distancia coseno
  D <- proxy::dist(as.matrix(X), method = "cosine")
  D <- as.matrix(D)
  
  # Definir número de clusters a partir de la cardinalidad real
  E <- target_cardinality
  k <- length(E)
  
  # Obtener medoids iniciales con PAM
  pam_result <- pam(X, k)
  C <- pam_result$id.med
  
  # Medir exclusivamente la ejecución de la función SC_medoids
  start_SC <- Sys.time()
  result <- SC_medoids(D, k, E, C)
  end_SC <- Sys.time()
  SC_time <- as.numeric(difftime(end_SC, start_SC, units = "secs"))
  
  label_pred <- result$labels
  
  # Calcular el coeficiente de silueta
  silhouette_values <- silhouette(x = label_pred, dist = as.dist(D))
  mean_silhouette <- mean(silhouette_values[, "sil_width"])
  
  # Calcular ARI, AMI y NMI
  ARI_value <- ARI(y, label_pred)
  AMI_value <- AMI(y, label_pred)
  NMI_value <- NMI(y, label_pred)
  
  # Obtener la cardinalidad predicha (número de elementos por cluster)
  cardinality_pred <- as.integer(table(label_pred))
  
  # Calcular incumplimientos en las restricciones de tamaño
  violations <- sum(abs(cardinality_pred - target_cardinality))
  
  # Número de clusters, instancias y características
  num_clusters <- length(unique(label_pred))
  num_instances <- nrow(X)
  num_features <- ncol(X) + 1
  
  return(list(
    dataset_name = dataset_name,
    ARI = ARI_value,
    AMI = AMI_value,
    NMI = NMI_value,
    Mean_Silhouette = mean_silhouette,
    Clusters = num_clusters,
    number_features = num_features,
    number_instances = num_instances,
    cardinality_pred = cardinality_pred,
    cardinality_real = target_cardinality,
    Violations = violations,
    Execution_Time = SC_time  # Nombre de columna unificado
  ))
}

# -----------------------------------------------------------------------------
# Data frame global para almacenar los resultados de KmedoidsSC
# -----------------------------------------------------------------------------
global_results_KmedoidsSC <- data.frame(
  name = character(),
  ARI = numeric(),
  AMI = numeric(),
  NMI = numeric(),
  Mean_Silhouette = numeric(),
  Clusters = integer(),
  number_features = integer(),
  number_instances = integer(),
  cardinality_pred = I(list()),
  cardinality_real = I(list()),
  Violations = numeric(),
  Execution_Time = numeric(),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Ejecución del algoritmo KmedoidsSC sobre la lista de datasets (odatasets_unique)
# -----------------------------------------------------------------------------
# Se asume que odatasets_unique está disponible y tiene las siguientes columnas:
# - dataset: lista con los datasets (cada elemento es accesible como odatasets_unique[i]$dataset[[1]])
# - name: nombre del dataset
# - class_distribution_vector: vector con la cardinalidad real

start_time_total <- Sys.time()

for (i in 1:nrow(odatasets_unique)) {
  cat("\n\n--- Ejecutando para dataset en la posición:", i, "---\n")
  
  tryCatch({
    # Extraer dataset, nombre y cardinalidad real
    dataset <- odatasets_unique[i, ]$dataset[[1]]
    dataset_name <- odatasets_unique[i, ]$name
    target_cardinality <- odatasets_unique[i, ]$class_distribution_vector[[1]]
    
    if (is.null(dataset) || is.null(target_cardinality)) {
      cat("Dataset o cardinalidad no disponibles en la posición", i, ". Saltando...\n")
      next
    }
    
    # Ejecutar KmedoidsSC y obtener el tiempo de ejecución de SC_medoids
    result <- run_KmedoidsSC(dataset, target_cardinality, dataset_name)
    
    # Mostrar en consola el incumplimiento de restricciones
    cat("Incumplimiento en las restricciones de tamaño (Violations):", result$Violations, "\n")
    
    # Agregar resultados al data frame global, incluyendo la columna Execution_Time
    global_results_KmedoidsSC <- rbind(global_results_KmedoidsSC, data.frame(
      name = result$dataset_name,
      ARI = result$ARI,
      AMI = result$AMI,
      NMI = result$NMI,
      Mean_Silhouette = result$Mean_Silhouette,
      Clusters = result$Clusters,
      number_features = result$number_features,
      number_instances = result$number_instances,
      cardinality_pred = I(list(result$cardinality_pred)),
      cardinality_real = I(list(result$cardinality_real)),
      Violations = result$Violations,
      Execution_Time = result$Execution_Time,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    cat("Error procesando dataset en la posición", i, ":", e$message, "\n")
  })
}

end_time_total <- Sys.time()
total_execution_time <- as.numeric(difftime(end_time_total, start_time_total, units = "secs"))
cat("\nTiempo total de ejecución (incluyendo preparación e iteración sobre datasets):", total_execution_time, "segundos.\n")

# -----------------------------------------------------------------------------
# Preparar los resultados para visualización y guardar en CSV
# -----------------------------------------------------------------------------
global_results_KmedoidsSC$cardinality_pred <- sapply(global_results_KmedoidsSC$cardinality_pred, paste, collapse = ", ")
global_results_KmedoidsSC$cardinality_real <- sapply(global_results_KmedoidsSC$cardinality_real, paste, collapse = ", ")

# Obtener la intersección de nombres (datasets comunes)
common_names <- intersect(global_results_total$name, global_results_KmedoidsSC$name)

# Filtrar global_results_KmedoidsSC para que sólo contenga los datasets con nombres comunes
global_results_KmedoidsSC <- global_results_KmedoidsSC[global_results_KmedoidsSC$name %in% common_names, ]

# Escribir los resultados en un archivo CSV
write.csv(global_results_KmedoidsSC, "results_KmedoidsSC.csv", row.names = FALSE)
