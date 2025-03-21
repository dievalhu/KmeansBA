# -----------------------------------------------------------------------------
# Cargar bibliotecas necesarias
# -----------------------------------------------------------------------------
if (!require("cluster")) install.packages("cluster")
library(cluster)
if (!require("proxy")) install.packages("proxy")
library(proxy)
if (!require("aricode")) install.packages("aricode")
library(aricode)
if (!require("lpSolve")) install.packages("lpSolve")
library(lpSolve)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
# -----------------------------------------------------------------------------
# Función CSCLP: clustering con restricciones de tamaño mediante programación lineal
# -----------------------------------------------------------------------------
run_CSCLP <- function(dataset, target_cardinality, dataset_name) {
  # Preparar datos
  prepared <- prepare_data(dataset)
  if (is.null(prepared)) stop("El dataset no está en un formato válido.")
  X <- prepared$X
  y <- prepared$y
  
  # Calcular la matriz de distancias usando la distancia coseno
  D_full <- proxy::dist(as.matrix(X), method = "cosine")
  D_full <- as.matrix(D_full)
  distancia <- D_full  # Se usará para la silueta
  
  r <- nrow(D_full)  # Número total de documentos
  k <- length(target_cardinality)  # Número de clusters deseados
  
  # Para la formulación, se usa n2 = (cardinalidad real - 1), pues se excluyen los centroides
  n2 <- matrix(target_cardinality - 1, ncol = 1)
  
  # Obtener índices de centroides mediante PAM (similar a KMedoidsSC)
  pam_result <- pam(X, k)
  c <- pam_result$id.med  # Índices de los centroides
  
  # Obtener el vector de índices de documentos que NO son centroides
  v <- 1:r
  v <- v[!v %in% c]
  
  # Ajustar la matriz de distancias:
  # Eliminar filas de D_full correspondientes a documentos NO centroides
  D <- D_full[-v, ]
  # Convertir las filas en vector (aunque f1 no se utiliza posteriormente)
  f1 <- as.vector(t(D))
  # Eliminar las columnas correspondientes a los centroides
  D <- D[, -c]
  # Convertir la matriz resultante en vector para la programación lineal
  f <- as.vector(t(D))
  
  # Construcción de la matriz de restricciones (Aeq)
  ncolD <- ncol(D)
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
  
  # Configuración del modelo de programación lineal binaria
  A <- Aeq
  dir <- rep("==", nrow(Aeq))
  rhs <- beq
  
  # Medir únicamente el tiempo de resolución del modelo LP y reconstrucción de asignación
  start_alg <- Sys.time()
  result_lp <- lp(direction = "min", objective.in = f, const.mat = A, 
                  const.dir = dir, const.rhs = rhs, all.bin = TRUE)
  end_alg <- Sys.time()
  alg_time <- as.numeric(difftime(end_alg, start_alg, units = "secs"))
  
  x <- result_lp$solution
  fval <- result_lp$objval
  
  s1 <- matrix(x, nrow = k, byrow = TRUE)
  s2 <- diag(k)
  st <- matrix(0, nrow = k, ncol = r)
  st[, c] <- s2
  st[, v] <- s1
  
  label_pred <- numeric(r)
  for (i in 1:k) {
    ind <- which(st[i, ] == 1)
    label_pred[ind] <- i
  }
  
  # Cálculo del coeficiente de silueta utilizando la matriz de distancias original
  silhouette_values <- silhouette(x = label_pred, dist = as.dist(distancia))
  mean_silhouette <- mean(silhouette_values[, "sil_width"])
  
  # Cálculo de las métricas ARI, AMI y NMI
  ARI_value <- ARI(y, label_pred)
  AMI_value <- AMI(y, label_pred)
  NMI_value <- NMI(y, label_pred)
  
  # Obtener la cardinalidad predicha (número de elementos por cluster)
  cardinality_pred <- as.integer(table(label_pred))
  
  # Calcular los incumplimientos en las restricciones de tamaño
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
    Execution_Time = alg_time  # Tiempo medido del algoritmo CSCLP
  ))
}

# -----------------------------------------------------------------------------
# Data frame global para almacenar los resultados de CSCLP
# -----------------------------------------------------------------------------
global_results_CSCLP <- data.frame(
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
# Ejecución del algoritmo CSCLP sobre la lista de datasets (odatasets_unique)
# -----------------------------------------------------------------------------
# Se asume que odatasets_unique está disponible y contiene las siguientes columnas:
# - dataset: lista con los datasets (cada uno accesible como odatasets_unique[i]$dataset[[1]])
# - name: nombre del dataset
# - class_distribution_vector: vector con la cardinalidad real

start_time_total <- Sys.time()

for (i in 1:nrow(odatasets_unique)) {
  cat("\n\n--- Ejecutando CSCLP para dataset en la posición:", i, "---\n")
  start_time <- Sys.time()
  
  tryCatch({
    # Extraer dataset, nombre y cardinalidad real
    dataset <- odatasets_unique[i, ]$dataset[[1]]
    dataset_name <- odatasets_unique[i, ]$name
    target_cardinality <- odatasets_unique[i, ]$class_distribution_vector[[1]]
    
    if (is.null(dataset) || is.null(target_cardinality)) {
      cat("Dataset o cardinalidad no disponibles en la posición", i, ". Saltando...\n")
      next
    }
    
    # Ejecutar CSCLP y obtener el tiempo de ejecución del algoritmo (medido en run_CSCLP)
    result <- run_CSCLP(dataset, target_cardinality, dataset_name)
    
    # Mostrar en consola el incumplimiento de restricciones
    cat("Incumplimiento en las restricciones (Violations):", result$Violations, "\n")
    
    # Agregar resultados al data frame global, incluyendo la columna Execution_Time
    global_results_CSCLP <- rbind(global_results_CSCLP, data.frame(
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
    
    end_time <- Sys.time()
    cat("Tiempo de ejecución para la posición", i, ":", as.numeric(difftime(end_time, start_time, units = "secs")), "segundos\n")
    
  }, error = function(e) {
    cat("Error procesando dataset en la posición", i, ":", e$message, "\n")
  })
}

end_time_total <- Sys.time()
total_execution_time <- as.numeric(difftime(end_time_total, start_time_total, units = "secs"))
cat("\nTiempo total de ejecución:", total_execution_time, "segundos.\n")

# -----------------------------------------------------------------------------
# Preparar los resultados para visualización y guardar en CSV
# -----------------------------------------------------------------------------
global_results_CSCLP$cardinality_pred <- sapply(global_results_CSCLP$cardinality_pred, paste, collapse = ", ")
global_results_CSCLP$cardinality_real <- sapply(global_results_CSCLP$cardinality_real, paste, collapse = ", ")

# (Opcional) Filtrar por nombres comunes si es necesario:
common_names <- intersect(global_results_total$name, global_results_CSCLP$name)
global_results_CSCLP <- global_results_CSCLP[global_results_CSCLP$name %in% common_names, ]

write.csv(global_results_CSCLP, "results_CSCLP.csv", row.names = FALSE)
