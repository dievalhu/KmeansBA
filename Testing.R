# Install and load necessary packages
if (!requireNamespace("rstudioapi", quietly = TRUE)) {
  install.packages("rstudioapi")
}
if (!requireNamespace("peakRAM", quietly = TRUE)) {
  install.packages("peakRAM")
}

library(rstudioapi)
library(parallel)
library(peakRAM)

# Function to run an algorithm and record its resource usage
run_algoritmo <- function(numalt) {
  log_file <- paste0("log_algoritmo_", numalt, "_", Sys.getpid(), ".txt")
  write(paste("Running algorithm", numalt, "in PID:", Sys.getpid(), Sys.time()),
        file = log_file, append = TRUE)
  
  # Names of algorithms and associated scripts
  algorithm_name <- switch(
    as.character(numalt),
    "1" = "K-MeansBA",
    "2" = "K-MedoidsSC",
    "3" = "CSCLP",    
    "4" = "MILP-KM",
    NA
  )
  
  algorithm_file <- switch(
    as.character(numalt),
    "1" = "Bat.R",
    "2" = "KmedoidsSC.R",
    "3" = "CSCLP.R",
    "4" = "MILP-KM.R",
    NA
  )
  
  if (!is.na(algorithm_file)) {
    # Measure time and RAM with peakRAM
    result <- peakRAM(source(algorithm_file, local = TRUE))
    
    # Add extra info
    result$Algoritmo <- algorithm_name
    result$PID <- Sys.getpid()
    result$Fecha <- Sys.time()
    
    # Save in CSV
    write.table(result, file = "results_peakRAM.csv",
                sep = ",", row.names = FALSE,
                col.names = !file.exists("results_peakRAM.csv"),
                append = TRUE)
    
    write(paste("Algorithm completed", numalt, "PID:", Sys.getpid(), Sys.time()),
          file = log_file, append = TRUE)
    
  } else {
    warning(paste("Invalid option:", numalt))
  }
}

# Main Function
Execute_Test <- function(numalt) {
  if (interactive()) {
    script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
    print(script_dir)
    if (file.exists(script_dir)) {
      setwd(script_dir)
    } else {
      cat("The directory does not exist..\n")
    }
  }
  
  # Load datasets
  source("Openml.R")
  
  # Algorithm selection or parallel execution
  switch(
    as.character(numalt),
    "1" = { source("Bat.R") },
    "2" = { source("KmedoidsSC.R") },
    "3" = { source("CSCLP.R") },
    "4" = { source("MILP-KM.R")},
    "5" = {
      cat("Running algorithms 1 to 4 in parallel...\n")
      opciones <- 1:4
      num_cores <- min(length(opciones), detectCores() - 1)
      cl <- makeCluster(num_cores)
      
      # 1. Load 'peakRAM' on all nodes in the cluster
      clusterEvalQ(cl, {
        if (!requireNamespace("peakRAM", quietly = TRUE)) {
          install.packages("peakRAM", repos = "https://cloud.r-project.org")
        }
        library(peakRAM)
      })
      
      # 2. Export the necessary functions and variables
      clusterExport(cl, varlist = c("odatasets_unique", "run_algoritmo"), envir = environment())
      
      # 3. Run in parallel
      parLapply(cl, opciones, run_algoritmo)
      
      stopCluster(cl)
      cat("Parallel execution completed.\n")
    },
    {
      warning("Invalid option.")
    }
  )
  
}

# Run all
Execute_Test(4)
