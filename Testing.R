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
run_algorithm <- function(numalt) {
  log_file <- paste0("log_algorithm_", numalt, "_", Sys.getpid(), ".txt")
  
  write(
    paste("Running algorithm", numalt, "in PID:", Sys.getpid(),
          "WD:", getwd(), Sys.time()),
    file = log_file, append = TRUE
  )
  
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
    
    # Use withCallingHandlers so warnings are logged but do NOT interrupt execution
    # tryCatch with warning handler was stopping source() at the first package warning
    withCallingHandlers(
      tryCatch({
        result <- peakRAM(source(algorithm_file, local = FALSE))
        
        result$Algorithm <- algorithm_name
        result$PID <- Sys.getpid()
        result$Date <- Sys.time()
        
        write.table(
          result,
          file = "results_peakRAM.csv",
          sep = ",",
          row.names = FALSE,
          col.names = !file.exists("results_peakRAM.csv"),
          append = TRUE
        )
        
        write(
          paste("Algorithm completed", numalt, "PID:", Sys.getpid(),
                "WD:", getwd(), Sys.time()),
          file = log_file, append = TRUE
        )
      }, error = function(e) {
        write(paste("  ERROR in algorithm", numalt, ":", e$message),
              file = log_file, append = TRUE)
      }),
      warning = function(w) {
        write(paste("  WARNING (non-fatal):", w$message),
              file = log_file, append = TRUE)
        invokeRestart("muffleWarning")  # log it but keep going
      }
    )
    
  } else {
    warning(paste("Invalid option:", numalt))
  }
}

# Main Function
Execute_Test <- function(numalt) {
  script_dir <- getwd()
  
  if (interactive()) {
    active_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    
    if (nzchar(active_path)) {
      script_dir <- dirname(active_path)
    }
  }
  
  if (!dir.exists(script_dir)) {
    stop("The working directory does not exist.")
  }
  
  setwd(script_dir)
  cat("Working directory:", script_dir, "\n")
  
  # Load datasets
  source(file.path(script_dir, "Openml.R"), local = TRUE)
  
  switch(
    as.character(numalt),
    "1" = {
      source(file.path(script_dir, "Bat.R"), local = TRUE)
    },
    "2" = {
      source(file.path(script_dir, "KmedoidsSC.R"), local = TRUE)
    },
    "3" = {
      source(file.path(script_dir, "CSCLP.R"), local = TRUE)
    },
    "4" = {
      source(file.path(script_dir, "MILP-KM.R"), local = TRUE)
    },
    "5" = {
      cat("Running algorithms 1 to 4 in parallel...\n")
      
      options <- 1:4
      num_cores <- min(length(options), max(1, detectCores() - 1))
      cl <- makeCluster(num_cores)
      on.exit(stopCluster(cl), add = TRUE)
      
      clusterExport(
        cl,
        varlist = c("odatasets_unique", "run_algorithm", "script_dir"),
        envir = environment()
      )
      
      clusterEvalQ(cl, {
        if (!requireNamespace("peakRAM", quietly = TRUE)) {
          install.packages("peakRAM", repos = "https://cloud.r-project.org")
        }
        library(peakRAM)
        setwd(script_dir)
        NULL
      })
      
      parLapply(cl, options, run_algorithm)
      
      cat("Parallel execution completed.\n")
    },
    {
      warning("Invalid option.")
    }
  )
}

# Run all
Execute_Test(5)
