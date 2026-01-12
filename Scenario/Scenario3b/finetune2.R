mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Scenario/Scenario3b/TriangleLassoHD2.R'))
source(paste0(mywd,'/Scenario/Scenario3b/MarsParallelHD.R'))
source(paste0(mywd,'/Functions Dai/knockoff.R'))
source(paste0(mywd,'/Functions Dai/analysis.R'))
source(paste0(mywd,'/Functions Dai/MBHq.R'))
source(paste0(mywd,'/Functions Dai/DS.R'))
source(paste0(mywd,'/Functions Dai/fdp_power.R'))

delta <- 10

evaluate_mars_config <- function(X, y, sample_index1, sample_index2, 
                                 signal_index, q,
                                 degree = 1, nprune = NULL, minspan = 2, 
                                 thresh = 0, nfold = 5,
                                 permR2Mars_func, SelectFeatures_func, 
                                 CalculateFDP_Power_func) {

  p <- ncol(X)
  dataTrain <- data.frame(y = y, X)
  
  # Fit MARS model with specified hyperparameters
  mars_model <- tryCatch({
    earth(
      y ~ .,
      pmethod = "cv",
      nfold = nfold,
      degree = degree,
      nprune = nprune,
      minspan = minspan,
      thresh = thresh,
      data = dataTrain
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(mars_model)) {
    return(list(fdp = NA, power = NA, R2_1 = NA, R2_2 = NA, 
                n_selected = NA, converged = FALSE))
  }
  
  # --- R^2 on train/test halves ---
  pred1 <- predict(mars_model, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
  pred2 <- predict(mars_model, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
  y1 <- y[sample_index1]
  y2 <- y[sample_index2]
  
  R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
  R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
  
  # --- Permutation-based drops ---
  Rnew1 <- sapply(seq_len(p), function(j)
    permR2Mars_func(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = mars_model))
  Rnew2 <- sapply(seq_len(p), function(j)
    permR2Mars_func(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = mars_model))
  
  beta1 <- R2orig1 - Rnew1
  beta2 <- R2orig2 - Rnew2
  mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
  
  # Feature selection
  selected_index <- SelectFeatures_func(mirror, abs(mirror), q)
  res <- CalculateFDP_Power_func(selected_index, signal_index)
  
  return(list(
    fdp = res$fdp,
    power = res$power,
    R2_1 = R2orig1,
    R2_2 = R2orig2,
    n_selected = length(selected_index),
    n_terms = length(mars_model$selected.terms),
    converged = TRUE,
    mirror_stats = list(
      mean = mean(mirror[-signal_index]),
      sd = sd(mirror[-signal_index])
    )
  ))
}



### algorithmic settings
num_split <- 5
n <-400
p <- 500
p0 <- 25
q <- 0.1
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)

n1 <- floor(n/2); n2 <- n - n1
X1 <- matrix(rnorm(n1*p, mean= -1), n1, p)
X2 <- matrix(rnorm(n2*p, mean= 1), n2, p)
X  <- rbind(X1, X2)
beta_star <- numeric(p)
beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))
y <- (X %*% beta_star + rnorm(n))
amountTrain=0.333
amountTest=1-amountTrain
data<-data.frame(cbind(y,X))
n <- dim(X)[1]
p <- dim(X)[2]
# --- indices ---
train_index <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
remaining_index <- setdiff(seq_len(n), train_index)
data<-data.frame(cbind(y,X))
remaining_index<-c(setdiff(c(1:n),train_index))
sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
sample_index2 <- setdiff(remaining_index, sample_index1)



tune_mars_two_stage <- function(
    X, y, sample_index1, sample_index2, signal_index, q,
    permR2Mars_func, SelectFeatures_func, CalculateFDP_Power_func,
    fdr_threshold = 0.10,
    verbose = TRUE) {
  
  if (verbose) cat("=== Stage 1: Coarse Grid Search ===\n")
  
  # Stage 1: Coarse grid
  stage1_results <- tune_mars_hyperparameters(
    X = X, y = y,
    sample_index1 = sample_index1,
    sample_index2 = sample_index2,
    signal_index = signal_index,
    q = q,
    permR2Mars_func = permR2Mars_func,
    SelectFeatures_func = SelectFeatures_func,
    CalculateFDP_Power_func = CalculateFDP_Power_func,
    fdr_threshold = fdr_threshold,
    degree_grid = 1,
    nprune_grid = c(NULL, 10,20,30, 40),
    minspan_grid = c(1,2, 3,4,5),
    thresh_grid = c(0, 0.01,0.025,0.1),
    verbose = verbose
  )
  
  if (is.null(stage1_results$best_params)) {
    return(stage1_results)
  }
  
  best1 <- stage1_results$best_params
  
  if (verbose) cat("\n=== Stage 2: Fine Grid Search ===\n")
  
  # Stage 2: Fine grid around best from stage 1
  degree_fine <- unique(c(max(1, best1$degree - 1), best1$degree, best1$degree + 1))
  degree_fine <- degree_fine[degree_fine >= 1 & degree_fine <= 3]
  
  minspan_fine <- unique(c(max(1, best1$minspan - 1), best1$minspan, best1$minspan + 1))
  
  nprune_base <- if(is.na(best1$nprune)) 30 else best1$nprune
  nprune_fine <- unique(c(NULL, nprune_base - 10, nprune_base, nprune_base + 10))
  nprune_fine <- nprune_fine[is.null(nprune_fine) | (nprune_fine > 0)]
  
  stage2_results <- tune_mars_hyperparameters(
    X = X, y = y,
    sample_index1 = sample_index1,
    sample_index2 = sample_index2,
    signal_index = signal_index,
    q = q,
    permR2Mars_func = permR2Mars_func,
    SelectFeatures_func = SelectFeatures_func,
    CalculateFDP_Power_func = CalculateFDP_Power_func,
    fdr_threshold = fdr_threshold,
    degree_grid = degree_fine,
    nprune_grid = nprune_fine,
    minspan_grid = minspan_fine,
    thresh_grid = c(0, 0.001, 0.005, 0.01),
    verbose = verbose
  )
  
  return(stage2_results)
}
# =============================================================================
# MARS Hyperparameter Tuning for Power Maximization with FDR Control ?earth
# =============================================================================

library(earth)

# =============================================================================
# Helper Functions
# =============================================================================

# Function to evaluate a single hyperparameter configuration
evaluate_mars_config <- function(X, y, sample_index1, sample_index2, 
                                 signal_index, q,
                                 degree = 1, nprune = NULL, minspan = 2, 
                                 thresh = 0, nfold = 5,
                                 permR2Mars_func, SelectFeatures_func, 
                                 CalculateFDP_Power_func) {
  
  p <- ncol(X)
  dataTrain <- data.frame(y = y, X)
  
  # Fit MARS model with specified hyperparameters
  mars_model <- tryCatch({
    earth(
      y ~ .,
      pmethod = "cv",
      nfold = nfold,
      degree = degree,
      nprune = nprune,
      minspan = minspan,
      thresh = thresh,
      data = dataTrain
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(mars_model)) {
    return(list(fdp = NA, power = NA, R2_1 = NA, R2_2 = NA, 
                n_selected = NA, converged = FALSE))
  }
  
  # --- R^2 on train/test halves ---
  pred1 <- predict(mars_model, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
  pred2 <- predict(mars_model, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
  y1 <- y[sample_index1]
  y2 <- y[sample_index2]
  
  R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
  R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
  
  # --- Permutation-based drops ---
  Rnew1 <- sapply(seq_len(p), function(j)
    permR2Mars_func(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = mars_model))
  Rnew2 <- sapply(seq_len(p), function(j)
    permR2Mars_func(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = mars_model))
  
  beta1 <- R2orig1 - Rnew1
  beta2 <- R2orig2 - Rnew2
  mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
  
  # Feature selection
  selected_index <- SelectFeatures_func(mirror, abs(mirror), q)
  res <- CalculateFDP_Power_func(selected_index, signal_index)
  
  return(list(
    fdp = res$fdp,
    power = res$power,
    R2_1 = R2orig1,
    R2_2 = R2orig2,
    n_selected = length(selected_index),
    n_terms = length(mars_model$selected.terms),
    converged = TRUE,
    mirror_stats = list(
      mean = mean(mirror[-signal_index]),
      sd = sd(mirror[-signal_index])
    )
  ))
}

# =============================================================================
# Main Hyperparameter Tuning Function
# =============================================================================

tune_mars_hyperparameters <- function(
    X, y, sample_index1, sample_index2, signal_index, q,
    permR2Mars_func, SelectFeatures_func, CalculateFDP_Power_func,
    fdr_threshold = 0.10,
    degree_grid = c(1, 2, 3),
    nprune_grid = c(NULL, 10, 20, 30, 40, 50),
    minspan_grid = c(1, 2, 3, 5),
    thresh_grid = c(0, 0.001, 0.01),
    nfold = 5,
    verbose = TRUE) {
  
  # Create hyperparameter grid
  param_grid <- expand.grid(
    degree = degree_grid,
    nprune = nprune_grid,
    minspan = minspan_grid,
    thresh = thresh_grid,
    stringsAsFactors = FALSE
  )
  
  n_configs <- nrow(param_grid)
  if (verbose) {
    cat(sprintf("Testing %d hyperparameter configurations...\n", n_configs))
  }
  
  # Storage for results
  results <- vector("list", n_configs)
  
  # Evaluate each configuration
  for (i in 1:n_configs) {
    if (verbose && i %% 10 == 0) {
      cat(sprintf("  Progress: %d/%d configurations tested\n", i, n_configs))
    }
    
    config <- param_grid[i, ]
    
    res <- evaluate_mars_config(
      X = X, y = y,
      sample_index1 = sample_index1,
      sample_index2 = sample_index2,
      signal_index = signal_index,
      q = q,
      degree = config$degree,
      nprune = config$nprune,
      minspan = config$minspan,
      thresh = config$thresh,
      nfold = nfold,
      permR2Mars_func = permR2Mars_func,
      SelectFeatures_func = SelectFeatures_func,
      CalculateFDP_Power_func = CalculateFDP_Power_func
    )
    
    results[[i]] <- c(config, res)
  }
  
  # Convert to data frame
  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      degree = x$degree,
      nprune = ifelse(is.null(x$nprune), NA, x$nprune),
      minspan = x$minspan,
      thresh = x$thresh,
      fdp = x$fdp,
      power = x$power,
      R2_1 = x$R2_1,
      R2_2 = x$R2_2,
      n_selected = x$n_selected,
      n_terms = x$n_terms,
      converged = x$converged,
      stringsAsFactors = FALSE
    )
  }))
  
  # Filter for valid results
  valid_results <- results_df[results_df$converged & !is.na(results_df$fdp), ]
  
  if (nrow(valid_results) == 0) {
    warning("No valid configurations found!")
    return(list(
      best_params = NULL,
      all_results = results_df,
      fdr_controlled_results = NULL
    ))
  }
  
  # Filter for FDR control
  fdr_controlled <- valid_results[valid_results$fdp <= fdr_threshold, ]
  
  if (nrow(fdr_controlled) == 0) {
    warning(sprintf("No configurations achieved FDR <= %.2f", fdr_threshold))
    # Return best overall FDR instead
    best_idx <- which.min(valid_results$fdp)
    best_params <- valid_results[best_idx, ]
  } else {
    # Among FDR-controlled configs, find max power
    best_idx <- which.max(fdr_controlled$power)
    best_params <- fdr_controlled[best_idx, ]
  }
  
  if (verbose) {
    cat("\n=== Tuning Results ===\n")
    cat(sprintf("Total configurations tested: %d\n", n_configs))
    cat(sprintf("Valid configurations: %d\n", nrow(valid_results)))
    cat(sprintf("FDR-controlled (FDR <= %.2f): %d\n", fdr_threshold, nrow(fdr_controlled)))
    cat("\n=== Best Configuration ===\n")
    cat(sprintf("Degree: %d\n", best_params$degree))
    cat(sprintf("Nprune: %s\n", ifelse(is.na(best_params$nprune), "NULL (auto)", best_params$nprune)))
    cat(sprintf("Minspan: %d\n", best_params$minspan))
    cat(sprintf("Thresh: %.4f\n", best_params$thresh))
    cat(sprintf("\nFDP: %.4f\n", best_params$fdp))
    cat(sprintf("Power: %.4f\n", best_params$power))
    cat(sprintf("R2 (fold 1): %.4f\n", best_params$R2_1))
    cat(sprintf("R2 (fold 2): %.4f\n", best_params$R2_2))
    cat(sprintf("# Selected: %d\n", best_params$n_selected))
    cat(sprintf("# MARS terms: %d\n", best_params$n_terms))
  }
  
  # Sort results by power (descending) for FDR-controlled
  if (nrow(fdr_controlled) > 0) {
    fdr_controlled <- fdr_controlled[order(-fdr_controlled$power), ]
  }
  
  return(list(
    best_params = best_params,
    all_results = results_df,
    fdr_controlled_results = fdr_controlled,
    valid_results = valid_results,
    fdr_threshold = fdr_threshold
  ))
}


myresults=tune_mars_hyperparameters(X, y, sample_index1, sample_index2, signal_index, q,
                       permR2Mars_func=permR2Mars, SelectFeatures_func=SelectFeatures, 
                       CalculateFDP_Power_func=CalculateFDP_Power,
                       fdr_threshold = 0.10)
library(openxlsx)
write.xlsx(myresults$valid_results,file='C:/Users/mde4023/Downloads/FDR_Datasplitting/Scenario/Scenario3b/mars_tuning_validresults.xlsx',sheetName = 'results',row.names=FALSE)
