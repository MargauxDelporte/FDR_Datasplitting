num_split <- 50
n <-1500
p <- 250
p0 <- 25
q <- 0.1

mywd='C:/Users/mraga/Downloads/FDR_Datasplitting'
setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Scenario/Scenario6/MarsParallelNL.R'))

EvaluateMarsParams <- function(
    degree,
    nk,
    penalty,
    minspan,
    thresh,
    nsim = 100,
    delta = 13,
    f = 5,
    n = 1500,
    p = 250,
    p0 = 25,
    q = 0.10,
    num_split = 50,
    signal_index,
    seeds = seq_len(nsim),
    verbose = TRUE) {
  stopifnot(length(signal_index) == p0)
  
  results <- matrix(
    NA_real_,
    nrow = nsim,
    ncol = 4,
    dimnames = list(
      NULL,
      c("DS_fdp", "DS_power", "MDS_fdp", "MDS_power")
    )
  )
  
  for (sim in seq_len(nsim)) {
    
    s <- seeds[sim]
    set.seed(s)
    
    # -----------------------------
    # Simulate covariance matrix
    # -----------------------------
    n_df <- f * p
    
    D <- matrix(
      rWishart(
        n = 1,
        df = n_df,
        Sigma = diag(p)
      ),
      nrow = p
    )
    
    # -----------------------------
    # Simulate X
    # -----------------------------
    n1 <- floor(n / 2)
    n2 <- n - n1
    
    X1 <- mvtnorm::rmvnorm(
      n1,
      mean = rep(0, p),
      sigma = D
    )
    
    X2 <- mvtnorm::rmvnorm(
      n2,
      mean = rep(0, p),
      sigma = D
    )
    
    X <- rbind(X1, X2)
    
    # -----------------------------
    # Generate true coefficients
    # -----------------------------
    beta_star <- numeric(p)
    
    beta_star[signal_index] <- rnorm(
      p0,
      mean = 0,
      sd = delta * sqrt(log(p) / n)
    ) * 10
    
    # -----------------------------
    # Outcome
    # -----------------------------
    y <- as.numeric(
      X^2 %*% beta_star + rnorm(n)
    )
    
    # -----------------------------
    # Run MARS procedure
    # -----------------------------
    fit <- ApplyMarsTrain_parallel(
      X = X,
      y = y,
      q = q,
      myseed = s + 10000,
      num_split = num_split,
      signal_index = signal_index,
      degree = degree,
      nk = nk,
      penalty = penalty,
      minspan = minspan,
      thresh = thresh
    )
    
    results[sim, ] <- c(
      fit$DS_fdp,
      fit$DS_power,
      fit$MDS_fdp,
      fit$MDS_power
    )
    
    if (verbose) {
      cat(
        "\rSimulation", sim, "/", nsim,
        "| MDS FDR:", round(mean(results[1:sim, "MDS_fdp"]), 3),
        "| MDS Power:", round(mean(results[1:sim, "MDS_power"]), 3)
      )
      flush.console()
    }
  }
  
  if (verbose) cat("\n")
  
  # -----------------------------
  # Summarize simulations
  # -----------------------------
  summary_df <- data.frame(
    degree = degree,
    nk = nk,
    penalty = penalty,
    minspan = minspan,
    thresh = thresh,
    
    DS_FDR = mean(results[, "DS_fdp"], na.rm = TRUE),
    DS_Power = mean(results[, "DS_power"], na.rm = TRUE),
    
    MDS_FDR = mean(results[, "MDS_fdp"], na.rm = TRUE),
    MDS_Power = mean(results[, "MDS_power"], na.rm = TRUE),
    
    MDS_FDR_SE =
      sd(results[, "MDS_fdp"], na.rm = TRUE) / sqrt(nsim),
    
    MDS_Power_SE =
      sd(results[, "MDS_power"], na.rm = TRUE) / sqrt(nsim)
  )
  
  # -----------------------------
  # Return BOTH full per-sim results and summary
  # -----------------------------
  list(
    summary = summary_df,
    all_results = as.data.frame(results)
  )
}

set.seed(456)

signal_index <- sample(
  seq_len(p),
  size = p0,
  replace = FALSE
)



res1 <- EvaluateMarsParams(
  degree = 1,
  nk = 50,
  penalty = 3,
  minspan = 10,
  thresh = 0.001,
  
  nsim = 100,
  delta = 13,
  f = 5,
  
  n = n,
  p = p,
  p0 = p0,
  q = q,
  num_split = num_split,
  
  signal_index = signal_index
)
mars_grid <- expand.grid(
  degree = c(1, 2),
  nk = c(30, 50),
  penalty = c(2, 4),
  minspan = c(5, 20),
  thresh = c(0.001, 0.005)
)

res1
res1$full

res2 <- EvaluateMarsParams(
  degree = 2,
  nk = 75,
  penalty = 2,
  minspan = 5,
  thresh = 0.0005,
  
  nsim = 100,
  delta = 13,
  f = 5,
  
  n = n,
  p = p,
  p0 = p0,
  q = q,
  num_split = num_split,
  
  signal_index = signal_index
)

# All 100 simulation results are preserved here:
res2$all_results
View(res2$all_results)
# Summary (means/SEs) is here:
res2$summary

write.csv(res2$all_results,file=paste0(mywd,'/Scenario/Scenario6/results'))