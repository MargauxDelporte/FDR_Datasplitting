# =========================
# Setup
# =========================
suppressPackageStartupMessages({
  library(earth)
  library(foreach)
  library(doParallel)
  library(dplyr)
  library(purrr)
})

# ---- Your finetune(), patched for data/split seeds & arg duplication ----
finetune <- function(
    myseed, mynk,
    # --- data generation ---
    n = 400, p = 500, p0 = 25,
    delta = 10,
    noise_sd = 1,
    # --- split ---
    amountTrain = 0.5,
    # --- MARS (earth) options ---
    penalty = 2,             # (unused directly — you pass mypenalty below)
    mypmethod = "cv",
    nfold = 5,
    myfast.k = 50,
    trace = 0,
    myfast.beta = 1,
    myminspan,
    mythresh,
    mynprune,
    mypenalty
) {
  # Independent seeds for stages
  data_seed  <- myseed + 1L
  split_seed <- myseed + 2L
  
  # --- choose signal indices ---
  set.seed(myseed)
  signal_index <- sample.int(p, size = p0, replace = FALSE)
  
  # --- simulate data ---
  n1 <- floor(n/2); n2 <- n - n1
  set.seed(data_seed)
  X1 <- matrix(rnorm(n1 * p, mean = 0), n1, p)
  X2 <- matrix(rnorm(n2 * p, mean = 0), n2, p)
  X  <- rbind(X1, X2)
  
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta * sqrt(log(p) / n)) * 100
  
  # Nonlinear quadratic signal + Gaussian noise
  y <- as.vector((X^2 %*% beta_star) + rnorm(n, sd = noise_sd))
  
  # --- prepare frames & split ---
  amountTest <- 1 - amountTrain
  data <- data.frame(y = y, X)
  names(data) <- c("y", paste0("X", seq_len(p)))
  
  size_half <- floor((amountTest / 2) * n)
  
  set.seed(split_seed)
  train_index     <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
  remaining_index <- setdiff(seq_len(n), train_index)
  sample_index1   <- sample(remaining_index, size = size_half, replace = FALSE)
  sample_index2   <- setdiff(remaining_index, sample_index1)
  
  dataTrain <- data[train_index, , drop = FALSE]
  
  # --- fit MARS ---
  mars_fit <- earth(
    y ~ .,
    data      = dataTrain,
    degree    = 2,
    nk        = mynk,
    pmethod   = mypmethod,
    nfold     = 5,
    ncross    = 3,
    penalty   = mypenalty,
    fast.k    = myfast.k,
    fast.beta = myfast.beta,
    minspan   = myminspan,
    thresh    = mythresh,
    nprune    = mynprune,
    trace     = trace
  )
  
  # --- evaluate on two disjoint test halves ---
  pred1 <- predict(mars_fit, newdata = data[sample_index1, , drop = FALSE])
  pred2 <- predict(mars_fit, newdata = data[sample_index2, , drop = FALSE])
  
  y1 <- y[sample_index1]; y2 <- y[sample_index2]
  
  R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
  R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
  
  mean(c(R2orig1, R2orig2))
}

# =========================
# Successive Halving
# =========================
successive_halving <- function(
    grid,                              # data.frame of configs (cols: mynk, mypmethod, mynprune, mythresh, mypenalty, myfast.k, myfast.beta, myminspan)
    seeds_per_stage = c(1, 3, 5),      # # of seeds to average per stage
    keep_frac       = 0.25,            # fraction to retain each round
    n_cores         = max(1, parallel::detectCores() - 1),
    verbose         = TRUE
){
  stopifnot(is.data.frame(grid), nrow(grid) > 0)
  stages <- length(seeds_per_stage)
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(earth))
    NULL
  })
  clusterExport(cl, varlist = c("finetune"), envir = environment())
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  registerDoParallel(cl)
  
  # Helper to score one config with given vector of seeds
  score_config <- function(cfg_row, seed_vec){
    # Average over seeds (higher is better)
    mean(sapply(seed_vec, function(s){
      finetune(
        myseed      = s,
        mynk        = cfg_row$mynk,
        mypmethod   = cfg_row$mypmethod,
        mynprune    = cfg_row$mynprune,
        mythresh    = cfg_row$mythresh,
        mypenalty   = cfg_row$mypenalty,
        myfast.k    = cfg_row$myfast.k,
        myfast.beta = cfg_row$myfast.beta,
        myminspan   = cfg_row$myminspan
      )
    }))
  }
  
  survivors <- grid %>% mutate(cfg_id = seq_len(n()))
  history   <- list()
  
  for(stage in seq_len(stages)){
    k_seeds <- seeds_per_stage[stage]
    # pick deterministic but distinct seeds per stage
    seed_vecs <- lapply(seq_len(nrow(survivors)), function(i){
      base <- 100000 * stage + i
      base + seq_len(k_seeds) - 1L
    })
    
    if (verbose) message(sprintf("Stage %d/%d: evaluating %d configs with %d seed(s) each...",
                                 stage, stages, nrow(survivors), k_seeds))
    
    # Parallel evaluation
    scores <- foreach(i = seq_len(nrow(survivors)), .combine = c, .inorder = TRUE) %dopar% {
      cfg <- survivors[i, , drop = FALSE]
      seedv <- seed_vecs[[i]]
      score_config(cfg, seedv)
    }
    
    scored <- survivors %>% mutate(score = scores, stage = stage)
    history[[stage]] <- scored
    
    # Select top fraction (ties handled by row order)
    K <- max(1, ceiling(nrow(scored) * keep_frac))
    survivors <- scored %>% arrange(desc(score)) %>% slice(1:K) %>% select(-score, -stage)
    
    if (verbose) {
      best <- scored %>% arrange(desc(score)) %>% slice(1)
      message(sprintf("  Best score: %.4f | cfg_id: %d", best$score, best$cfg_id))
    }
  }
  
  results <- bind_rows(history) %>% arrange(stage, desc(score))
  list(
    leaderboard = results %>% arrange(desc(score)),
    final_survivors = survivors,
    stages = stages
  )
}
# =========================
# Example: build a reasonable grid, run SH
# =========================

# You can keep your original large sets, but here's a trimmed example:
nk_grid         <- seq(10, 70, by = 10)
pmethod_grid    <- c('backward', 'forward', 'seqrep', 'exhaustive', 'cv')  # keep 'cv' if you want internal CV
nprune_grid     <- c(10, 25, 40)
thresh_grid     <- c(1e-3, 1e-2)          # drop 1e-4 in practice; keep if needed
penalty_grid    <- c(-1, 2, 3)            # modest set
fastk_grid      <- c(0, 5, 10, 20)
fastb_grid      <- c(0, 1)
minspan_grid    <- c(-3, 0, 1)

grid0 <- expand.grid(
  mynk        = nk_grid,
  mypmethod   = pmethod_grid,
  mynprune    = nprune_grid,
  mythresh    = thresh_grid,
  mypenalty   = penalty_grid,
  myfast.k    = fastk_grid,
  myfast.beta = fastb_grid,
  myminspan   = minspan_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Constraint: fast.beta irrelevant when fast.k == 0 -> force 0 there
grid0$myfast.beta[grid0$myfast.k == 0] <- 0L
grid0 <- distinct(grid0)

# Optional: subsample initial pool (e.g., 1500) to start faster
set.seed(123)
if (nrow(grid0) > 1500) {
  grid0 <- grid0 %>% slice(sample.int(n(), 1500))
}

# Run successive halving
sh_out <- successive_halving(
  grid            = grid0,
  seeds_per_stage = c(1, 3, 5),   # escalate evaluation budget
  keep_frac       = 0.25,         # keep top 25% each round
  n_cores         = max(1, parallel::detectCores() - 1),
  verbose         = TRUE
)

# Inspect results
head(sh_out$leaderboard, 10)
sh_out$final_survivors
