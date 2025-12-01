library(readr)
library(randomForest)

# ==============================================================================
# 1. File Paths & Data Import
# ==============================================================================

base_dir <- "C:/Users/mde4023/Downloads/FDR_Datasplitting/Case study/yatsunenko/refseq"

task <- read.delim(
  file         = file.path(base_dir, "task-baby-age.txt"),
  header       = FALSE,
  comment.char = "#"
)

otu <- read_delim(
  file         = file.path(base_dir, "otutable.txt"),
  delim        = "\t",
  escape_double = FALSE,
  trim_ws       = TRUE
)

taxa <- read_delim(
  file         = file.path(base_dir, "taxatable.txt"),
  delim        = "\t",
  escape_double = FALSE,
  trim_ws       = TRUE
)

# ==============================================================================
# 2. OTU Pre-processing (prevalence + relative abundance filter)
# ==============================================================================

# Samples as rows, OTUs as columns
otu_t <- t(as.matrix(otu[, -1]))   # drop OTU ID column, transpose

# Prevalence per OTU
prev <- colSums(otu_t > 0) / nrow(otu_t)

# Relative abundances
otu_rel <- otu_t / rowSums(otu_t)

# Filtering thresholds
min_prev        <- 0.10   # present in ≥10% of samples
min_relab       <- 1e-4   # minimum relative abundance
min_relab_prev  <- 0.05   # abundance > min_relab in ≥5% of samples

keep_prev  <- prev >= min_prev
keep_abund <- colSums(otu_rel > min_relab) / nrow(otu_rel) >= min_relab_prev
keep       <- keep_prev & keep_abund

otu_filtered <- otu_t[, keep, drop = FALSE]

# ==============================================================================
# 3. Merge with phenotype (age) & basic setup
# ==============================================================================

names(task) <- c("SampleID", "Age")

# Check alignment
stopifnot(all(task$SampleID %in% rownames(otu_t)))

mydata <- merge(
  x     = task,
  y     = otu_filtered,
  by.x  = "SampleID",
  by.y  = "row.names"
)

# Drop SampleID after merge
mydata <- mydata[, -1]

n <- nrow(mydata)
p <- ncol(mydata) - 1   # number of OTUs

# Design matrix X and outcome y
X <- mydata[, -1, drop = FALSE]
colnames(X) <- paste0("X", seq_len(p))

y <- mydata$Age

# Combined data frame for modeling
dat <- data.frame(y = y, X)

# ==============================================================================
# 4. Train/Test Split
# ==============================================================================

set.seed(20251127)

train_fraction <- 0.5
train_index    <- sample.int(n, size = floor(train_fraction * n), replace = FALSE)
test_index     <- setdiff(seq_len(n), train_index)

trainDat <- dat[train_index, ]
testDat  <- dat[test_index, ]

# ==============================================================================
# 5. Helper: R² function
# ==============================================================================

R2_fun <- function(y, yhat) {
  1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}

# ==============================================================================
# 6. Random Forest Hyperparameter Grid
# ==============================================================================

p_rf <- ncol(trainDat) - 1  # number of predictors (should equal p)

mtry_grid     <- unique(pmax(1, round(c(sqrt(p_rf), p_rf / 4, p_rf / 2))))
ntree_grid    <- c(100,250,500, 1000)
nodesize_grid <- c(1, 5, 10,15)

grid <- expand.grid(
  mtry     = mtry_grid,
  ntree    = ntree_grid,
  nodesize = nodesize_grid
)

# ==============================================================================
# 7. Grid Search: Maximize Test R²
# ==============================================================================

results <- data.frame()
best_R2 <- -Inf
best_rf <- NULL

for (i in seq_len(nrow(grid))) {
  params <- grid[i, ]
  
  set.seed(123 + i)  # reproducible RF fits
  
  rf_mod <- randomForest(
    y ~ .,
    data     = trainDat,
    mtry     = params$mtry,
    ntree    = params$ntree,
    nodesize = params$nodesize
  )
  
  pred_test <- predict(rf_mod, newdata = testDat)
  R2_test   <- R2_fun(testDat$y, pred_test)
  
  results <- rbind(
    results,
    data.frame(
      mtry     = params$mtry,
      ntree    = params$ntree,
      nodesize = params$nodesize,
      R2_test  = R2_test
    )
  )
  
  if (R2_test > best_R2) {
    best_R2 <- R2_test
    best_rf <- rf_mod
  }
}

# ==============================================================================
# 8. Best Parameter Combination
# ==============================================================================

results_sorted <- results[order(-results$R2_test), ]
best_pm        <- results_sorted[1, ]
#  mtry ntree nodesize   R2_test
# 2   81   100        1 0.6978718
best_pm
best_R2
best_rf
