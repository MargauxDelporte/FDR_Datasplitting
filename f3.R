# Function to check missingness comprehensively
check_missingness <- function(data) {
  
  # Calculate missing statistics for each variable
  missing_stats <- data.frame(
    variable = names(data),
    n_missing = sapply(data, function(x) sum(is.na(x))),
    n_total = nrow(data),
    pct_missing = sapply(data, function(x) round(100 * sum(is.na(x)) / length(x), 2)),
    n_present = sapply(data, function(x) sum(!is.na(x))),
    data_type = sapply(data, class)
  )
  
  # Sort by percentage missing (descending)
  missing_stats <- missing_stats[order(-missing_stats$pct_missing), ]
  rownames(missing_stats) <- NULL
  
  return(missing_stats)
}

# Check missingness in your datasets
cat("=== Missingness in Training Data ===\n")
missing_train <- check_missingness(dataTrain)
print(missing_train)

cat("\n=== Missingness in Full X Matrix ===\n")
missing_X <- check_missingness(as.data.frame(X))
print(missing_X)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat(sprintf("Total variables: %d\n", nrow(missing_train)))
cat(sprintf("Variables with ANY missing: %d\n", sum(missing_train$n_missing > 0)))
cat(sprintf("Variables with >5%% missing: %d\n", sum(missing_train$pct_missing > 5)))
cat(sprintf("Variables with >25%% missing: %d\n", sum(missing_train$pct_missing > 25)))
cat(sprintf("Variables with >50%% missing: %d\n", sum(missing_train$pct_missing > 50)))

# Variables with high missingness
if (any(missing_train$pct_missing > 10)) {
  cat("\n=== Variables with >10% Missing ===\n")
  high_missing <- missing_train[missing_train$pct_missing > 10, ]
  print(high_missing)
}