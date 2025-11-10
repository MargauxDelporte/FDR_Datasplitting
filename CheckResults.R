# Set working directory
mywd <- "C:/Users/mde4023/Downloads/FDR_Datasplitting/Temp"

# List all CSV files in the folder
files <- list.files(path = mywd, pattern = "\\.csv$", full.names = TRUE)

# Read all CSVs and combine them into one data frame
data_list <- lapply(files, read.csv)
combined_df <- do.call(rbind, data_list)

# Identify columns that contain "Method" in their name
library(dplyr)

summary_stats <- combined_df %>%
  group_by(Method, Delta) %>%
  summarise(
    mean_FDP   = mean(FDP,   na.rm = TRUE),
    mean_Power = mean(Power, na.rm = TRUE),
    count      = n(),
    .groups    = "drop"
  )
table(summary_stats$Method)
summary_stats[summary_stats$Method=='Mars MS',]
summary_stats[summary_stats$Method=='Mars MS',]

# Calculate mean for each "Method" column (ignoring NAs)
method_means <- sapply(combined_df[method_cols], function(x) mean(x, na.rm = TRUE))

# Print results
print(method_means)


# other method
library(readr)
library(readxl)
ResultsNonlinearScenario2 <- read_csv("Results/Old/ResultsNonlinearScenario.csv")
summary_stats_comp <- ResultsNonlinearScenario2 %>%
  group_by(Method,  Delta) %>%
  summarise(
    mean_FDP   = mean(FDP,   na.rm = TRUE),
    mean_Power = mean(Power, na.rm = TRUE),
    count      = n(),
    .groups    = "drop"
  )
summary_stats_comp[summary_stats_comp$Method=='Boost MS',]
summary_stats_comp[summary_stats_comp$Method=='Boost MS',]
setwd(mywd)
