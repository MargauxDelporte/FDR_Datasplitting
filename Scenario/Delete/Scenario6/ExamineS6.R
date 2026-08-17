library(openxlsx)
library(dplyr)

# Set your working directory
mywd='C:/Users/mraga/Downloads/FDR_Datasplitting'
mywd2 <- paste0(mywd,'/Scenario/Scenario6/Temp2/')

# Get list of all Excel files (.xlsx and .xls)
excel_files <- list.files(path = mywd2, pattern = "\\.csv?$", full.names = TRUE)

# Read and combine all files
combined_data <- excel_files %>%
  lapply(read.csv) %>%
  bind_rows()
names(combined_data)[2:4]=c('SignalStrength','FDR','Power')
combined_data$Method <- dplyr::recode(
  combined_data$Method,
  "LinReg DS"             = "Delporte (single split)",
  "LinReg MS"             = "Delporte (50 splits)",
  "DataSplitting"         = "Dai (single split)",
  "MultipleDataSplitting" = "Dai (50 splits)",
  "Knockoff"              = "Knockoff",
  "BH"                    = "Benjamini–Hochberg (BH)"
)

write.xlsx(combined_data,paste0(mywd,'/Scenario/Scenario6/','Scenario6.xlsx'))
s6=combined_data


head(s6)
resultsagg <- s6 %>%
  group_by(Method, SignalStrength) %>%
  summarize(
    Avg_FDR = mean(FDR, na.rm = TRUE),
    Avg_Power = mean(Power, na.rm = TRUE),
    # Add empirical confidence intervals
    FDR_SE = sd(FDR, na.rm = TRUE) / sqrt(sum(!is.na(FDR))),
    Power_SE = sd(Power, na.rm = TRUE) / sqrt(sum(!is.na(Power))),
    FDR_Lower = Avg_FDR - 1.96 * FDR_SE,
    FDR_Upper = Avg_FDR + 1.96 * FDR_SE,
    Power_Lower = Avg_Power - 1.96 * Power_SE,
    Power_Upper = Avg_Power + 1.96 * Power_SE,
    # Alternative: use quantiles for empirical CI
    FDR_Q025 = quantile(FDR, 0.025, na.rm = TRUE),
    FDR_Q975 = quantile(FDR, 0.975, na.rm = TRUE),
    Power_Q025 = quantile(Power, 0.025, na.rm = TRUE),
    Power_Q975 = quantile(Power, 0.975, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  ) %>%
  mutate(Scenario = 'scenario_name')  # Add scenario label

print(resultsagg)
View(resultsagg)
