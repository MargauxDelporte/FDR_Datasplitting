library(dplyr)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)

# Path to your folder
csv_dir <- "C:/Users/mde4023/Downloads/FDR_Datasplitting/Scenario/Scenario3/Temp2"
csv_files <- list.files(
  path       = csv_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)

# Get full paths of all .csv files
# Read each file into a list of data.frames
data_list <- lapply(csv_files, read.csv, stringsAsFactors = FALSE)
all_data <- bind_rows(data_list, .id = "source_file")
all_data <- all_data %>%
  mutate(Method  = case_when(
    Method == "Benjamini–Hochberg (BH)" ~ "Benjamini–Hochberg (BH)",
    Method == "Benjamini-Hochberg (BH)" ~ "Benjamini–Hochberg (BH)",
    Method == "DataSplitting" ~ "Dai (single split)",
    Method == "MultipleDataSplitting" ~ "Dai (50 splits)",
    Method == "Knockoff" ~ "Knockoff",
    Method == "Boost DS" ~ "Delporte (single split)",
    Method == "Boost MS" ~ "Delporte (50 splits)",
    TRUE ~ Method  # default if none match
  ))
names(all_data)=c('seed','Method','SignalStrength','FDR','Power')
#write.xlsx(all_data,file='C:/Users/mde4023/Downloads/FDR_Datasplitting/Scenario/Scenario3/Scenario3.xlsx')

resultsagg <- all_data %>%
  group_by(Method, SignalStrength) %>%
  summarize(
    Avg_FDR = mean(FDR),
    Avg_Power = mean(Power),
    N=length(FDR)
  )
View(resultsagg)
