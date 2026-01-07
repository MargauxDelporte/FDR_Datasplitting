library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggpubr)

#import results
s1=read.xlsx('C:/Users/mde4023/Downloads/FDR_Datasplitting/Scenario/Scenario1/Scenario1.xlsx')
s2=read.xlsx('C:/Users/mde4023/Downloads/FDR_Datasplitting/Scenario/Scenario2/Scenario2.xlsx')
s3=read.xlsx('C:/Users/mde4023/Downloads/FDR_Datasplitting/Scenario/Scenario3/Scenario3.xlsx')
s4=read.xlsx('C:/Users/mde4023/Downloads/FDR_Datasplitting/Scenario/Scenario4/Scenario4.xlsx')

myresults <- list(
  Scenario1 = s1, 
  Scenario2 = s2, 
  Scenario3 = s3, 
  Scenario4 = s4
)

# Store aggregated results
all_agg <- list()
for(scenario_name in names(myresults)){
  cat("\n=====", scenario_name, "=====\n")
  
  resultsagg <- myresults[[scenario_name]] %>%
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
    mutate(Scenario = scenario_name)  # Add scenario label
  
  print(resultsagg)
  all_agg[[scenario_name]] <- resultsagg
}

# Combine all scenarios into one dataframe
combined_results <- bind_rows(all_agg)
print(combined_results)
View(combined_results)

combined_results$Signal_noisy <- as.numeric(combined_results$SignalStrength) + 
  runif(nrow(combined_results), -0.2, 0.2)

colors <- c("#000000","#FF00FF","#009900", "#99ccff", "#0000FF", "#FF0000")

# Power plot with facets and CI ribbons
PowerPlot_all <- ggplot(combined_results, aes(x = Signal_noisy, y = as.numeric(Avg_Power), 
                                              color = Method, fill = Method)) +
  geom_ribbon(aes(ymin = Power_Lower, ymax = Power_Upper), alpha = 0.2, color = NA) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = "Signal Strength", y = "Power") +
  scale_x_continuous(breaks = seq(from = 5, to = 13, by = 1)) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(-0.01, 1)) +
  facet_wrap(~ Scenario, ncol = 4) +
  theme_minimal() +
  theme(strip.text = element_blank())  
PowerPlot_all

# FDR plot with facets and CI ribbons
FDRPlot_all <- ggplot(combined_results, aes(x = Signal_noisy, y = as.numeric(Avg_FDR), 
                                            color = Method, fill = Method)) +
  geom_ribbon(aes(ymin = FDR_Lower, ymax = FDR_Upper), alpha = 0.2, color = NA) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = "", y = "FDR") +
  scale_x_continuous(breaks = seq(from = 5, to = 13, by = 1)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  facet_wrap(~ Scenario, ncol = 4) +
  theme_minimal()

# Combined view
ggarrange(
  FDRPlot_all, PowerPlot_all, 
  nrow = 2,
  common.legend = TRUE, 
  legend = "bottom"
)

# Alternative: Use error bars instead of ribbons
PowerPlot_errorbars <- ggplot(combined_results, aes(x = Signal_noisy, y = as.numeric(Avg_Power), 
                                                    color = Method)) +
  geom_errorbar(aes(ymin = Power_Lower, ymax = Power_Upper), width = 0.2, alpha = 0.5) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = "Signal Strength", y = "Power") +
  scale_x_continuous(breaks = seq(from = 5, to = 13, by = 1)) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_color_manual(values = colors) +
  coord_cartesian(ylim = c(-0.01, 1)) +
  facet_wrap(~ Scenario, ncol = 4) +
  theme_minimal() +
  theme(strip.text = element_blank())  

FDRPlot_errorbars <- ggplot(combined_results, aes(x = Signal_noisy, y = as.numeric(Avg_FDR), 
                                                  color = Method)) +
  geom_errorbar(aes(ymin = FDR_Lower, ymax = FDR_Upper), width = 0.2, alpha = 0.5) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = "", y = "FDR") +
  scale_x_continuous(breaks = seq(from = 5, to = 13, by = 1)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  scale_color_manual(values = colors) +
  facet_wrap(~ Scenario, ncol = 4) +
  theme_minimal()

ggarrange(
  FDRPlot_errorbars, PowerPlot_errorbars, 
  nrow = 2,
  common.legend = TRUE, 
  legend = "bottom"
)
