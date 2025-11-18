# Path to your folder
csv_dir <- "C:/Users/mde4023/Downloads/FDR_Datasplitting/Results/ResultsHDNL Scenario/Temp2"
setwd("C:/Users/mde4023/Downloads/FDR_Datasplitting/Results/ResultsHDNL Scenario")
csv_files <- list.files(
  path       = csv_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)
warnings()

# Get full paths of all .csv files
# Read each file into a list of data.frames
data_list <- lapply(csv_files, read.csv, stringsAsFactors = FALSE)
library(dplyr)
library(openxlsx)
all_data <- bind_rows(data_list, .id = "source_file")

Results=all_data
#View(subset(Results2,Results2$Method=='Mars MS'))
##########visualise the results###########
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)
mywd='C:/Users/marga/Downloads/FDR_Datasplitting'
mywd <- paste0(mywd,'/Results')

colors <- c("#000000","#FF00FF","#009900", "#99ccff", "#0000FF", "#FF0000")
Results2=Results
names(Results2)=c('seed','Method','SignalStrength','FDR','Power')
Results2$FDR=round(as.numeric(Results2$FDR),3)
Results2$Power=round(as.numeric(Results2$Power),2)
resultsagg <- Results2 %>%
  group_by(Method, SignalStrength) %>%
  summarize(
    Avg_FDR = mean(FDR),
    Avg_Power = mean(Power),
    N=length(FDR)
  )
resultsagg$Signal_noisy <- as.numeric(resultsagg$SignalStrength) + runif(nrow(resultsagg), -0.2, 0.2)
View(resultsagg)
resultsagg <- resultsagg %>%
  mutate(Method  = case_when(
    Method == "BH" ~ "Benjamini–Hochberg",
    Method == "DataSplitting" ~ "Dai (single split)",
    Method == "MultipleDataSplitting" ~ "Dai (50 splits)",
    Method == "Knockoff" ~ "Knockoff",
    Method == "Mars DS" ~ "Delporte (single split)",
    Method == "Mars MS" ~ "Delporte (50 splits)",
    TRUE ~ Method  # default if none match
  ))
PowerPlot <- ggplot(resultsagg, aes(x = Signal_noisy, y = as.numeric(Avg_Power), color = Method)) +
  geom_point(size = 3) +
  geom_line()+
  labs(x = "Signal", y = "Power") +
  scale_x_continuous(breaks = seq(from = 5, to = 13, by = 1)) +
  geom_hline(yintercept = 0.8) +
  scale_color_manual(values = colors)+
  coord_cartesian(ylim = c(-0.01, 1)) 
FDRPlot=ggplot(resultsagg, aes(x = Signal_noisy, y = as.numeric(Avg_FDR ), color = Method)) +
  geom_point(size = 3) +geom_line()+
  labs(x = "Signal", y = "FDR")+
  scale_x_continuous(breaks=seq(from=5,to=13,by=1))+
  geom_hline(yintercept=0.1)+
  scale_color_manual(values = colors)
PlotPermute=ggarrange(
  PowerPlot, FDRPlot,
  common.legend = TRUE, legend = "right"
)
PlotPermute
ggsave("NLScenario.png",
       plot   = PlotPermute,
       width  = 8,
       height = 8/18*8,
       units  = "in",
       dpi    = 100)