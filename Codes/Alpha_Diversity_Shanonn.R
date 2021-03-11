library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)
library(colorblindr)
library(tidyverse)
library(tikzDevice)
library(xtable)
library(reshape2)
library(scales)
library(dataMaid)
library(nortest)


# Set working directory to source file location
if(Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}else{
  setwd(utils::getSrcDirectory()[1])
}

# Load data
load("Tables.RData")

##### Statistical analyses using all.pipelines.detailed #####

### ALPHA ###

# Remove columns Reservoir, Season and Pipeline
statistic = all.pipelines.detailed[,! colnames(all.pipelines.detailed) %in% c("Reservoir", "Season", "Pipeline")]

# Alpha diversity: shannon index
all.pipelines.detailed$Shannon = diversity(statistic)

# Create column for Reservoir, Season and Pipeline together
all.pipelines.detailed = unite(all.pipelines.detailed, Reservoir_Season_Pipeline, c(Reservoir, Season, Pipeline), remove=F, sep=" ")

# Create column for Reservoir and Season together
all.pipelines.detailed = unite(all.pipelines.detailed, Reservoir_Season, c(Reservoir, Season), remove=F, sep=" ")

# Create column for Reservoir and Pipeline together
all.pipelines.detailed = unite(all.pipelines.detailed, Reservoir_Pipeline, c(Reservoir, Pipeline), remove=F, sep=" ")

# Create column for Season and Pipeline together
all.pipelines.detailed = unite(all.pipelines.detailed, Season_Pipeline, c(Season, Pipeline), remove=F, sep=" ")

# Create column with unique name
all.pipelines.detailed$All = "All"

### Shannon index ###

# Export ggplot to Latex
tikz(file = "Alpha_diversity_shannon_PRS.tex", width = 6, height = 3)
# Plot Reservoir_Season_Pipeline
ggplot(all.pipelines.detailed, aes(x = Pipeline, y = Shannon)) + 
  geom_point(aes(fill = Pipeline), shape = 21, size = 3) +
  xlab(NULL) +
  ylab("Shannon index") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, face = "plain"),
        legend.margin = margin(t = 0, unit='cm'),
        legend.key = element_rect(fill = NA, color = NA),
        strip.text.x = element_text(size = 6, face = "plain", color = "black"),
        strip.background = element_blank(),
        #plot.margin = margin(10, 10, 10, 50),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "lightgray", size = 0.3),
        panel.grid.minor = element_line(colour = "lightgray", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  ) + facet_grid(~Reservoir_Season, switch = "x") + 
  scale_fill_discrete(name = "Pipelines")
dev.off()

# Plot All
p = ggplot(all.pipelines.detailed, aes(y = Shannon)) +
  stat_boxplot(geom = "errorbar", lwd = 1, position = "dodge", show.legend = F, aes(color = All)) +
  geom_boxplot(coef = 1.5, show.legend = F, alpha = 0.5, outlier.alpha = 1, varwidth = T, lwd = 1, aes(color = All, fill = All)) +
  xlab(NULL) +
  ylab("Shannon index") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 10, face = "plain"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, face = "plain"),
        legend.margin = margin(t = 0, unit='cm'),
        legend.key = element_rect(fill = NA, color = NA),
        strip.text.x = element_text(size = 8, face = "plain", color = "black"),
        strip.background = element_blank(),
        #plot.margin = margin(10, 10, 10, 50),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'lightgray'),
        panel.grid.major = element_line(colour = "white", size = 0.3),
        panel.grid.minor = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  ) + facet_grid(~All, switch="x")

# Calculate median and quartiles
ggplot_build(p)$data
min(all.pipelines.detailed$Shannon)
max(all.pipelines.detailed$Shannon)
median(all.pipelines.detailed$Shannon)
quartiles(all.pipelines.detailed$Shannon)
all.pipelines.detailed[all.pipelines.detailed$Shannon == min(all.pipelines.detailed$Shannon),]["Shannon"]
all.pipelines.detailed[all.pipelines.detailed$Shannon == max(all.pipelines.detailed$Shannon),]["Shannon"]

## Pipelines ##

# Plot Pipeline
p = ggplot(all.pipelines.detailed, aes(y = Shannon)) +
  stat_boxplot(geom = "errorbar", lwd = 1, position = "dodge", show.legend = F, aes(color = Pipeline)) +
  geom_boxplot(coef = 1.5, show.legend = F, alpha = 0.5, outlier.alpha = 1, varwidth = T, lwd = 1, aes(color = Pipeline, fill = Pipeline)) +
  xlab(NULL) +
  ylab("Shannon index") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 10, face = "plain"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, face = "plain"),
        legend.margin = margin(t = 0, unit='cm'),
        legend.key = element_rect(fill = NA, color = NA),
        strip.text.x = element_text(size = 8, face = "plain", color = "black"),
        strip.background = element_blank(),
        #plot.margin = margin(10, 10, 10, 50),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "lightgray", size = 0.3),
        panel.grid.minor = element_line(colour = "lightgray", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  ) + facet_grid(~Pipeline, switch="x")

# Export ggplot to Latex
tikz(file = "Alpha_diversity_shannon_pipelines.tex", width = 6, height = 3)
p
dev.off()

# Get min, max, median, and quartiles for each pipeline
ggplot_build(p)$data
# Get pipelines
pipelines = unique(all.pipelines.detailed$Pipeline)
for (pipeline in pipelines) {
  print(pipeline)
  # Extract pipeline
  pipeline = all.pipelines.detailed[all.pipelines.detailed$Pipeline == pipeline,]
  # Get min, max, median, and quartiles for each pipeline
  print(round(min(pipeline$Shannon), digits = 3))
  print(round(max(pipeline$Shannon), digits = 3))
  print(median(pipeline$Shannon))
  print(quartiles(pipeline$Shannon))
  print(pipeline[pipeline$Shannon == min(pipeline$Shannon),]["Shannon"])
  print(pipeline[pipeline$Shannon == max(pipeline$Shannon),]["Shannon"])
  print("--------------------------------------------------------------")
  print("")
}

## Statistical analyses ##

# fit linear models
mod.Shannon = aov(Shannon~Pipeline, data=all.pipelines.detailed)
# ANOVA 
anova.test = anova(mod.Shannon)
anova.test
# Tukey
tukey.test = TukeyHSD(mod.Shannon)
tukey.test

# Check if p-value < 0.05
tukey.test = as.data.frame(tukey.test[["Pipeline"]])
tukey.test[tukey.test$`p adj` < 0.05,]

# Export as latex tables
print(xtable(anova.test, digits = c(0, 0, 4, 4, 4, 4)), booktabs=TRUE, file = "Alpha_Diversity_Shannon_ANOVA_Pipelines.tex")
print(xtable(tukey.test$Pipeline, digits = c(0, 4, 4, 4, 4)), booktabs=TRUE, file = "Alpha_Diversity_Shannon_Tukey_Pipelines.tex")

## Reservoirs ##

# Plot Reservoir
p = ggplot(all.pipelines.detailed, aes(y = Shannon)) +
  stat_boxplot(geom = "errorbar", lwd = 1, position = "dodge", show.legend = F, aes(color = Reservoir)) +
  geom_boxplot(coef = 1.5, show.legend = F, alpha = 0.5, outlier.alpha = 1, varwidth = T, lwd = 1, aes(color = Reservoir, fill = Reservoir)) +
  xlab(NULL) +
  ylab("Shannon index") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 10, face = "plain"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, face = "plain"),
        legend.margin = margin(t = 0, unit='cm'),
        legend.key = element_rect(fill = NA, color = NA),
        strip.text.x = element_text(size = 8, face = "plain", color = "black"),
        strip.background = element_blank(),
        #plot.margin = margin(10, 10, 10, 50),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "lightgray", size = 0.3),
        panel.grid.minor = element_line(colour = "lightgray", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  ) + facet_grid(~Reservoir, switch="x")

# Export ggplot to Latex
tikz(file = "Alpha_diversity_shannon_reservoirs.tex", width = 6, height = 3)
p
dev.off()

# Get min, max, median, and quartiles for each reservoir
ggplot_build(p)$data
# Get reservoirs
reservoirs = unique(all.pipelines.detailed$Reservoir)
for (reservoir in reservoirs) {
  print(reservoir)
  # Extract reservoir
  reservoir = all.pipelines.detailed[all.pipelines.detailed$Reservoir == reservoir,]
  # Get min, max, median, and quartiles for each reservoir
  print(round(min(reservoir$Shannon), digits = 3))
  print(round(max(reservoir$Shannon), digits = 3))
  print(median(reservoir$Shannon))
  print(quartiles(reservoir$Shannon))
  print(reservoir[reservoir$Shannon == min(reservoir$Shannon),]["Shannon"])
  print(reservoir[reservoir$Shannon == max(reservoir$Shannon),]["Shannon"])
  print("--------------------------------------------------------------")
  print("")
}

## Statistical analyses ##

# fit linear models
mod.Shannon = aov(Shannon~Reservoir, data=all.pipelines.detailed)
# ANOVA 
anova.test = anova(mod.Shannon)
anova.test
# Tukey
tukey.test = TukeyHSD(mod.Shannon)
tukey.test

# Check if p-value < 0.05
tukey.test = as.data.frame(tukey.test[["Reservoir"]])
tukey.test[tukey.test$`p adj` < 0.05,]

# Export as latex tables
print(xtable(anova.test, digits = c(0, 0, 4, 4, 4, 4)), booktabs=TRUE, file = "Alpha_Diversity_Shannon_ANOVA_Reservoirs.tex")

## Seasons ##

# Plot Season
p = ggplot(all.pipelines.detailed, aes(y = Shannon)) +
  stat_boxplot(geom = "errorbar", lwd = 1, position = "dodge", show.legend = F, aes(color = Season)) +
  geom_boxplot(coef = 1.5, show.legend = F, alpha = 0.5, outlier.alpha = 1, varwidth = T, lwd = 1, aes(color = Season, fill = Season)) +
  xlab(NULL) +
  ylab("Shannon index") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 10, face = "plain"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, face = "plain"),
        legend.margin = margin(t = 0, unit='cm'),
        legend.key = element_rect(fill = NA, color = NA),
        strip.text.x = element_text(size = 8, face = "plain", color = "black"),
        strip.background = element_blank(),
        #plot.margin = margin(10, 10, 10, 50),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "lightgray", size = 0.3),
        panel.grid.minor = element_line(colour = "lightgray", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  ) + facet_grid(~Season, switch="x")

# Export ggplot to Latex
tikz(file = "Alpha_diversity_shannon_seasons.tex", width = 6, height = 3)
p
dev.off()

# Get min, max, median, and quartiles for each season
ggplot_build(p)$data
# Get seasons
seasons = unique(all.pipelines.detailed$Season)
for (season in seasons) {
  print(season)
  # Extract season
  season = all.pipelines.detailed[all.pipelines.detailed$Season == season,]
  # Get min, max, median, and quartiles for each season
  print(round(min(season$Shannon), digits = 3))
  print(round(max(season$Shannon), digits = 3))
  print(median(season$Shannon))
  print(quartiles(season$Shannon))
  print(season[season$Shannon == min(season$Shannon),]["Shannon"])
  print(season[season$Shannon == max(season$Shannon),]["Shannon"])
  print("--------------------------------------------------------------")
  print("")
}

# Fences
quartiles = quartiles(all.pipelines.detailed[all.pipelines.detailed$Season == "Summer",]$Shannon)
upperq = round(quartiles$value[["75%"]])
lowerq = round(quartiles$value[["25%"]])
iqr = upperq - lowerq
upper.fence = upperq + (1.5 * iqr)
lower.fence = lowerq - (1.5 * iqr)

## Statistical analyses ##

# fit linear models
mod.Shannon = aov(Shannon~Season, data=all.pipelines.detailed)
# ANOVA 
anova.test = anova(mod.Shannon)
anova.test
# T-test
t.test(Shannon~Season, data=all.pipelines.detailed, var.equal = TRUE)
# Tukey
tukey.test = TukeyHSD(mod.Shannon)
tukey.test

# Check if p-value < 0.05
tukey.test = as.data.frame(tukey.test[["Season"]])
tukey.test[tukey.test$`p adj` < 0.05,]

# Export as latex tables
print(xtable(anova.test, digits = c(0, 0, 4, 4, 4, 4)), booktabs=TRUE, file = "Alpha_Diversity_Shannon_ANOVA_Seasons.tex")