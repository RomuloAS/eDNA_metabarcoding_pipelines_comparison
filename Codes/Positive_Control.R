library(tidyr)
library(ggplot2)
library(rstatix)
library(EnvStats)
library(normtest)
library(tikzDevice)
library(RVAideMemoire)

# Set working directory to source file location
if(Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}else{
  setwd(utils::getSrcDirectory()[1])
}

# Load data
load("Tables.RData")

# Create column for Reservoir and Season together
all.pipelines.only.positive.control = unite(all.pipelines.only.positive.control, Reservoir_Season, c(Reservoir, Season), remove=F, sep=" ")

# Export ggplot to Latex
tikz(file = "Positive_Control_PRS.tex", width = 6, height = 3)
# Plot Reservoir_Season_Pipeline
ggplot(all.pipelines.only.positive.control, aes(x = Pipeline, y = Percentage_of_the_total)) + 
  geom_point(aes(fill = Pipeline), shape = 21, size = 3) +
  xlab(NULL) +
  ylab("Percentage of assigned reads") +
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

# Minimum percentage
min(all.pipelines.only.positive.control$Percentage_of_the_total)
all.pipelines.only.positive.control[all.pipelines.only.positive.control$Percentage_of_the_total == min(all.pipelines.only.positive.control$Percentage_of_the_total),]

# Maximum percentage
max(all.pipelines.only.positive.control$Percentage_of_the_total)
all.pipelines.only.positive.control[all.pipelines.only.positive.control$Percentage_of_the_total == max(all.pipelines.only.positive.control$Percentage_of_the_total),]

# Mean and median of the percentage list
median(all.pipelines.only.positive.control$Percentage_of_the_total)
mean(all.pipelines.only.positive.control$Percentage_of_the_total)

# Check normality
byf.shapiro(Percentage_of_the_total~Pipeline, data = all.pipelines.only.positive.control)

# fit linear models
mod.percentage = aov(Percentage_of_the_total~Pipeline, data=all.pipelines.only.positive.control)
# ANOVA 
anova.test = anova(mod.percentage)
anova.test
# Tukey
tukey.test = TukeyHSD(mod.percentage)
tukey.test
# Check if p-value < 0.05
tukey.test = as.data.frame(tukey.test[["Pipeline"]])
tukey.test[tukey.test$`p adj` < 0.05,]

# Extract pipelines
pipelines = unique(all.pipelines.only.positive.control$Pipeline)

# For each pipeline do t.test, wilcox, median, median, min and max
for (pipeline in pipelines) {
  dat = all.pipelines.only.positive.control[all.pipelines.only.positive.control$Pipeline == pipeline,]
  
  # Print pipeline name and data
  print(pipeline)
  print(dat$Percentage_of_the_total)
  
  # Statistical tests
  print(paste("chi square: ", round(varTest(dat$Percentage_of_the_total, alternative = "greater", sigma.squared = 5)$p.value, digits = 3)))
  print(paste("t.test: ", t.test(dat$Percentage_of_the_total)$p.value))
  print(paste("wilcox.test: ", wilcox.test(dat$Percentage_of_the_total)$p.value))
  
  # Mean and median of the percentage list
  print(paste("Median: ", median(dat$Percentage_of_the_total)))
  print(paste("Mean: ", mean(dat$Percentage_of_the_total)))
  
  # Min and max of the percentage list
  print(paste("Min: ", min(dat$Percentage_of_the_total)))
  print(paste("Max: ", max(dat$Percentage_of_the_total)))
  
  print("")
}