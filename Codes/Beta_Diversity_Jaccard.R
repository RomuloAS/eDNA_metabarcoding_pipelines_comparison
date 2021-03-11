library(dplyr)
library(tidyr)
library(vegan)
library(xtable)
library(ggforce)
library(ranacapa)
library(tikzDevice)
library(concaveman)

# Set working directory to source file location
if(Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}else{
  setwd(utils::getSrcDirectory()[1])
}

# Load data
load("Tables.RData")

##### Statistical analyses using all.pipelines.detailed #####

### Pipeline, Reservoir, and Season

# Remove columns Reservoir, Season and Pipeline
statistic = all.pipelines.detailed[,! colnames(all.pipelines.detailed) %in% c("Reservoir", "Season", "Pipeline")]

### BETA

# Compute dissimilarity indices
beta.jaccard = vegdist(statistic, method = "jaccard")

# Minimum and Maximum values
min(beta.jaccard)
which(as.matrix(beta.jaccard) == min(beta.jaccard), arr.ind=TRUE)
max(beta.jaccard)
which(as.matrix(beta.jaccard) == max(beta.jaccard), arr.ind=TRUE)

# Calculate PCoA principal coordinates analysis
pc.jaccard <- as.data.frame(cmdscale(beta.jaccard, k = 2))

# Create new columns
pc.jaccard$Reservoir_Season_Pipeline = rownames(pc.jaccard)
pc.jaccard = pc.jaccard %>% separate(Reservoir_Season_Pipeline, c("Reservoir", "Season", "Pipeline"), " ", remove = F)

# Multivariate homogeneity of groups dispersions (variances)
beta.disp = betadisper(beta.jaccard, pc.jaccard$Pipeline)

# Tukey Honest Significant Differences
tukey = TukeyHSD(beta.disp)
# Export pair Tukey table
print(xtable(tukey$group, digits = c(1, 3, 3, 3, 3)), booktabs=TRUE, file = "Beta_Diversity_Jaccard_Pipelines_Tukey.tex")

# Pipeline
# Permutational Multivariate Analysis of Variance Using Distance Matrices
permanova =  adonis(as.formula("beta.jaccard~Pipeline"), data = pc.jaccard)
# Export permanova table
print(xtable(permanova$aov.tab), booktabs=TRUE, file = "Beta_Diversity_Jaccard_Pipelines_Permanova.tex")

# Pairwise multilevel comparison using adonis
pair.permanova = pairwise_adonis(statistic, pc.jaccard$Pipeline, sim_method = "jaccard")
# Export pair permanova table
print(xtable(pair.permanova, digits = c(0, 1, 2, 2, 3, 2, 1)), booktabs=TRUE, file = "Beta_Diversity_Jaccard_Pipelines_Pair_Permanova.tex")

# Reservoir
# Permutational Multivariate Analysis of Variance Using Distance Matrices
permanova = adonis(as.formula("beta.jaccard~Reservoir"), data = pc.jaccard)

# Pairwise multilevel comparison using adonis
pair.permanova = pairwise_adonis(statistic, pc.jaccard$Reservoir, sim_method = "jaccard")

# Season
# Permutational Multivariate Analysis of Variance Using Distance Matrices
permanova = adonis(as.formula("beta.jaccard~Season"), data = pc.jaccard)

# Pairwise multilevel comparison using adonis
pair.permanova = pairwise_adonis(statistic, pc.jaccard$Season, sim_method = "jaccard")

# Calculate eigenvalue percentage
eigenvalues = summary(eigenvals(beta.disp, model = "all"))
eigenvalues.percentage = eigenvalues["Proportion Explained", ]
eigenvalues.percentage = round(100 * eigenvalues.percentage, 2)

pc.jaccard = unite(pc.jaccard, "Reservoir Season", Reservoir, Season, sep = " ", remove = F)

# Export plot (Remove % of the plot and put it after export table)
tikz(file = "Beta_Diversity_Jaccard_PCoA_PRS.tex", width = 6, height = 3.7)
# Plot PCoA
ggplot(pc.jaccard, aes(x = V1, y = V2)) + 
  geom_point(aes(shape = Pipeline, fill = `Reservoir Season`, color = `Reservoir Season`), size = 1, stroke = 1) +
  # stat_ellipse(aes(fill = Pipeline, color = Pipeline), geom = "polygon", alpha = 0.3, size = 0.5) +
  #geom_mark_hull(aes(fill = Reservoir_Season), alpha = 0.3, expand = unit(3, "mm")) +
  xlab(paste("PCoA1[", eigenvalues.percentage[["PCoA1"]], "]", sep = "")) +
  ylab(paste("PCoA2[", eigenvalues.percentage[["PCoA2"]], "]", sep = "")) +
  theme(axis.text.x = element_text(angle = 45, size = 6, color = "black",
                                   face = "plain", vjust = 1, hjust = 1),
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
  ) + scale_shape_manual(values = c("Anacapa" = 21,
                                    "Barque" = 22,
                                    "MetaBEAT" = 23,
                                    "MiFish" = 24,
                                    "SEQme" = 25)) +
  scale_fill_manual(values = c("Klíčava Autumn" = "white",
                                "Klíčava Summer" = "#B79F00",
                                "Římov Autumn" = "white",
                                "Římov Summer" = "#00BFC4",
                                "Žlutice Autumn" = "white",
                                "Žlutice Summer" = "#F564E3")) +
  scale_color_manual(values = c("Klíčava Autumn" = "#B79F00",
                               "Klíčava Summer" = "black",
                               "Římov Autumn" = "#00BFC4",
                               "Římov Summer" = "black",
                               "Žlutice Autumn" = "#F564E3",
                               "Žlutice Summer" = "black")) +
  guides(fill=guide_legend(override.aes=list(shape=21)))
dev.off()

### Pipeline

# Remove columns Reservoir, Season
statistic = all.pipelines.detailed[,! colnames(all.pipelines.detailed) %in% c("Reservoir", "Season")]
statistic = aggregate(. ~ Pipeline, data = statistic, FUN = sum)
rownames(statistic) = statistic$Pipeline
statistic$Pipeline = NULL

# Computes dissimilarity indices for pipelines
beta.jaccard = vegdist(statistic, method = "jaccard")
min(beta.jaccard)
max(beta.jaccard)
beta.jaccard

# Export jaccard table
print(xtable(as.matrix(beta.jaccard)), booktabs=TRUE, file = "Beta_Diversity_Jaccard_Pipelines.tex")

### Reservoir

# Remove columns Pipeline, Season
statistic = all.pipelines.detailed[,! colnames(all.pipelines.detailed) %in% c("Pipeline", "Season")]
statistic = aggregate(. ~ Reservoir, data = statistic, FUN = sum)
rownames(statistic) = statistic$Reservoir
statistic$Reservoir = NULL

# Computes dissimilarity indices for Reservoirs
beta.jaccard = vegdist(statistic, method = "jaccard")
min(beta.jaccard)
max(beta.jaccard)
beta.jaccard

# Export jaccard table
print(xtable(as.matrix(beta.jaccard)), booktabs=TRUE, file = "Beta_Diversity_Jaccard_Reservoirs.tex")

### Season

# Remove columns Reservoir, Pipeline
statistic = all.pipelines.detailed[,! colnames(all.pipelines.detailed) %in% c("Reservoir", "Pipeline")]
statistic = aggregate(. ~ Season, data = statistic, FUN = sum)
rownames(statistic) = statistic$Season
statistic$Season = NULL

# Computes dissimilarity indices for Seasons
beta.jaccard = vegdist(statistic, method = "jaccard")
beta.jaccard