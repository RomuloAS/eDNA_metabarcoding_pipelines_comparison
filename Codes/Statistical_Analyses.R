# Vector with package names
packages = c("vegan", "hdd", "tidyr", "dplyr", "ggplot2", "colorblindr", "tidyverse",
             "tikzDevice", "xtable", "reshape2", "scales", "dataMaid", "nortest", "xtable",
             "ggforce", "ranacapa", "concaveman")

# Stop scientific notation
options(scipen=10000)

# Load packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(packages, character.only = T)

# Set working directory to source file location
if(Sys.getenv("RSTUDIO") == "1"){
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}else{
    setwd(utils::getSrcDirectory()[1])
}

load("../Taxonomic_Identification/datasets.RData")

# Create folder to save the results
#datetime = format(Sys.time(), '_%Y%m%d_%H%M%S')
#RESULT = paste0("Statistic", datetime)
#dir.create(RESULT)

# Read file
#path = file.choose(new = FALSE)
#delim = guess_delim(path)
#all.pipelines.detailed.rarefied = read.csv(path, sep = delim,
#                                header = T, check.names = F,
#                                row.names = 1,
#                                stringsAsFactors = FALSE)

#statistic = data.frame( t(all.pipelines.detailed.rarefied) )
statistic = all.pipelines.detailed.rarefied[,! colnames(all.pipelines.detailed.rarefied) %in% c("Reservoir", "Season", "Pipeline")]

### Alpha diversity: observed species
all.pipelines.richness = all.pipelines.detailed.rarefied
all.pipelines.richness$Richness = rowSums(statistic > 0)
all.pipelines.richness = unite(all.pipelines.richness, Reservoir_Season_Pipeline, c(Reservoir, Season, Pipeline), remove=F, sep=" ")
all.pipelines.richness = unite(all.pipelines.richness, Reservoir_Season, c(Reservoir, Season), remove=F, sep=" ")
all.pipelines.richness = unite(all.pipelines.richness, Reservoir_Pipeline, c(Reservoir, Pipeline), remove=F, sep=" ")
all.pipelines.richness = unite(all.pipelines.richness, Season_Pipeline, c(Season, Pipeline), remove=F, sep=" ")
all.pipelines.richness$All = "All"

ggplot(all.pipelines.richness, aes(x = Pipeline, y = Richness)) + 
  geom_point(aes(fill = Pipeline), shape = 21, size = 3) +
  xlab(NULL) +
  ylab("Species richness") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, face = "plain"),
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
  ) + facet_grid(~Reservoir_Season, switch = "x") + 
  scale_fill_discrete(name = "Pipelines")


# Plot All
p = ggplot(all.pipelines.richness, aes(y = Richness)) +
  stat_boxplot(geom = "errorbar", lwd = 1, position = "dodge", show.legend = F, aes(color = All)) +
  geom_boxplot(coef = 1.5, show.legend = F, alpha = 0.5, outlier.alpha = 1, varwidth = T, lwd = 1, aes(color = All, fill = All)) +
  xlab(NULL) +
  ylab("Species richness") +
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
min(all.pipelines.richness$Richness)
max(all.pipelines.richness$Richness)
median(all.pipelines.richness$Richness)
quartiles(all.pipelines.richness$Richness, maxDecimals = 0)
all.pipelines.richness[all.pipelines.richness$Richness == min(all.pipelines.richness$Richness),]["Richness"]
all.pipelines.richness[all.pipelines.richness$Richness == max(all.pipelines.richness$Richness),]["Richness"]

# fit linear models
mod.richness = aov(Richness~Pipeline, data=all.pipelines.richness)
# ANOVA 
anova.test = anova(mod.richness)
anova.test
# Tukey
tukey.test = TukeyHSD(mod.richness)
tukey.test
# Check if p-value < 0.05
tukey.test = as.data.frame(tukey.test[["Pipeline"]])
tukey.test[tukey.test$`p adj` < 0.05,]

### END

### Alpha diversity: shannon index
all.pipelines.shannon = all.pipelines.detailed.rarefied
all.pipelines.shannon$Shannon = diversity(statistic)
all.pipelines.shannon = unite(all.pipelines.shannon, Reservoir_Season_Pipeline, c(Reservoir, Season, Pipeline), remove=F, sep=" ")
all.pipelines.shannon = unite(all.pipelines.shannon, Reservoir_Season, c(Reservoir, Season), remove=F, sep=" ")
all.pipelines.shannon = unite(all.pipelines.shannon, Reservoir_Pipeline, c(Reservoir, Pipeline), remove=F, sep=" ")
all.pipelines.shannon = unite(all.pipelines.shannon, Season_Pipeline, c(Season, Pipeline), remove=F, sep=" ")
all.pipelines.shannon$All = "All"

ggplot(all.pipelines.shannon, aes(x = Pipeline, y = Shannon)) + 
  geom_point(aes(fill = Pipeline), shape = 21, size = 3) +
  xlab(NULL) +
  ylab("Shannon index") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, face = "plain"),
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
  ) + facet_grid(~Reservoir_Season, switch = "x") + 
  scale_fill_discrete(name = "Pipelines")

p = ggplot(all.pipelines.shannon, aes(y = Shannon)) +
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
min(all.pipelines.shannon$Shannon)
max(all.pipelines.shannon$Shannon)
median(all.pipelines.shannon$Shannon)
quartiles(all.pipelines.shannon$Shannon)
all.pipelines.shannon[all.pipelines.shannon$Shannon == min(all.pipelines.shannon$Shannon),]["Shannon"]
all.pipelines.shannon[all.pipelines.shannon$Shannon == max(all.pipelines.shannon$Shannon),]["Shannon"]

# fit linear models
mod.Shannon = aov(Shannon~Pipeline, data=all.pipelines.shannon)
# ANOVA 
anova.test = anova(mod.Shannon)
anova.test
# Tukey
tukey.test = TukeyHSD(mod.Shannon)
tukey.test
# Check if p-value < 0.05
tukey.test = as.data.frame(tukey.test[["Pipeline"]])
tukey.test[tukey.test$`p adj` < 0.05,]

### END

### BETA Bray
beta.bray = vegdist(statistic, method = "bray")

min(beta.bray)
which(as.matrix(beta.bray) == min(beta.bray), arr.ind=TRUE)
max(beta.bray)
which(as.matrix(beta.bray) == max(beta.bray), arr.ind=TRUE)

pc.bray <- as.data.frame(cmdscale(beta.bray, k = 2))
pc.bray$Reservoir_Season_Pipeline = rownames(pc.bray)
pc.bray = pc.bray %>% separate(Reservoir_Season_Pipeline, c("Reservoir", "Season", "Pipeline"), " ")

permanova =  adonis(as.formula("beta.bray~Season"), data = pc.bray)
pair.permanova = pairwise_adonis(statistic, pc.bray$Reservoir, sim_method = "bray")
beta.disp = betadisper(beta.bray, pc.bray$Pipeline)
tukey = TukeyHSD(beta.disp)

eigenvalues = summary(eigenvals(beta.disp, model = "all"))
eigenvalues.percentage = eigenvalues["Proportion Explained", ]
eigenvalues.percentage = round(100 * eigenvalues.percentage, 2)

# Plot PCoA
ggplot(pc.bray, aes(x = V1, y = V2)) + 
  geom_point(aes(fill = Pipeline), shape = 21, size = 3) +
  geom_mark_hull(aes(fill = Pipeline, color = Pipeline), alpha = 0.3, expand = unit(3, "mm")) +
  xlab(paste("PCoA1[", eigenvalues.percentage[["PCoA1"]], "%]", sep = "")) +
  ylab(paste("PCoA2[", eigenvalues.percentage[["PCoA2"]], "%]", sep = "")) +
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
        panel.background = element_rect(fill = 'lightgray'),
        panel.grid.major = element_line(colour = "white", size = 0.3),
        panel.grid.minor = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

### END

### Beta Jaccard
beta.jaccard = vegdist(statistic, method = "jaccard", binary = T)

min(beta.jaccard)
which(as.matrix(beta.jaccard) == min(beta.jaccard), arr.ind=TRUE)
max(beta.jaccard)
which(as.matrix(beta.jaccard) == max(beta.jaccard), arr.ind=TRUE)

pc.jaccard <- as.data.frame(cmdscale(beta.jaccard, k = 2))
pc.jaccard$Reservoir_Season_Pipeline = rownames(pc.jaccard)
pc.jaccard = pc.jaccard %>% separate(Reservoir_Season_Pipeline, c("Reservoir", "Season", "Pipeline"), " ", remove = F)
beta.disp = betadisper(beta.jaccard, pc.jaccard$Pipeline)
tukey = TukeyHSD(beta.disp)
permanova =  adonis(as.formula("beta.jaccard~Pipeline"), data = pc.jaccard)
pair.permanova = pairwise_adonis(statistic, pc.jaccard$Pipeline, sim_method = "jaccard")

eigenvalues = summary(eigenvals(beta.disp, model = "all"))
eigenvalues.percentage = eigenvalues["Proportion Explained", ]
eigenvalues.percentage = round(100 * eigenvalues.percentage, 2)

pc.jaccard = unite(pc.jaccard, "Reservoir Season", Reservoir, Season, sep = " ", remove = F)
# Plot PCoA
ggplot(pc.jaccard, aes(x = V1, y = V2)) + 
  geom_point(aes(shape = Pipeline, fill = `Reservoir Season`, color = `Reservoir Season`), size = 2, stroke = 1) +
  xlab(paste("PCoA1[", eigenvalues.percentage[["PCoA1"]], "%]", sep = "")) +
  ylab(paste("PCoA2[", eigenvalues.percentage[["PCoA2"]], "%]", sep = "")) +
  theme(axis.text.x = element_text(angle = 45, size = 6, color = "black",
                                   face = "plain", vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, face = "plain"),
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

### END

### Mantel

# Extract pipelines
statistic_anacapa <- statistic[c(1,6,11,16,21,26),]
statistic_barque <- statistic[c(1,6,11,16,21,26)+1,]
statistic_metaBeat <- statistic[c(1,6,11,16,21,26)+2,]
statistic_mifish <- statistic[c(1,6,11,16,21,26)+3,]
statistic_seqME <- statistic[c(1,6,11,16,21,26)+4,]

# Beta
beta.jaccard.anacapa = vegdist(statistic_anacapa, method = "jaccard", binary = T)
beta.jaccard.barque = vegdist(statistic_barque, method = "jaccard", binary = T)
beta.jaccard.metaBeat = vegdist(statistic_metaBeat, method = "jaccard", binary = T)
beta.jaccard.mifish = vegdist(statistic_mifish, method = "jaccard", binary = T)
beta.jaccard.seqME = vegdist(statistic_seqME, method = "jaccard", binary = T)

# Plot
plot(beta.jaccard.anacapa, beta.jaccard.barque)
mantel(beta.jaccard.anacapa, beta.jaccard.barque)
plot(beta.jaccard.anacapa, beta.jaccard.metaBeat)
mantel(beta.jaccard.anacapa, beta.jaccard.metaBeat)
plot(beta.jaccard.anacapa, beta.jaccard.mifish)
mantel(beta.jaccard.anacapa, beta.jaccard.mifish)
plot(beta.jaccard.anacapa, beta.jaccard.seqME)
mantel(beta.jaccard.anacapa, beta.jaccard.seqME)
plot(beta.jaccard.barque, beta.jaccard.metaBeat)
mantel(beta.jaccard.barque, beta.jaccard.metaBeat)
plot(beta.jaccard.barque, beta.jaccard.mifish)
mantel(beta.jaccard.barque, beta.jaccard.mifish)
plot(beta.jaccard.barque, beta.jaccard.seqME)
mantel(beta.jaccard.barque, beta.jaccard.seqME)
plot(beta.jaccard.metaBeat, beta.jaccard.mifish)
mantel(beta.jaccard.metaBeat, beta.jaccard.mifish)
plot(beta.jaccard.metaBeat, beta.jaccard.seqME)
mantel(beta.jaccard.metaBeat, beta.jaccard.seqME)
plot(beta.jaccard.mifish, beta.jaccard.seqME)
mantel(beta.jaccard.mifish, beta.jaccard.seqME)
