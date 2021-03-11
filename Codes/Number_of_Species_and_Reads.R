library(tidyverse)
library(ggplot2)
library(colorblindr)
library(tikzDevice)
library(xtable)
library(reshape2)
library(scales)
library(ggpubr)
library(gridExtra)

# Load R object
load("Tables.RData")

# Ignore scientific notation
options(scipen=10000)

# Read table
# dat = read.csv("Tables/All_Pipelines_detailed.csv", sep = '\t', header=T, row.names = 1, check.names=F)
dat = all.pipelines.detailed

# Species order
species_order = c("Lampetra planeri", "Acipenser-sp.", "Anguilla anguilla", "Phoxinus phoxinus", "Rutilus rutilus", "Chondrostoma nasus", "Squalius cephalus", "Alburnus alburnus", "Blicca+Vimba", "Abramis brama", "L.idus+leuciscus", "Aspius+Scardinius", "Pseudorasbora parva", "Gobio gobio", "Rhodeus amarus", "Tinca tinca", "Hypophthalmichthys molitrix", "Hypophthalmichthys nobilis", "Ctenopharyngodon idella", "Barbus barbus", "Cyprinus carpio", "Carassius auratus", "Carassius carassius", "Barbatula barbatula", "Esox lucius", "Gasterosteus aculeatus", "Sander+Perca", "Gymnocephalus cernua", "Lepomis gibbosus", "Thymallus thymallus", "Oncorhynchus mykiss", "Salvelinus fontinalis", "Salmo trutta", "Coregonus-sp.", "Cottus gobio", "Cottus poecilopus", "Silurus glanis")
species_order = rev(species_order)

#### Count number of species/reads for each pipeline, reservoir, and season.

# Number of reads

# Average of all pipelines 
ignored = dat[ , !(colnames(dat) %in% c("Season", "Reservoir"))]
ignored = aggregate(. ~ Pipeline, data=ignored, FUN=sum)
ignored$Total = rowSums(ignored[, !(colnames(ignored) %in% c("Pipeline"))])
sum(ignored$Total) / length(unique(dat$Pipeline))


# Number of reads

# Pipeline 
ignored = dat[ , !(colnames(dat) %in% c("Season", "Reservoir"))]
ignored = aggregate(. ~ Pipeline, data=ignored, FUN=sum)
ignored$Total = rowSums(ignored[, !(colnames(ignored) %in% c("Pipeline"))])
ignored[c("Pipeline", "Total")]

# Export as latex table
print(xtable(ignored[c("Pipeline", "Total")]), booktabs=TRUE, file = "Number_of_reads_by_Pipeline.tex")

# Reservoir 
ignored = dat[ , !(colnames(dat) %in% c("Pipeline", "Season"))]
ignored = aggregate(. ~ Reservoir, data=ignored, FUN=sum)
ignored$Total = rowSums(ignored[, !(colnames(ignored) %in% c("Reservoir"))])
# Use the code below to calculate the average for Reservoir or Season
ignored$Total = ignored$Total / length(unique(dat$Pipeline))
ignored[c("Reservoir", "Total")]

# Season 
ignored = dat[ , !(colnames(dat) %in% c("Pipeline", "Reservoir"))]
ignored = aggregate(. ~ Season, data=ignored, FUN=sum)
ignored$Total = rowSums(ignored[, !(colnames(ignored) %in% c("Season"))])
# Use the code below to calculate the average for Reservoir and Season
ignored$Total = ignored$Total / length(unique(dat$Pipeline))
ignored[c("Season", "Total")]

# Number of reads for pipeline, reservoir, and season together
total = dat
total$total = rowSums(total[, !(colnames(total) %in% c("Pipeline", "Reservoir", "Season"))])
total = total[c("Pipeline", "Reservoir", "Season", "total")]
total
total[total$total == min(total$total),]
total[total$total == max(total$total),]

# Number of species

# Total
dim(dat[! colnames(dat) %in% c("Pipeline", "Reservoir", "Season")])

# Pipeline 
ignored = dat[ , !(colnames(dat) %in% c("Reservoir", "Season"))]
ignored = aggregate(. ~ Pipeline, data=ignored, FUN=sum)
ignored$Total = rowSums(ignored[, !(colnames(ignored) %in% c("Pipeline"))] != 0)
ignored[c("Pipeline", "Total")]

# Export as latex table
#print(xtable(ignored[c("Pipeline", "Total")]), booktabs=TRUE, file = "Number_of_Species.tex")

# Reservoir 
ignored = dat[ , !(colnames(dat) %in% c("Pipeline", "Season"))]
ignored = aggregate(. ~ Reservoir, data=ignored, FUN=sum)
ignored$Total = rowSums(ignored[, !(colnames(ignored) %in% c("Reservoir"))] != 0)
ignored[c("Reservoir", "Total")]

# Season 
ignored = dat[ , !(colnames(dat) %in% c("Pipeline", "Reservoir"))]
ignored = aggregate(. ~ Season, data=ignored, FUN=sum)
ignored$Total = rowSums(ignored[, !(colnames(ignored) %in% c("Season"))] != 0)
ignored[c("Season", "Total")]

# Number of species for pipeline, reservoir, and season together
total = dat
total$total = rowSums(total[, !(colnames(total) %in% c("Pipeline", "Reservoir", "Season"))] != 0)
total = total[c("Pipeline", "Reservoir", "Season", "total")]
total
total[total$total == min(total$total),]
total[total$total == max(total$total),]

### Average number of reads per species

number_of_reads = colSums(dat[ , !(colnames(dat) %in% c("Reservoir", "Season", "Pipeline"))])
# Where 1 means to apply FUN to each row of df, 2 would mean to apply FUN to columns.
# min_of_reads = apply(dat[ , !(colnames(dat) %in% c("Reservoir", "Season", "Pipeline"))], 2, FUN=min)
st_dev = sapply(dat[ , !(colnames(dat) %in% c("Reservoir", "Season", "Pipeline"))], sd)
number_of_reads = round(number_of_reads / 5)
number_of_reads = data.frame(number_of_reads, st_dev)
number_of_reads$Species = rownames(number_of_reads)
number_of_reads[order(number_of_reads$number_of_reads),]

maximum = max(number_of_reads$number_of_reads)

# Create plots
p = ggplot(number_of_reads, aes(x = factor(Species, level = species_order), y = number_of_reads)) +
  geom_errorbar(aes(ymin = number_of_reads-st_dev, ymax = number_of_reads+st_dev), alpha = 0.75) + 
  geom_point(aes(fill = number_of_reads), shape = 21, size = 2.5) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_text(angle = 0, size = 8, color = "black",
                                   vjust = 1, hjust = 0.5),
        axis.text.y = element_text(size = 8, color = "black", face = "italic"),
        axis.title = element_text(size = 10, face = "plain"),
        plot.margin = margin(10, 10, 10, 15),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "lightgray", size = 0.3),
        panel.grid.minor = element_line(colour = "lightgray", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.position = "none"
  ) + coord_flip() +
  scale_fill_gradientn(limits = c(1, maximum),
                       expand = c(0, 0),
                       colors = rev(c("darkred", "red", "orange", "yellow", "green", "lightgreen", "lightblue", "darkblue")),
                       name = NULL)

plog = ggplot(number_of_reads, aes(x = factor(Species, level = species_order), y = number_of_reads)) +
  geom_errorbar(aes(ymin = number_of_reads-st_dev, ymax = number_of_reads+st_dev), alpha = 0.75) +
  geom_point(aes(fill = number_of_reads), shape = 21, size = 2.5) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_text(angle = 0, size = 8, color = "black",
                                   vjust = 1, hjust = 0.5),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 10, face = "plain"),
        plot.margin = margin(10, 10, 10, 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "lightgray", size = 0.3),
        panel.grid.minor = element_line(colour = "lightgray", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.position = "none"
  ) + scale_y_continuous(trans='log2') + coord_flip() +
  scale_fill_gradientn(limits = c(1, maximum),
                       expand = c(0, 0),
                       colors=rev(c("darkred", "red", "orange", "yellow", "green", "lightgreen", "lightblue", "darkblue")),
                       name = NULL)

# Add label on the right side of the plots
p = grid.arrange(p, top = text_grob("No Transformation", hjust = -0.25, vjust = 1.5, just = "centre", size = 8, face = "plain"))
plog = grid.arrange(plog, top = text_grob("Logarithm 2", hjust = 0.5, vjust = 1.5, just = "centre", size = 8, face = "plain"))

# Join plots
figure = ggarrange(p, plog, ncol = 2, nrow = 1, widths=c(1.65,1))
#figure = annotate_figure(figure, bottom = text_grob("Number of reads", hjust = 0, size = 10, face = "plain"))

# Export ggplot to Latex
# factor(Species, level = species_order) to change the order based on phylogenetic tree
tikz(file = "Number_of_reads_average_by_species.tex", width = 6, height = 5.5)
figure
dev.off()


### Number of reads per species considering Pipeline
ignored = dat[ , !(colnames(dat) %in% c("Reservoir", "Season"))]
ignored = aggregate(. ~ Pipeline, data=ignored, FUN=sum)

# For each pipeline
rownames(ignored) = ignored$Pipeline
pipelines = unique(ignored$Pipeline)
#ignored$Pipeline = NULL
for (pipeline in pipelines) {
  print(sort(ignored[pipeline, !(colnames(ignored) %in% "Pipeline")]))
  print("################")
}


# Export ggplot to Latex

# Convert to format of GGplot
number_of_reads = melt(ignored, id.vars='Pipeline')

# Calculate breaks
maximum = max(number_of_reads$value)
breaks = c(maximum)
while (tail(breaks, n=1) != 0) {
  breaks = c(breaks, round(tail(breaks, n=1) / 4))
}
breaks = breaks[-length(breaks)]

# Convert zero to NA
number_of_reads[number_of_reads == 0] = NA

tikz(file = "Number_of_reads_by_pipelines_and_species.tex", width = 6.25, height = 5.5)
ggplot(number_of_reads, aes(y=Pipeline, x=factor(variable, level = species_order), fill=value)) + 
  geom_tile(color = "black", size = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black", face = "italic"),
        axis.title = element_text(size = 10, face = "plain"),
        axis.ticks.x = element_blank(),
        plot.margin = margin(10, 10, 10, 20),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'lightgray'),
        panel.grid.major = element_line(colour = "white", size = 0.3),
        panel.grid.minor = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key.size = unit(1, "cm"),
        legend.position = "right",
        legend.text = element_text(angle = 0, size = 6, color = "black",
                                   face = "plain", vjust = 1, hjust = 1),
        strip.text = element_text(size = 8, face = "plain", color = "black"),
        strip.background = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(limits=c(1,maximum),
                      breaks = breaks,
                      expand = c(0,0),
                      colors=rev(c("darkred", "red", "orange", "yellow", "green", "lightgreen", "lightblue", "darkblue")),
                      na.value = 'white',
                      name = NULL,
                      trans = pseudo_log_trans(base = 2)) +
  coord_flip() + facet_grid(~Pipeline, scales='free')
dev.off()
# Because trans = is applying log2, so breaks values are 2**number (exponential of 2)


### Number of reads by species considering Reservoir
ignored = dat[ , !(colnames(dat) %in% c("Pipeline", "Season"))]
ignored = aggregate(. ~ Reservoir, data=ignored, FUN=sum)

# For each reservoir
rownames(ignored) = ignored$Reservoir
reservoirs = unique(ignored$Reservoir)
#ignored$Reservoir = NULL
for (reservoir in reservoirs) {
  print(sort(round(ignored[reservoir, !(colnames(ignored) %in% "Reservoir")] / length(unique(dat$Pipeline)))))
  print("################")
}


# Export ggplot to Latex

# Convert to format of GGplot
number_of_reads = melt(ignored, id.vars='Reservoir')
number_of_reads$value = round(number_of_reads$value / 5)

# Calculate breaks
maximum = max(number_of_reads$value)
breaks = c(maximum)
while (tail(breaks, n=1) != 0) {
  breaks = c(breaks, round(tail(breaks, n=1) / 4))
}
breaks = breaks[-length(breaks)]

# Convert zero to NA
number_of_reads[number_of_reads == 0] = NA

tikz(file = "Number_of_reads_by_reservoirs_and_species.tex", width = 6.25, height = 5.5)
ggplot(number_of_reads, aes(y=Reservoir, x=factor(variable, level = species_order), fill=value)) + 
  geom_tile(color = "black", size = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black", face = "italic"),
        axis.title = element_text(size = 10, face = "plain"),
        axis.ticks.x = element_blank(),
        plot.margin = margin(10, 10, 10, 20),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'lightgray'),
        panel.grid.major = element_line(colour = "white", size = 0.3),
        panel.grid.minor = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key.size = unit(1, "cm"),
        legend.position = "right",
        legend.text = element_text(angle = 0, size = 6, color = "black",
                                   face = "plain", vjust = 1, hjust = 1),
        strip.text = element_text(size = 8, face = "plain", color = "black"),
        strip.background = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(limits=c(1,maximum),
                       breaks = breaks,
                       expand = c(0,0),
                       colors=rev(c("darkred", "red", "orange", "yellow", "green", "lightgreen", "lightblue", "darkblue")),
                       na.value = 'white',
                       name = NULL,
                       trans = pseudo_log_trans(base = 2)) +
  coord_flip() + facet_grid(~Reservoir, scales='free')
dev.off()
# Because trans = is applying log2, so breaks values are 2**number (exponential of 2)


### Number of reads per species considering Season
ignored = dat[ , !(colnames(dat) %in% c("Reservoir", "Pipeline"))]
ignored = aggregate(. ~ Season, data=ignored, FUN=sum)

# For each Season
rownames(ignored) = ignored$Season
seasons = unique(ignored$Season)
#ignored$Season = NULL
for (season in seasons) {
  print(sort(round(ignored[season, !(colnames(ignored) %in% "Season")] / length(unique(dat$Pipeline)))))
  print("################")
}


# Export ggplot to Latex

# Convert to format of GGplot
number_of_reads = melt(ignored, id.vars='Season')
number_of_reads$value = round(number_of_reads$value / 5)

# Calculate breaks
maximum = max(number_of_reads$value)
breaks = c(maximum)
while (tail(breaks, n=1) != 0) {
  breaks = c(breaks, round(tail(breaks, n=1) / 4))
}
breaks = breaks[-length(breaks)]

# Convert zero to NA
number_of_reads[number_of_reads == 0] = NA

tikz(file = "Number_of_reads_by_seasons_and_species.tex", width = 6.25, height = 5.5)
ggplot(number_of_reads, aes(y=Season, x=factor(variable, level = species_order), fill=value)) + 
  geom_tile(color = "black", size = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black", face = "italic"),
        axis.title = element_text(size = 10, face = "plain"),
        axis.ticks.x = element_blank(),
        plot.margin = margin(10, 10, 10, 20),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'lightgray'),
        panel.grid.major = element_line(colour = "white", size = 0.3),
        panel.grid.minor = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key.size = unit(1, "cm"),
        legend.position = "right",
        legend.text = element_text(angle = 0, size = 6, color = "black",
                                   face = "plain", vjust = 1, hjust = 1),
        strip.text = element_text(size = 8, face = "plain", color = "black"),
        strip.background = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(limits=c(1,maximum),
                       breaks = breaks,
                       expand = c(0,0),
                       colors=rev(c("darkred", "red", "orange", "yellow", "green", "lightgreen", "lightblue", "darkblue")),
                       na.value = 'white',
                       name = NULL,
                       trans = pseudo_log_trans(base = 2)) +
  coord_flip() + facet_grid(~Season, scales='free')
dev.off()
# Because trans = is applying log2, so breaks values are 2**number (exponential of 2)


### Number of reads by species considering pipeline, reservoir, Season
df = dat

# Remove columns and create a new column with row names
df$Pipeline = NULL
df$Reservoir = NULL
df$Season = NULL
df$Pipeline_Reservoir_Season = rownames(df)

# Export ggplot to Latex

# Convert to format of GGplot
number_of_reads = melt(df, id.vars='Pipeline_Reservoir_Season')

# Sort and check in the file created the largest and smallest value
write.csv(number_of_reads[order(number_of_reads$value),], 
          file = "LS.csv", row.names = FALSE)

# Calculate breaks
maximum = max(number_of_reads$value)
breaks = c(maximum)
while (tail(breaks, n=1) != 0) {
  breaks = c(breaks, round(tail(breaks, n=1) / 4))
}
breaks = breaks[-length(breaks)]

# Convert zero to NA
number_of_reads[number_of_reads == 0] = NA
number_of_reads = number_of_reads %>% separate(Pipeline_Reservoir_Season, c("Pipeline", "Reservoir", "Season"), " ")

tikz(file = "Number_of_reads_by_PRS_and_species.tex", width = 6, height = 9.5)
ggplot(number_of_reads, aes(y=Season, x=factor(variable, level = species_order), fill=value)) + 
  geom_tile(color = "black", size = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, size = 8, color = "black",
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8, color = "black", face = "italic"),
        axis.title = element_text(size = 10, face = "plain"),
        plot.margin = margin(10, 10, 10, 20),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold",
                                  margin = margin(10, 0, 10, 0)),
        panel.background = element_rect(fill = 'lightgray'),
        panel.grid.major = element_line(colour = "white", size = 0.3),
        panel.grid.minor = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key.size = unit(1, "cm"),
        legend.position = "right",
        legend.text = element_text(angle = 0, size = 6, color = "black",
                                   face = "plain", vjust = 1, hjust = 1),
        strip.text = element_text(size = 8, face = "plain", color = "black")
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(limits=c(1,maximum),
                       breaks = breaks,
                       expand = c(0,0),
                       colors=rev(c("darkred", "red", "orange", "yellow", "green", "lightgreen", "lightblue", "darkblue")),
                       na.value = 'white',
                       name = NULL,
                       trans = pseudo_log_trans(base = 2)) +
  coord_flip() + facet_grid(Reservoir~Pipeline)
dev.off()
# Because trans = is applying log2, so breaks values are 2**number (exponential of 2)


### Species detected only once or in all pipelines

# Get data and ignore Reservoir and Season columns
ignored = dat[ , !(colnames(dat) %in% c("Reservoir", "Season"))]
# Sum up values by Pipeline
ignored = aggregate(. ~ Pipeline, data=ignored, FUN=sum)
# Make rownames the values in Pipeline column
rownames(ignored) = ignored$Pipeline
# Check if values are different from zero
ignored = ignored[, !(colnames(ignored) %in% c("Pipeline"))] != 0

# keep only columns summing up 1
names_only_once = colnames(ignored[,colSums(ignored) == 1])
# keep only columns summing up 5 (all pipelines)
names_in_all = colnames(ignored[,colSums(ignored) == 5])

# Show which pipeline detected species found just once
rowSums(ignored[,names_only_once])


### Species detected only once or in all reservoirs

# Get data and ignore Pipeline and Season columns
ignored = dat[ , !(colnames(dat) %in% c("Pipeline", "Season"))]
# Sum up values by Reservoir
ignored = aggregate(. ~ Reservoir, data=ignored, FUN=sum)
# Make rownames the values in Reservoir column
rownames(ignored) = ignored$Reservoir
# Check if values are different from zero
ignored = ignored[, !(colnames(ignored) %in% c("Reservoir"))] != 0

# keep only columns summing up 1
names_only_once = colnames(ignored[,colSums(ignored) == 1])
# keep only columns summing up 3 (all reservoirs)
names_in_all = colnames(ignored[,colSums(ignored) == 3])

# Show which reservoir detected species found just once
rowSums(ignored[,names_only_once])


### Species detected only once or in all seasons

# Get data and ignore Pipeline and Reservoir columns
ignored = dat[ , !(colnames(dat) %in% c("Pipeline", "Reservoir"))]
# Sum up values by Season
ignored = aggregate(. ~ Season, data=ignored, FUN=sum)
# Make rownames the values in Season column
rownames(ignored) = ignored$Season
# Check if values are different from zero
ignored = ignored[, !(colnames(ignored) %in% c("Season"))] != 0

# keep only columns summing up 1
names_only_once = colnames(ignored[,colSums(ignored) == 1])
# keep only columns summing up 2 (all seasons)
names_in_all = colnames(ignored[,colSums(ignored) == 2])

# Show which season detected species found just once
rowSums(ignored[,names_only_once])


### Species detected only once or in pipelines, reservoirs and seasons together

df = dat

# Remove columns
df$Pipeline = NULL
df$Reservoir = NULL
df$Season = NULL

df = df != 0

# keep only columns summing up 1
names_only_once = colnames(df[,colSums(df) == 1])
# keep only columns summing up 30 (all combination)
names_in_all = colnames(df[,colSums(df) == 30])

# Show which pipeline,reservoir and season detected species found just once
rowSums(df[,names_only_once])
