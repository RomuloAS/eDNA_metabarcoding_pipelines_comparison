library(ranacapa)
library(microbiome)

# Use ranacapa app to get the sample_size value
#ranacapa::runRanacapaApp()

ranacapa <- read.csv('/path/ranacapa.txt', sep = "\t")
ranacapa_metadata <- read.csv('/path/ranacapa_metadata.txt', sep = "\t")
physeq_object <- convert_anacapa_to_phyloseq(ranacapa, ranacapa_metadata)
rarefied = custom_rarefaction(physeq_object, sample_size = 113122, replicates = 20)

otus <- otu_table(rarefied)
write.csv(otus, file='rarefied.csv')
