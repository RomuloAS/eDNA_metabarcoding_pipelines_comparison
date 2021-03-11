source("Create_Tables_Utils.R")
library(rChoiceDialogs)

# Number of the folder
folder.number = 3

# Folder variables
TABLES.JOINED = paste(auto.increment(folder.number), "_Joined", sep = "")
TABLES.JOINED.REDUCED = paste(auto.increment(folder.number), "_Joined_Reduced", sep = "")
#TABLES.RANACAPA = paste(auto.increment(folder.number), "_Ranacapa", sep = "")
TABLES.CONTROLS = paste(auto.increment(folder.number), "_No_Controls", sep = "")
TABLES.ONLY.CONTROLS = paste(auto.increment(folder.number), "_Only_Controls", sep = "")
TABLES.DETAILED = paste(auto.increment(folder.number), "_Detailed", sep = "")
TABLES.CUMULATIVE.SPECIES = paste(auto.increment(folder.number), "_Cumulative_Species", sep = "")
TABLES.CUMULATIVE.READS = paste(auto.increment(folder.number), "_Cumulative_Reads", sep = "")

# Set working directory to source file location
if(Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}else{
  setwd(utils::getSrcDirectory()[1])
}

# Read files
# Get table file names
tables.path = list.files(rchoose.dir(caption = "Tables directory"), pattern = "*.tsv", full.names = TRUE ,recursive = TRUE)

# Read tables
tables = read.tables(tables.path)

# Create folder to save the new file
dir.create(TABLES.JOINED, showWarnings = FALSE)
# Create table joined
create.table.joined(tables)

# Create folder to save the new files after removing controls
dir.create(TABLES.ONLY.CONTROLS, showWarnings = FALSE)
# Keep controls
# for (table in tables) {
#   keep.only.controls(table)
# }
# Keep only positive control
create.table.positive.control()

# Create folder to save the new file
dir.create(TABLES.JOINED.REDUCED, showWarnings = FALSE)
# Create table joined reduced
create.table.joined.reduced(tables)

# Before the execution of ranacapa code below, modify the line:
# new.name = paste(info, collapse = " ")
# by
# new.name = paste(info, collapse = ".")

# Create folder to save the new file
#dir.create(TABLES.RANACAPA, showWarnings = FALSE)
# Create tables ranacapa
#create.tables.ranacapa(all.pipelines.reduced.no.controls, "ranacapa", T)

# Create folder to save the new files after removing controls
dir.create(TABLES.CONTROLS, showWarnings = FALSE)
# Remove controls
for (table in tables) {
  remove.controls(table)
}

# Create folder to save the new files
dir.create(TABLES.DETAILED, showWarnings = FALSE)
# Create detailed version of the tables
for (table.name in tables) {
  create.table.detailed(table.name)
}

# Create folder to save the new file
dir.create(TABLES.CUMULATIVE.SPECIES, showWarnings = FALSE)
# Create tables cumulative species
create.table.cumulative.species(tables)

# Create folder to save the new file
dir.create(TABLES.CUMULATIVE.READS, showWarnings = FALSE)
# Create tables cumulative reads
create.table.cumulative.reads(tables)

# Create table all pipelines detailed
create.tables.all.pipelines.detailed()

# Save R object
save.image("Tables.RData")
