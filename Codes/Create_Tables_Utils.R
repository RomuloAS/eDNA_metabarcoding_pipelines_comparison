library(docstring)
library(data.table)
library(tidyverse)
library(collections)
library(taxize)
library(zeallot)

auto.increment = function(variable){
  #' "Auto" increment 1 to the variable
  #' 
  #' @description This function adds 1 to the variable.
  #' 
  #' @param variable numeric.
  #' @usage auto.increment(variable)
  #' @return It returns the variable after adding 1 to it.
  #' @details The input must be a numeric.
  #' @examples
  #' auto.increment(variable)
  
  # Get name of the variable
  name = deparse(substitute(variable))
  
  # Get value of the variable
  value = get(name)
  
  # Sum 1 to the value
  value = value + 1
  
  # Assign new value to a global variable with the same name
  # Auto-increment
  assign(name, value, envir = .GlobalEnv)
  
  # Return new value
  return(value)
}

modify.column.names <- function(column.names){
  #' Modify column names
  #' 
  #' @description This function modifies the column names.
  #' 
  #' It adds the name of the pipeline at the end, and it formats the 
  #' reservoir and season names.
  #' 
  #' @param column.names list
  #' @usage modify.column.names(column.names)
  #' @return It returns the new name.
  #' @details The input must be a list.
  #' @examples
  #' modify.column.name(column.names)
  
  # Create a dictionary for reservoirs
  reservoirs <- dict()
  reservoirs$set("K", "Klíčava")
  reservoirs$set("R", "Římov")
  reservoirs$set("Z", "Žlutice")
  
  # Create a dictionary for seasons
  seasons <- dict()
  seasons$set("S", "Summer")
  seasons$set("W", "Autumn")
  
  new.column.names <- lapply(column.names, function(column.name) {
    # Get reservoir
    reservoir = reservoirs$get(toupper(substring(column.name, 1, 1)))
    
    # Get season
    season = seasons$get(toupper(substring(column.name, 2, 2)))
    
    # Get extra info
    extra.info = substring(column.name, 3)
    
    # Join all info in a vector
    info = c(reservoir, season, extra.info, table.name)
    
    # Remove empty element in the vector
    info = info[info != ""]
    
    # Join reservoir, season, extra info, and pipeline name.
    new.name = paste(info, collapse = " ")
    
    # Return the new name
    return(new.name)
  })
  
  # Return new column names
  return(new.column.names)
}

create.tables.all.pipelines.detailed <- function(){
  #' Create Table All Pipelines Detailed
  #' 
  #' @description This function creates a detailed version of the
  #' table using all pipelines.
  #' 
  #' @usage create.tables.all.pipelines.detailed()
  #' @return It does not return anything.
  #' @examples
  #' create.tables.all.pipelines.detailed()
  
  # Remove TAX column and rows Total, Assigned and Unassigned
  table = all.pipelines.reduced.no.controls
  table$TAX = NULL
  table = table[!(row.names(table) %in% c("Total", "Assigned", "Unassigned")), ]
  table = as.data.frame(t(table))
  
  # Transform the name of the rows as the first column
  table = tibble::rownames_to_column(table, "Reservoir_Season_Pipeline")
  rownames(table) = table$Reservoir_Season_Pipeline 
  table = table %>% separate(Reservoir_Season_Pipeline, c("Reservoir", "Season", "Pipeline"), " ")
  
  # Assign table name
  assign("all.pipelines.detailed", table, envir = .GlobalEnv)
  
  # Write tsv file
  file = paste(TABLES.DETAILED, "/All_Pipelines_detailed.csv", sep = "")
  fwrite(x=table, file=file, sep = "\t", row.names = T)
  
}

create.tables.ranacapa.metadata <- function(table, table.name, reduced){
  #' Create Tables Ranacapa Metadata
  #' 
  #' @description This function creates a metadata table
  #' to be used with Ranacapa statistical analyses.
  #' 
  #' @param table dataframe.
  #' @param table.name name for the new table.
  #' @param reduced True if it is reduced or False if not.
  #' @usage create.tables.ranacapa.metadata(table, table.name, reduced)
  #' @return It does not return anything.
  #' @examples
  #' create.tables.ranacapa.metadata(ranacapa, "metadata", T)
  
  # Create new dataframes
  columns <- c("Sample", "Sample_or_Control", "Reservoir", "Season", "Pipeline")
  ###
  metadata <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(metadata) <- columns
  
  # Get the name of the columns
  columns = colnames(table)[colnames(table) != "sum.taxonomy"]
  
  # For each column name
  # Get reservoir, season, pipeline and sample_or_control
  for (column in columns) {
    
    column.info = str_split(column, "\\.", n = Inf, simplify = FALSE)
    
    # Extract info
    reservoir = column.info[[1]][1]
    season = column.info[[1]][2]
    info = column.info[[1]][3]
    pipeline = column.info[[1]][4]
    
    # Check if it is sample or control
    sample_or_control = "control"
    if (str_detect(info, ".*[0-9].*")) {
      sample_or_control = "sample"
    }
    
    if(reduced){
      pipeline = column.info[[1]][3]
      sample_or_control = "sample"
    }
    
    # Add row to the dataframe
    metadata = rbind(metadata, list(Sample = column,
                                    Sample_or_Control = sample_or_control,
                                    Reservoir = reservoir,
                                    Season = season,
                                    Pipeline = pipeline),
                     stringsAsFactors = FALSE)
    
  }
  
  # Write tsv file
  file = paste(TABLES.RANACAPA, "/", table.name, ".txt", sep = "")
  fwrite(x=metadata, file=file, sep = "\t", row.names = F)
}

create.tables.ranacapa <- function(table, table.name, reduced){
  #' Create Tables Ranacapa
  #' 
  #' @description This function creates a table to be used with
  #' Ranacapa statistical analyses.
  #' 
  #' @param table dataframe.
  #' @param table.name name for the new table.
  #' @param reduced True if it is reduced or False if not.
  #' @usage create.tables.ranacapa(table, table.name, reduced)
  #' @return It does not return anything.
  #' @examples
  #' create.tables.ranacapa(all.pipelines, "my_new_table", T)
  
  # Create initial dataframes
  ranacapa = table[!(row.names(table) %in% c("Total", "Assigned", "Unassigned")), ]
  
  # Get row names
  tax = rownames(ranacapa)
  
  # Change column name
  colnames(ranacapa)[which(names(ranacapa) == "TAX")] <- "sum.taxonomy"
  
  # Get taxonomic information
  tax.info = classification(tax, db='gbif', rows = 1)
  
  # Change tax from only species by phylum to species
  for (species in ranacapa$sum.taxonomy) {
    info = tax.info[[species]]
    phylum.to.species = paste(info[info$rank == "phylum",]$name,
                              info[info$rank == "class",]$name,
                              info[info$rank == "order",]$name,
                              info[info$rank == "family",]$name,
                              info[info$rank == "genus",]$name,
                              info[info$rank == "species",]$name,
                              sep = ";")
    
    ranacapa$sum.taxonomy[ranacapa$sum.taxonomy == species] = phylum.to.species
  }
  
  # Create metadata
  file.name = paste(table.name, "_metadata", sep = "")
  create.tables.ranacapa.metadata(ranacapa, file.name, reduced)
  
  # Write tsv file
  file = paste(TABLES.RANACAPA, "/", table.name, ".txt", sep = "")
  fwrite(x=ranacapa, file=file, sep = "\t", row.names = F)
  
}

create.table.joined.reduced <- function(tables){
  #' Create Table Joined Reduced for all Pipelines
  #' 
  #' @description This function creates a table joining all pipelines
  #' in a unique table for each tuple reservoir and season. It also 
  #' removes controls.
  #' 
  #' @param tables vector with name each table variable.
  #' @usage create.table.joined.reduced(tables)
  #' @return It does not return anything.
  #' @details The input must be a vector with the names
  #' of the tables.
  #' @examples
  #' create.table.joined.reduced(tables)
  
  # Create new dataframe
  columns <- c("TAX")
  all.pipelines.reduced.no.controls <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(all.pipelines.reduced.no.controls) <- columns
  
  # For each table
  for (table.name in tables) {
    
    # Assign table name
    assign("table.name", table.name, envir = .GlobalEnv)
    
    # Get table
    table = get(table.name)
    
    # Remove control Mayland zebra
    table = table[!(row.names(table) %in% "Maylandia zebra"),]
    
    # Collect control names
    control.names = colnames(table[,! grepl("\\d", colnames(table))])
    
    # Remove control columns
    table = table[,!(colnames(table) %in% control.names)]
    
    # Create vector with the combination of the reservoir and the season
    reservoirs.seasons = c("KS", "KW", "RS", "RW", "ZS", "ZW")
    
    for (reservoir.season in reservoirs.seasons){
      
      # Sum up columns matching reservoir.season variable pattern
      row.sums = rowSums(table[ , grepl(reservoir.season, colnames(table))])
      
      # Remove columns matching reservoir.season variable pattern
      table = table[ , ! grepl(reservoir.season, colnames(table))]
      
      # Add column where the name is the pattern from reservoir.season variable
      table[, reservoir.season] = row.sums
    }
    
    # Rename table columns
    table = table %>% dplyr::rename_all(modify.column.names)
    
    # Transform the name of the rows as the first column
    table = tibble::rownames_to_column(table, "TAX")
    
    # Merge table with the dataframe
    all.pipelines.reduced.no.controls = merge(all.pipelines.reduced.no.controls,
                                              table, by="TAX", all=TRUE)
  }
  
  # Set NA as zero
  all.pipelines.reduced.no.controls[is.na(all.pipelines.reduced.no.controls)] <- 0
  
  # Set indexes as tax
  rownames(all.pipelines.reduced.no.controls) <- all.pipelines.reduced.no.controls$TAX
  
  # Get row names
  tax = rownames(all.pipelines.reduced.no.controls)
  
  # Put Total, Assigned and Unassigned to the end
  tax = tax[which(!tax %in% c("Total", "Assigned", "Unassigned"))]
  tax = c(tax, c("Total", "Assigned", "Unassigned"))
  all.pipelines.reduced.no.controls = all.pipelines.reduced.no.controls[tax,]
  
  # Sort columns and keep TAX at the beginning
  columns = sort(colnames(all.pipelines.reduced.no.controls))
  columns = columns[which(!columns %in% "TAX")]
  columns = c("TAX", columns)
  all.pipelines.reduced.no.controls = all.pipelines.reduced.no.controls[ , columns]
  
  # Assign to a global variable
  assign("all.pipelines.reduced.no.controls", all.pipelines.reduced.no.controls,
         envir = .GlobalEnv)
  
  # Write tsv file
  fwrite(x=all.pipelines.reduced.no.controls, file=paste(TABLES.JOINED.REDUCED,
                                                         "/all_pipelines_reduced_no_control.tsv",
                                                         sep = ""),
         sep = "\t",
         row.names = F)
  
}

create.table.joined <- function(tables) {
  #' Create Table Joined for all Pipelines
  #' 
  #' @description This function creates a table joining all pipelines
  #' in a unique table.
  #' 
  #' @param tables vector with name each table variable.
  #' @usage create.table.joined(tables)
  #' @return It does not return anything.
  #' @details The input must be a vector with the names
  #' of the tables.
  #' @examples
  #' create.table.joined(tables)
  
  # Create new dataframe
  columns <- c("TAX")
  all.pipelines <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(all.pipelines) <- columns
  
  # For each table
  for (table.name in tables) {
    
    # Assign table name
    assign("table.name", table.name, envir = .GlobalEnv)
    
    # Get table
    table = get(table.name)
    
    # Rename table columns
    table = table %>% dplyr::rename_all(modify.column.names)
    
    # Transform the name of the rows as the first column
    table = tibble::rownames_to_column(table, "TAX")
    
    # Merge table with the dataframe
    all.pipelines = merge(all.pipelines, table, by="TAX", all=TRUE)
  }
  
  # Set NA as zero
  all.pipelines[is.na(all.pipelines)] <- 0
  
  # Set indexes as tax
  rownames(all.pipelines) <- all.pipelines$TAX
  
  # Get row names
  tax = rownames(all.pipelines)
  
  # Put Total, Assigned and Unassigned to the end
  tax = tax[which(!tax %in% c("Total", "Assigned", "Unassigned"))]
  tax = c(tax, c("Total", "Assigned", "Unassigned"))
  all.pipelines = all.pipelines[tax,]
  
  # Sort columns and keep TAX at the beginning
  columns = sort(colnames(all.pipelines))
  columns = columns[which(!columns %in% "TAX")]
  columns = c("TAX", columns)
  all.pipelines = all.pipelines[ , columns]
  
  # Assign to a global variable
  assign("all.pipelines", all.pipelines, envir = .GlobalEnv)
  
  # Write tsv file
  fwrite(x=all.pipelines, file=paste(TABLES.JOINED,
                                     "/all_pipelines.tsv",
                                     sep = ""),
         sep = "\t",
         row.names = F)
  
}

create.table.cumulative.reads <- function(tables){
  #' Create Table Cumulative Reads
  #' 
  #' @description This function creates a table showing the
  #' cumulative summary by reads.
  #'
  #' It counts the number of reads for each tuple reservoir and season.
  #' 
  #' @param tables vector with name each table variable.
  #' @usage create.table.cumulative.reads(tables)
  #' @return It does not return anything.
  #' @details The input must be a dataframe table.
  #' @examples
  #' create.table.cumulative.reads(tables)
  
  # Create new dataframe
  columns <- c("Reservoir", "Season", "Assigned_Reads", "Pipeline")
  cumulative.reads <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(cumulative.reads) <- columns
  
  # For  each table
  for (table.name in tables) {
    
    # Get table
    table = get(table.name)
    
    # Count the number of reads detected for each pipeline
    total = rowSums(table[,! colnames(table) %in% c("Reservoir", "Season")])
    
    # Add row to the dataframe
    cumulative.reads = rbind(cumulative.reads, list(Reservoir = table[,"Reservoir"],
                                                    Season = table[,"Season"],
                                                    Assigned_Reads = total,
                                                    Pipeline = rep(table.name,
                                                                   times = 6,
                                                                   length.out = NA,
                                                                   each = 1)
    ),
    stringsAsFactors = FALSE)
    
  }
  
  # Assign to a global variable
  assign("cumulative.reads", cumulative.reads, envir = .GlobalEnv)
  
  # Reset indexes (row names)
  rownames(cumulative.reads) <- NULL
  
  # Write tsv file
  fwrite(x=cumulative.reads, file=paste(TABLES.CUMULATIVE.READS,
                                        "/","Cumulative_reads",
                                        ".tsv", sep = ""),
         sep = "\t",
         row.names = F)
  
}

create.table.cumulative.species <- function(tables){
  #' Create Table Cumulative Species
  #' 
  #' @description This function creates a table
  #' showing the cumulative summary by species.
  #'
  #' It counts the number of species for each tuple reservoir and season.
  #' 
  #' @param tables vector with name each table variable.
  #' @usage create.table.cumulative.species(tables)
  #' @return It does not return anything.
  #' @details The input must be a dataframe table.
  #' @examples
  #' create.table.cumulative.species(tables)
  
  # Create new dataframe
  columns <- c("Reservoir", "Season", "N_Species", "Pipeline")
  cumulative.species <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(cumulative.species) <- columns
  
  # For  each table
  for (table.name in tables) {
    
    # Get table
    table = get(table.name)
    
    # Convert to logical
    logical = table[,! colnames(table) %in% c("Reservoir", "Season")] %>% mutate_all(as.logical)
    
    # Count the number of species detected for each pipeline
    total = rowSums(logical %>% mutate_all(as.numeric))
    
    # Add row to the dataframe
    cumulative.species = rbind(cumulative.species, list(Reservoir = table[,"Reservoir"],
                                                        Season = table[,"Season"],
                                                        N_Species = total,
                                                        Pipeline = rep(table.name,
                                                                       times = 6,
                                                                       length.out = NA,
                                                                       each = 1)
    ),
    stringsAsFactors = FALSE)
    
  }
  
  # Assign to a global variable
  assign("cumulative.species", cumulative.species, envir = .GlobalEnv)
  
  # Reset indexes (row names)
  rownames(cumulative.species) <- NULL
  
  # Write tsv file
  fwrite(x=cumulative.species, file=paste(TABLES.CUMULATIVE.SPECIES,
                                          "/","Cumulative_species",
                                          ".tsv", sep = ""),
         sep = "\t",
         row.names = F)
  
}

create.table.detailed <- function(table.name){
  #' Create Tables detailed
  #' 
  #' @description This function creates a detailed version
  #' of the original one.
  #'
  #' It creates a table where reservoir, season and species are the columns.
  #' 
  #' @param table.name dataframe Where the indexes are species and 
  #' columns are reservoirs.
  #' @usage create.table.detailed(table.name)
  #' @return It does not return anything.
  #' @details The input must be a dataframe table.
  #' @examples
  #' create.table.detailed(table.name)
  
  # Get table without Unassigned, Assigned and Total
  table = get(table.name)[! row.names(get(table.name)) %in% c("Unassigned", "Assigned", "Total"),]
  
  # Create vector with the combination of the reservoir and the season
  reservoirs.seasons = c("KS", "KW", "RS", "RW", "ZS", "ZW")
  
  for (reservoir.season in reservoirs.seasons){
    
    # Sum up columns matching reservoir.season variable pattern
    row.sums = rowSums(table[ , grepl(reservoir.season, colnames(table))])
    
    # Remove columns matching reservoir.season variable pattern
    table = table[ , ! grepl(reservoir.season, colnames(table))]
    
    # Add column where the name is the pattern from reservoir.season variable
    table[, reservoir.season] = row.sums
  }
  
  # Transpose Dataframe
  table = as.data.frame(t(table))
  
  # Reorder dataframe by columns
  table = table[,order(colnames(table))]
  
  # Create new column season
  seasons = c("Summer", "Autumn")
  table = cbind(Season = seasons, table)
  
  # Create new column Reservoir
  reservoirs = c("Klíčava", "Klíčava", "Římov", "Římov", "Žlutice", "Žlutice")
  table = cbind(Reservoir = reservoirs, table)
  
  # Assign to a global variable
  assign(table.name, table, envir = .GlobalEnv)
  
  # Write tsv file
  fwrite(x=table, file=paste(TABLES.DETAILED, "/",table.name, ".tsv", sep = ""),
         sep = "\t", row.names = F)
  
}

keep.only.controls <- function(table){
  #' Keep only Controls
  #' 
  #' @description This function removes no controls from table
  #' and keep only controls.
  #'
  #' It removes not control and keep D, E, F, N, P, and de.
  #' 
  #' @param table dataframe where the indexes are species and 
  #' columns are reservoirs.
  #' @usage keep.only.controls(table)
  #' @return It does not return anything.
  #' @details The input must be a dataframe table.
  #' @examples
  #' keep.only.controls(table)
  
  # Collect control names
  control.names = colnames(get(table)[, grepl("\\d", colnames(get(table)))])
  
  # Remove control columns
  assign(table, get(table)[,!(colnames(get(table)) %in% control.names)], envir = .GlobalEnv)
  
  # Write tsv file
  fwrite(x=get(table), file=paste(TABLES.ONLY.CONTROLS, "/",table, ".tsv", sep = ""),
         sep = "\t", row.names = T)
}

create.table.positive.control <- function(){
  #' Create Table All Pipelines only positive control
  #' 
  #' @description This function creates a table with only
  #' positive control using all pipelines.
  #' 
  #' @usage create.table.positive.control()
  #' @return It does not return anything.
  #' @examples
  #' create.table.positive.control()
  
  # Remove TAX column and rows Total, Assigned and Unassigned
  table = all.pipelines
  table$TAX = NULL
  table = table[!(row.names(table) %in% c("Assigned", "Unassigned")), ]
  table = as.data.frame(t(table))
  
  # Keep only Mayland Zebra
  table = table[table["Maylandia zebra"] > 0, c("Total", "Maylandia zebra")]
  
  # Correction of row names
  rownames(table) = gsub(" P ", " ", rownames(table))
  
  # Calculate percentage
  table["Percentage_of_the_total"] = round(table["Maylandia zebra"] / table["Total"] * 100, digits = 2)
  
  # Transform the name of the rows as the first column
  table = tibble::rownames_to_column(table, "Reservoir_Season_Pipeline")
  rownames(table) = table$Reservoir_Season_Pipeline 
  table = table %>% separate(Reservoir_Season_Pipeline, c("Reservoir", "Season", "Pipeline"), " ")
  
  # Assign table name
  assign("all.pipelines.only.positive.control", table, envir = .GlobalEnv)
  
  # Write tsv file
  file = paste(TABLES.ONLY.CONTROLS, "/All_Pipelines_only_positive_control.csv", sep = "")
  fwrite(x=table, file=file, sep = "\t", row.names = T)
  
}

remove.controls <- function(table){
  #' Remove Controls
  #' 
  #' @description This function removes controls from table.
  #'
  #' It removes index Mayland zebra and columns D, E, F, N, P, and de.
  #' 
  #' @param table dataframe Where the indexes are species and 
  #' columns are reservoirs.
  #' @usage remove.controls(table)
  #' @return It does not return anything.
  #' @details The input must be a dataframe table.
  #' @examples
  #' remove.controls(table)
  
  # Remove Maylandia zebra
  assign(table, get(table)[!(row.names(get(table)) %in% "Maylandia zebra"),], envir = .GlobalEnv)
  
  # Collect control names
  control.names = colnames(get(table)[,! grepl("\\d", colnames(get(table)))])
  
  # Remove control columns
  assign(table, get(table)[,!(colnames(get(table)) %in% control.names)], envir = .GlobalEnv)
  
  # Write tsv file
  fwrite(x=get(table), file=paste(TABLES.CONTROLS, "/",table, ".tsv", sep = ""),
         sep = "\t", row.names = T)
}

read.tables <- function(tables.path){
  #' Read Tables
  #' 
  #' @description This function reads all tables from the vector
  #'
  #' From a vector of table paths, It reads each tsv file
  #' and creates a dataframe for each one. 
  #' 
  #' The name of each variable is the name of each file without
  #' extension.
  #' 
  #' @param tables.path vector. A vector with paths to tables.
  #' @usage read.tables(tables.path)
  #' @return Return a vector with the names of the dataframes
  #' for each table.
  #' @details The input must be a vector with the paths to the tables
  #' in tsv format.
  #' @examples
  #' read.tables(tables.path)
  #' read.tables("PATH_TO_TABLE")
  
  # Initialize tables variable
  tables = c()
  
  # For each file
  for (file.name in tables.path) {
    
    # Parse base name and dir name
    base.name = basename(file.name)
    dir.name = dirname(file.name)
    
    # File name without extension
    variable.name = tools::file_path_sans_ext(base.name)
    
    # Read table and assign it to the variable.name
    assign(variable.name, read.csv(file.name, sep = '\t', header=T, row.names = 1, check.names=F), envir = .GlobalEnv)
    
    # Add new element to the list of tables
    tables = c(tables, variable.name)
  }
  
  return(tables)
}