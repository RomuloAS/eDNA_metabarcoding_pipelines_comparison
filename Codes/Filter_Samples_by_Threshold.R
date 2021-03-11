library(rChoiceDialogs)
library(data.table)

# Threshold of 0.1% 
# It needs to be divided by 100 in R
threshold = 0.001

# Get table file names
tables.path = list.files(rchoose.dir(caption = "Choose tables directory"), pattern = "*.tsv", full.names = TRUE ,recursive = TRUE)

for (file.name in tables.path) {

  # Parse base name and dir name
  base.name = basename(file.name)
  dir.name = dirname(file.name)

  # Read file
  dat = read.csv(file.name, sep = '\t', header=T, row.names = 1, check.names=F)
  # Remove last column "taxomomy" if you did not remove it#
  # dat = dat[, !colnames(dat) %in% "taxomomy",]

  # Transpose data
  dat = as.data.frame(t(dat))
  
  # Create 'datt': a proportion reads data frame#
  datt = dat / dat$Total
  datt[is.na(datt)] <- 0
  
  # Apply threshold to proportions of datt to dat#
  dat[datt < threshold] = 0
  
  # Calculate assigned 
  dat$Assigned = rowSums(dat[, !colnames(dat) %in% c("Unassigned", "Assigned", "Total")])
  
  # Calculate unassigned
  dat$Unassigned = dat$Total - dat$Assigned

  # Transpose back
  dat = as.data.frame(t(dat))

  # Remove species that sum zero
  dat = dat[rowSums(dat) > 0,]

  # Export dataframe to CSV
  file.result.name <- paste(dir.name, '/Filtered_', base.name, sep = "")
  fwrite(x=dat, file=file.result.name, sep = "\t", row.names = T)
}
  
