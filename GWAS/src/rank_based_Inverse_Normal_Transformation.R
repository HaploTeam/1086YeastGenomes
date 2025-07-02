#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: rank_based_Inverse_Normal_Transformation.r
##
## Purpose of script: Perform the rank-based INT from a phenotype file
##
## Author: Victor Loegler
##
## Date Created: 2021-10-11
##
## ---------------------------
##
## Notes: Takes as argument the phenotypes in .tsv format (but .phen extension), 
##        that is a tsv file with for each row:
##        Strain      Condition
##        XXX         1563
##        XXX         4563
##        ....
##
##        The pheno file can contain NAs, which will be kept in the output file. 
##        The output file is a tsv file containing 3 columns:
##        Strain  Strain  Normalized value
##
## ---------------------------
library(stringr)
## ---------------------------

rank.based.INT <- function(x, c=3/8, method="average")
{
  # This function performs the rank-based inverse normal transformation (INT)
  # If method is "average" ties will share the same average value. 
  # If method is "random", ties are given rank randomly
  # Formula found in Beasley 2009 with an offset of c=3/8 as recommended in Blom 1958
  
  r <- rank(x, ties.method = method)
  r[is.na(x)] <- NA # reput NA values in the vector because rank() gives a rank to NAs
  N <- length(x[!is.na(x)])
  qnorm((r-c)/(N-2*c+1))
}


# Get paths of input files
args = commandArgs(trailingOnly=TRUE)
phenotypesPath <- args[1]

# Get Cond name
cond <- unlist(strsplit(x = phenotypesPath, split = "/"))
cond <- cond[length(cond)] # get file name (last position of path)
cond <- str_remove(cond, ".phen") # remove .phen

# Convert the phenotypes to a single vector
phenotypes <- read.csv(file = phenotypesPath, header = T, sep = "\t")[,2]

# Get the strains list in a vector
strainList <- read.csv(file = phenotypesPath, header = T, sep = "\t")[,1]

# Normalization with the rank-based INT function
normalized.phenotypes <- rank.based.INT(x = phenotypes)

# Converting to dataframe
normalized.phenotypes.df <- data.frame(strainList, normalized.phenotypes)
colnames(normalized.phenotypes.df) <- c("Strain", cond)

# removing NA values for the final file
normalized.phenotypes.df.noNA <- na.omit(normalized.phenotypes.df)

# Exporting to tsv file
newPath <- paste(substr(phenotypesPath, start = 1, stop = nchar(phenotypesPath)-5), ".norm.phen", sep = "")
write.table(normalized.phenotypes.df.noNA , file = newPath, row.names = F, col.names = T, sep = "\t", quote = F)


