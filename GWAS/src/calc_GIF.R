args <- commandArgs(trailingOnly = TRUE)

fileIn <- paste(args[1])
#print(fileIn)
directory = getwd()
setwd(dirname(fileIn))
assoc_df <- read.table(basename(fileIn), header=TRUE, stringsAsFactors = FALSE, sep='\t')
assoc_df <- assoc_df[complete.cases(assoc_df[,6]),]
chisq <- qchisq(1-assoc_df$PValue,1)
#lambda GC
lgc <- median(chisq)/qchisq(0.5,1)
#print(lgc)

df<-data.frame(fileIn, lgc)

out <- paste(fileIn, ".lgc", sep="")
write.table(df, sep='\t', file=out, row.names=FALSE, quote = FALSE)

setwd(directory)

