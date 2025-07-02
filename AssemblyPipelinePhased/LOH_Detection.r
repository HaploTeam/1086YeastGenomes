#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: LOH_Detection.r
##
## Purpose of script: Detect LOH region from SNP data (VCF file)
##                    A LOH region is defines as a 50kb region that contains
##                    less than 10 SNPs. 
##
## Author: Victor Loegler
##
## Date Created: 2022-10-10
##
## ---------------------------
##
## Notes: The script takes as argument
##        -b --bed Reference.bed
##        -v --vcf SNPs.vcf/.vcf.gz/.bcf
##        
##        Note that bcftools must be in the path to make this script work, and
##        that the VCF file must have a AC INFO (equivalent of AD for Longshot
##        SNP caller, can be replaced by AD). 
##
## ---------------------------
library(optparse)
library(ggplot2)
library(dplyr)
## ---------------------------

option_list = list(
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="BED file of the reference genome", metavar="character"),
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help="SNP matrix at the vcf, vcf.gz of bcf format", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

refBED <- read.table(file = opt$bed, header = F, sep = "\t")
colnames(refBED) <- c("Chr", "Start", "End")
refChr <- refBED$Chr
refLen <- refBED$End

prefix <- sub(".vcf", "", sub(".vcf.gz", "", sub(".bcf", "", basename(opt$vcf))))

# Filter VCF to keep only biallelic heterozygous SNPs
# Any other filtering (quality, depth) must be done before using this script
system(command = paste0(
  "bcftools view -g het -v snps -m2 -M2 ", opt$vcf, " | bcftools query -f \"%CHROM\t%POS\t%AC\n\" > ", prefix, "_HeterozygousSNPs_tmp.tsv"
))

SNPdata <- read.table(file = paste0(prefix, "_HeterozygousSNPs_tmp.tsv"), header = F, sep = "\t")
colnames(SNPdata) <- c("Chr", "Pos", "AC")
SNPdata$AB <- sapply(SNPdata$AC, function(AC) as.numeric(unlist(strsplit(AC, ",")))[1] / sum(as.numeric(unlist(strsplit(AC, ",")))))
SNPdata$PosKb <- ceiling(SNPdata$Pos / 1000)

# Get number of SNP in each 1kb window
SNPdistribution <- SNPdata %>%
  group_by(Chr,PosKb) %>%
  summarise(HetSNPs = n())

# Get each 50kb window with less than 10 SNPs in total
LOH <- data.frame()
for (chr in refChr){
  print(paste("Running on chr:", chr))
  for (pos in c(1:(ceiling(refLen[match(chr, refChr)]/1000)-49))){
    nbSNPs <- sum(SNPdistribution$HetSNPs[SNPdistribution$Chr == chr & SNPdistribution$PosKb %in% c(pos:(pos+49))])
    if (nbSNPs < 10){
      start = 1+(1000*(pos-1))
      end = 1000*(pos + 49)
      if (end > refLen[match(chr, refChr)]) {end <- refLen[match(chr, refChr)]} # If end is larger than chr length
      LOH <- rbind(LOH, data.frame(Chr = chr, Start = start, End = end))
    }
  }
}

if (nrow(LOH > 0)){
  # Merge overlapping 50kb LOH regions
  mergeLOH <- data.frame()
  for (chr in unique(sort(LOH$Chr))){
    tmpLOH <- LOH[LOH$Chr == chr,]
    start = tmpLOH$Start[1]
    end = tmpLOH$End[1]
    if (nrow(tmpLOH) == 1) { # If only a single LOH region in the chromosome
      mergeLOH <- rbind(mergeLOH, data.frame(Chr = chr, Start = start, End = end))
    } else {
      for (i in c(2:nrow(tmpLOH))){
        if (tmpLOH$Start[i] <= tmpLOH$End[i-1]){
          end = tmpLOH$End[i]
          if (i == nrow(tmpLOH)) {
            mergeLOH <- rbind(mergeLOH, data.frame(Chr = chr, Start = start, End = end))
          }
        } else {
          mergeLOH <- rbind(mergeLOH, data.frame(Chr = chr, Start = start, End = end))
          start = tmpLOH$Start[i]
          end = tmpLOH$End[i]
        }
      }
    }
  }
  # Order by chromosome
  LOH <- mergeLOH
  LOH$Chr <- factor(LOH$Chr, levels = refChr)
  LOH <- LOH[order(LOH$Chr, LOH$Start),]
  # Convert 1-based inclusive to 0-based exclusive bed
  LOH$Start <- LOH$Start -1
}

# Write output file
write.table(LOH, file = paste0(prefix, "_LOH.bed"), quote = F, 
            row.names = F, col.names = F, sep = "\t")

# plot Allele ballance and LOH
SNPdata$Chr <- factor(SNPdata$Chr, levels = refChr)
if (nrow(LOH > 0)){
  ggplot() +
    geom_rect(data = LOH, mapping = aes(xmin = Start, xmax = End, ymin = 0, ymax = 1), 
              fill = "indianred1", alpha = 0.5) +
    geom_point(data = SNPdata, mapping = aes(x = Pos, y = AB), 
               size = 0.2, alpha = 0.5) +
    facet_wrap(~Chr, ncol = 2) +
    theme_minimal()
} else {
  ggplot() +
    geom_point(data = SNPdata, mapping = aes(x = Pos, y = AB), 
               size = 0.2, alpha = 0.5) +
    facet_wrap(~Chr, ncol = 2) +
    theme_minimal()
}
ggsave(filename = paste0(prefix, "_LOH.jpg"), bg="white")

# Remove temporary file
file.remove(paste0(prefix, "_HeterozygousSNPs_tmp.tsv"))

