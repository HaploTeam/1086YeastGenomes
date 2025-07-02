library(SNPRelate)
library(ape)

# Get paths of input files
args = commandArgs(trailingOnly=TRUE)
my.vcf <- args[1]
prefix <- basename(my.vcf)

snpgdsVCF2GDS(my.vcf, paste0(prefix, ".gds"), ignore.chr.prefix="chromosome")
genofile <- snpgdsOpen(paste0(prefix, ".gds"))
snpgdsSummary(genofile)
dissMatrix <- snpgdsDiss(genofile, sample.id = NULL, snp.id = NULL, autosome.only = FALSE,remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, num.thread = 4, verbose = TRUE)
saveRDS(dissMatrix, paste0(prefix, ".rds"))
colnames(dissMatrix$diss) <- dissMatrix$sample.id
tr <- bionjs(dissMatrix$diss)
write.tree(tr, file = paste0(prefix, ".newick"), append = FALSE,digits = 10, tree.names = FALSE)


