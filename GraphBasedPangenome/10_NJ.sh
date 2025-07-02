#!/bin/bash
#MSUB -r MinigraphCactus
#MSUB -n 1
#MSUB -c 16
#MSUB -Q long
#MSUB -T 259200
#MSUB -q milan
#MSUB -m scratch,work,store
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load bcftools
module load r
module load bioconductor

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis

cd $WORKDIR/10_PopulationVCF

PREFIX=GraphGenotyping.2874Samples.Reheader.DP2.TrimAlt.Atomize.VariantsAnnotated
# Convert to plink
plink --vcf $PREFIX.vcf.gz --out $PREFIX.plink --make-bed -aec --vcf-half-call 'missing'
# Change genotypes to 1 and 2 to prevent problematic multicharacter genotypes
cut -f 1-4 $PREFIX.plink.bim | sed 's/$/\t1\t2/g' > $PREFIX.plink.bim2
mv $PREFIX.plink.bim $PREFIX.plink.bim.ExplicitAlleles
mv $PREFIX.plink.bim2 $PREFIX.plink.bim

# Make NJ tree
Rscript $SCRIPTDIR/NJTree_SV.r $PREFIX.plink

