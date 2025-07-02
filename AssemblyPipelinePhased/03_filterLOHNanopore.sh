#!/bin/bash
#MSUB -r filterLOH
#MSUB -n 92
#MSUB -c 1
#MSUB -T 86400
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load r
module load bcftools
module load bedtools
module load python
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased/03_NanoporeSNPs_noLOH
VCFDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased/02_NanoporeSNPs
REFBED=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.bed

CMD=$WORKDIR/filterLOH.cmd
rm -f $CMD*

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	# Find LOH regions
	echo "Rscript LOH_Detection.r -b $REFBED -v $VCFDIR/$S.LongshotCall.vcf.gz" >> $CMD"1"
	# Get nonLOH regions
	echo "bedtools subtract -a $REFBED -b $WORKDIR/$S.LongshotCall_LOH.bed > $WORKDIR/$S.LongshotCall_NonLOH.bed" >> $CMD"2" 
	# Filter out SNPs
	echo "bcftools view --regions-file $WORKDIR/$S.LongshotCall_NonLOH.bed -Oz -o $WORKDIR/$S.LongshotCall.noLOH.vcf.gz $VCFDIR/$S.LongshotCall.vcf.gz" >> $CMD"3"
done

ccc_mprun glost_launch $CMD"1"
ccc_mprun glost_launch $CMD"2"
ccc_mprun glost_launch $CMD"3"
rm -f $CMD*
