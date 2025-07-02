#!/bin/bash
#MSUB -r whatshap
#MSUB -n 92
#MSUB -c 1
#MSUB -T 86400
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load whatshap
module load bcftools
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
WORKDIR=$BATCHDIR/05_phasedVCF
INBAM=$BATCHDIR/04_LRmapping
INVCF=$BATCHDIR/03_NanoporeSNPs_noLOH
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna

CMD=$WORKDIR/phaseVCF.cmd
rm -f $CMD

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	BAM=$INBAM/$S.NanoporeReads.bam
	VCF=$INVCF/$S.LongshotCall.noLOH.vcf.gz
	outVCF=$WORKDIR/$S.LongshotCall.noLOH.phased.vcf.gz

	#bcftools index $VCF
	#whatshap phase -o $outVCF --reference $REF $VCF $BAM
	#bcftools index $outVCF
	echo "bcftools index -f $VCF && whatshap phase -o $outVCF --reference $REF $VCF $BAM && bcftools index $outVCF" >> $CMD
done

ccc_mprun glost_launch $CMD
rm -f $CMD
