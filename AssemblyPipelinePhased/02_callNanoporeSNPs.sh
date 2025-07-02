#!/bin/bash
#MSUB -r Longshot
#MSUB -n 64
#MSUB -c 4
#MSUB -T 86400
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load minimap2
module load gatk
module load samtools
module load longshot
module load bcftools
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased/02_NanoporeSNPs
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna

CMD=$WORKDIR/NanoporeCalling.cmd
rm -f $CMD

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	NANOPOREREADS=$BATCHDIR/01_SeqData/$S.NanoporeReads.fastq.gz
	BAM=$WORKDIR/$S.NanoporeReads.bam
	PREFIX=$S.NanoporeReadsMapping
	# Minimap2 read mapping
	echo "minimap2 -ax map-ont -t 20 $REF $NANOPOREREADS | samtools sort -o $BAM -T reads.$PREFIX.tmp && gatk AddOrReplaceReadGroups -I $BAM -O $BAM.ReadGroups --RGID $S --RGLB $S --RGPL NANOPORE --RGPU $S --RGSM $S && rm -f $BAM && mv $BAM.ReadGroups $BAM && samtools index $BAM" >> ${CMD}1
	# Longshot SNP calling
	OUT=$WORKDIR/$S.LongshotCall.vcf
    echo "longshot --bam $BAM --ref $REF --no_haps --min_cov 7 --min_alt_count 7 --min_alt_frac 0.2 --sample_id $S --out $OUT" >> ${CMD}2
    echo "bcftools view -O z -o $OUT.gz $OUT && bcftools index $OUT.gz" >> ${CMD}3
done

ccc_mprun glost_launch ${CMD}1
ccc_mprun glost_launch ${CMD}2
ccc_mprun glost_launch ${CMD}3

rm -f $CMD*
rm -f $WORKDIR/*.bam
rm -f $WORKDIR/*.bam.bai
rm -f $WORKDIR/*.vcf
