#!/bin/bash
#MSUB -r minimap2
#MSUB -n 64
#MSUB -c 4
#MSUB -T 86400
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load minimap2
module load samtools
module load gatk
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
WORKDIR=$BATCHDIR/04_LRmapping
INNANO=$BATCHDIR/01_SeqData
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna

CMD=$WORKDIR/mapLR.cmd
rm -f $CMD

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	NANOPOREREADS=$INNANO/$S.NanoporeReads.fastq.gz
	BAM=$WORKDIR/$S.NanoporeReads.bam
	PREFIX=$S.NanoporeReadsMapping

	# Map LR against reference genome
	#minimap2 -ax map-ont -t 20 $REF $NANOPOREREADS | samtools sort -o $BAM -T reads.$PREFIX.tmp
	#gatk AddOrReplaceReadGroups -I $BAM -O $BAM.ReadGroups --RGID $S --RGLB $S --RGPL ILLUMINA --RGPU $S --RGSM $S
	#rm -f $BAM
	#mv $BAM.ReadGroups $BAM
	#samtools index $BAM
	echo "minimap2 -ax map-ont -t 20 $REF $NANOPOREREADS | samtools sort -o $BAM -T reads.$PREFIX.tmp && gatk AddOrReplaceReadGroups -I $BAM -O $BAM.ReadGroups --RGID $S --RGLB $S --RGPL NANOPORE --RGPU $S --RGSM $S && rm -f $BAM && mv $BAM.ReadGroups $BAM && samtools index $BAM" >> $CMD
done

ccc_mprun glost_launch $CMD
rm -f $CMD
