#!/bin/bash
#MSUB -r CNVcalling
#MSUB -n 256
#MSUB -c 2
#MSUB -T 86400
#MSUB -q milan
#MSUB -m scratch,work,store
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load bwa-mem2
module load python
module load samtools
module load glost

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/04-analysis

PANGENOME=$DATADIR/Pangenome.NoRedundancy.cds.fna
READSDIR=$DATADIR/SequencingReads

if false; then
	bwa-mem2 index $PANGENOME
	cd $WORKDIR
	tail -n +2  $DATADIR/SupTab2_AsmStatistics.csv | cut -d ',' -f 1 | tr -d '"' | sort -u | grep -v AIS > Strains.txt
	#sed 's/^/,/g' Strains.txt | sed 's/$/,/g' > StrainsRegEx.txt
	#tail -n +2 $DATADIR/operationalTable3039Sace.csv | grep -v NotPublishedInLucia | grep -v -f StrainsRegEx.txt | cut -d "," -f 7 >> Strains.txt
	#rm -f StrainsRegEx.txt

	# Reads Mapping ======================================
	rm -rf $WORKDIR/1_Mapping
	mkdir $WORKDIR/1_Mapping
	cd $WORKDIR/1_Mapping

	CMD=Mapping.cmd
	rm -f $CMD

	for S in $(cat $WORKDIR/Strains.txt)
	do
		READS=$READSDIR/${S}_[1-2].fastq.gz
		if [ ! $(ls -1 $READS | wc -l) -eq 2 ]; then
			READS=$READSDIR/${S}.fastq.gz
		fi
		echo "bash $SCRIPTDIR/readsMapping.sh $S $PANGENOME $S $READS" >> $CMD
	done

	ccc_mprun glost_launch $CMD
	rm -f $CMD
fi

# Gene CN Calling ====================================
rm -rf $WORKDIR/2_CNV_Calling
mkdir $WORKDIR/2_CNV_Calling
cd $WORKDIR/2_CNV_Calling

CMD=Calling.cmd
rm -f $CMD

# Get ploidy and replace unknown ploidy by 1
#tail -n +2  $DATADIR/SupTab2_AsmStatistics.csv | cut -d ',' -f 1 | grep -v AIS | tr -d '"' | sort -u | sed 's/^/,/g' | sed 's/$/,/g' > StrainsRegEx.txt
tail -n +2  $DATADIR/SupTab2_AsmStatistics.csv | cut -d ',' -f 1,5 | grep -v AIS | tr -d '"' | sed 's/,NA$/,1/g' | tr ',' '\t' | sort -u > Ploidy.tsv
#tail -n +2 $DATADIR/operationalTable_Full3039Sace_Clades_Ploidy_Aneuploidy.csv | grep -v NotPublishedInLucia | grep -v -f StrainsRegEx.txt | cut -d "," -f 7,16 | sed 's/\?/1/g' | tr ',' '\t' >> Ploidy.tsv
#rm -f StrainsRegEx.txt

for BAM in $WORKDIR/1_Mapping/*.bam
do
	S=$(basename $BAM .bam)
	PLOIDY=$(awk -v S=$S '{if ($1 == S) print $2}' Ploidy.tsv)
	echo "python $SCRIPTDIR/callGeneCNVs.py -B $BAM -s $S -p $PLOIDY -i $DATADIR/Pangenome.NoRedundancy.cds.IdenticalRegions.tsv -x 0.15 -t 0.3 -o $S.geneCN.tsv" >> $CMD
done

ccc_mprun glost_launch $CMD
rm -f $CMD

# Matrix building =====================================
rm -rf $WORKDIR/3_CNVMatrix
mkdir $WORKDIR/3_CNVMatrix
cd $WORKDIR/3_CNVMatrix

NbS=$(ls -1 $WORKDIR/2_CNV_Calling/*.geneCN.tsv | wc -l)

echo -e "Strain\tPloidy\tGene\tDepth\tNormDepth\tCN" > full${NbS}Sace_MatrixCNV.tsv
cat $WORKDIR/2_CNV_Calling/*.geneCN.tsv | grep -v NormDepth >> full${NbS}Sace_MatrixCNV.tsv
gzip full${NbS}Sace_MatrixCNV.tsv
