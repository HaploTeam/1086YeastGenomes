#!/bin/bash
#MSUB -r CNVcalling
#MSUB -n 256
#MSUB -c 2
#MSUB -T 86400
#MSUB -q milan
#MSUB -m scratch,work,store
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load python
module load samtools
module load star
module load glost

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/04-analysis

PANGENOME=$DATADIR/Pangenome.NoRedundancy.cds.fna
READSDIR=$DATADIR/SequencingReads

rm -rf $WORKDIR/4_RNAmapping
mkdir $WORKDIR/4_RNAmapping
cd $WORKDIR/4_RNAmapping

mkdir PangenomeDir
zcat $WORKDIR/3_CNVMatrix/full1086Sace_MatrixCNV.tsv.gz | tail -n +2 | awk '{if ($6 > 0) print $3}' | sort -u > GenesPresent.txt
python $CCCWORKDIR/GenomeAssemblyTools/extractFragment.py -f $PANGENOME -F GenesPresent.txt -o PangenomeDir/Pangenome.$(cat GenesPresent.txt | wc -l)Genes.cds.fasta

PANGENOME=PangenomeDir/Pangenome.$(cat GenesPresent.txt | wc -l)Genes.cds.fasta
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./PangenomeDir --genomeFastaFiles $PANGENOME --genomeSAindexNbases 10

CMD=RNAmapping.cmd
rm -f $CMD*

for R in $CCCSTOREDIR/../fg0006/rawdata/1002RNAseq/run*/*_sequence.txt.gz
do
	S=$(basename $R _sequence.txt.gz | rev | cut -d '_' -f 1 | rev | sed 's/lane1//g')
	ln -s $R
	echo "STAR --runThreadN 4 --limitBAMsortRAM 2307054565 --genomeDir ./PangenomeDir --readFilesIn $(basename $R) --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $S." >> ${CMD}1
	echo "samtools index $S.Aligned.sortedByCoord.out.bam" >> ${CMD}2
	echo "echo -e 'Gene\tLength\tNbMappedReads' > $S.ReadCounts.tsv && samtools idxstats $S.Aligned.sortedByCoord.out.bam | cut -f 1-3 >> $S.ReadCounts.tsv" >> ${CMD}3
done

ccc_mprun glost_launch ${CMD}1
ccc_mprun glost_launch ${CMD}2
ccc_mprun glost_launch ${CMD}3

python $SCRIPTDIR/callTPM.py -i *.ReadCounts.tsv -o Pangenome.TPMmatrix.tsv
gzip Pangenome.TPMmatrix.tsv

rm -f *.ReadCounts.tsv *.bam *.bai *_sequence.txt.gz GenesPresent.txt ${CMD}*


