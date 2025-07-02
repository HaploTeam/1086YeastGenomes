#!/bin/bash
#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 16
#SBATCH -c 4                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=400              # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o slurm.%N.%j.out      # File to which STDOUT will be written
#SBATCH -e slurm.%N.%j.err      # File to which STDERR will be written

module load ncbi-blast+

SCRIPTDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts
DATADIR=/shared/home/vloegler/Pangenome_Sace1000ONT/03-data
WORKDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis

rm -rf $WORKDIR/08_BlastOnRefSeq/BlastResults
mkdir $WORKDIR/08_BlastOnRefSeq/BlastResults
cd $WORKDIR/08_BlastOnRefSeq/BlastResults

# Subset Pangenome.pep.faa to parallelize the blast search
# Only non reference genes (genes not beginning by Y or Q) are kept
PANGENOME=$WORKDIR/07_GeneClustering_HQStrains/Pangenome.NoRedundancy.pep.faa
grep '>' $PANGENOME | grep ">YX" | cut -d ' ' -f 1 | tr -d '>' > PangenomeNonRef.ID.txt
i=1
mkdir PangenomeNonRef_Subsets
for ID in $(cat PangenomeNonRef.ID.txt)
do
	echo $ID >> PangenomeNonRef_Subsets/PangenomeNonRef_Subset$i.txt
	let i++
	if [ $i == 17 ]; then
		i=1
	fi
done
rm -f PangenomeNonRef.ID.txt

for F in PangenomeNonRef_Subsets/PangenomeNonRef_Subset*.txt
do
	python $SCRIPTDIR/GenomeAssemblyTools/extractFragment.py -f $PANGENOME -F $F -o PangenomeNonRef_Subsets/$(basename $F .txt).fa
done
rm -f PangenomeNonRef_Subsets/PangenomeNonRef_Subset*.txt

# Download Blast taxDB
export BLASTDB=$DATADIR/BlastDatabases
#export BLASTDB=$WORKDIR/08_BlastOnRefSeq/BlastDatabases
for F in PangenomeNonRef_Subsets/*.fa
do
	srun --exclusive -N1 -n1 blastp -query $F -db $WORKDIR/08_BlastOnRefSeq/BlastDatabases/RefSeq_Shen2018_Yue2017 -out $(basename $F .fa).blastp -num_threads 20 -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore length qlen slen qcovs sscinames scomnames staxids" &
done
wait

cat PangenomeNonRef_Subset*.blastp > PangenomeNonRef_on_RefSeq_Shen2018_Yue2017.blastp
rm -f PangenomeNonRef_Subset*.blastp
