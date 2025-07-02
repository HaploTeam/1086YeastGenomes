#!/bin/bash
#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem=400              # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o slurm.%N.%j.out      # File to which STDOUT will be written
#SBATCH -e slurm.%N.%j.err      # File to which STDERR will be written

module load ncbi-blast+


SCRIPTDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts
DATADIR=/shared/home/vloegler/Pangenome_Sace1000ONT/03-data
WORKDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis

rm -rf $WORKDIR/08_BlastOnRefSeq
mkdir $WORKDIR/08_BlastOnRefSeq
mkdir $WORKDIR/08_BlastOnRefSeq/BlastDatabases
cd $WORKDIR/08_BlastOnRefSeq/BlastDatabases

# Create Blast database with 332 Fungi species
# ============================================

# First change file name and sequences names
# Some species are not found in the NCBI taxonomy DB because of typos
# Let's change names for these species
# Write equivalence species name in a file
echo "ambrosiozyma kashinagacola	 ambrosiozyma kashinagicola	
ashbya aceri	Ashbya sp. FD-2008
blastobotrys nivea	Blastobotrys niveus
candida azyma	Wickerhamiella azyma
candida carpophila	Meyerozyma carpophila
hanseniaspora vinae	hanseniaspora vineae	
magnusiomyces tetrasperma	Magnusiomyces tetraspermus
martiniozyma abiesophila	martiniozyma abietophila
metschnikowia dekortum	metschnikowia dekortorum
metschnikowia lockheadii	Metschnikowia lochheadii
metschnikowia matae maris	Metschnikowia matae var maris
ogataea methylivora	Ogataea methylovora
ogataea populiabae	Ogataea populi-albae" > SpeciesEqName.tsv

rm -rf Shen2018_Pep_Renamed
mkdir Shen2018_Pep_Renamed
for F in $DATADIR/BlastDatabases/Shen2018_Pep/*.max.pep
do
	SPECIES=$(basename $F .max.pep | tr '_' '\n'| grep -v [0-9] | tr '\n' ' ' | rev | cut -c 2- | rev)
	# Check if species has to be changed
	NAMETOCHANGE=$(awk -v S="$SPECIES" -F "\t" '{if ($1 == S) print $2}' SpeciesEqName.tsv | wc -l)
	if [ $NAMETOCHANGE == 1 ]
	then
		SPECIES=$(awk -v S="$SPECIES" -F "\t" '{if ($1 == S) print $2}' SpeciesEqName.tsv)
	fi
	# Change protein names
	PREFIX=$(echo $SPECIES | tr ' ' '_' | tr '[:upper:]' '[:lower:]')
	python $SCRIPTDIR/PythonScripts/ChangeSeqID.py -f $F -p $PREFIX -o Shen2018_Pep_Renamed/$PREFIX.pep.fa
done

# Get taxid for each species
rm -f taxids.txt
for F in Shen2018_Pep_Renamed/*.pep.fa
do
	SPECIES=$(basename $F .pep.fa | tr '_' ' ')
	echo $SPECIES
	# Get Taxid for species
	TAXID=$(esearch -db taxonomy -query "$SPECIES" | efetch -format taxid)
	if [ -z $TAXID ]; then # If variable empty, rerun taxid search
		TAXID=$(esearch -db taxonomy -query "$SPECIES" | efetch -format taxid)
	fi
	if [ -z $TAXID ]; then # If variable empty, rerun taxid search
		TAXID=$(esearch -db taxonomy -query "$SPECIES" | efetch -format taxid)
	fi
	echo -e "$(basename $F)\t$TAXID" >> taxids.txt
done

# Create Taxid map file
rm -f TaxidMap_Shen2018
for F in Shen2018_Pep_Renamed/*.pep.fa
do
	# Get Taxid for species
	TAXID=$(grep $(basename $F) taxids.txt | cut -f 2)
	# Get gene names
	cat $F | grep '>' | cut -f 1 -d " " | tr -d ">" > geneID.txt
	# Add taxid at the end of each line
	sed -i "s/$/\t${TAXID}/g" geneID.txt
	cat geneID.txt >> TaxidMap_Shen2018
done
rm -f geneID.txt

cat Shen2018_Pep_Renamed/*.pep.fa > Shen2018.pep.fa
makeblastdb -in Shen2018.pep.fa -dbtype prot -title Shen2018_DB -out Shen2018_DB -taxid_map TaxidMap_Shen2018 -parse_seqids
rm -f Shen2018.pep.fa

# Create Blast database with 5 Sapa strains
# =========================================
rm -rf Yue2017_Pep_Renamed
mkdir Yue2017_Pep_Renamed
for F in $DATADIR/BlastDatabases/Yue2017_Pep/*.pep.fa
do
	STRAIN=$(basename $F .pep.fa)
	# Change protein names
	python $SCRIPTDIR/PythonScripts/ChangeSeqID.py -f $F -p $STRAIN -o Yue2017_Pep_Renamed/$STRAIN.pep.fa
done

cat Yue2017_Pep_Renamed/*.pep.fa > Yue2017.pep.fa
makeblastdb -in Yue2017.pep.fa -dbtype prot -title Yue2017_DB -out Yue2017_DB -taxid 27291 -parse_seqids
rm -f Yue2017.pep.fa


# Merge the 2 databases with RefSeq
# =================================
blastdb_aliastool -dblist "$DATADIR/BlastDatabases/refseq_protein Yue2017_DB Shen2018_DB" -dbtype prot -out RefSeq_Shen2018_Yue2017 -title "RefSeq_Shen2018_Yue2017"
