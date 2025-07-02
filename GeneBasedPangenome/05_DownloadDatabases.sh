#!/bin/bash
#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 10                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=400              # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o slurm.%N.%j.out      # File to which STDOUT will be written
#SBATCH -e slurm.%N.%j.err      # File to which STDERR will be written

module load ncbi-blast+/2.10.0

DATADIR=/shared/home/vloegler/Pangenome_Sace1000ONT/03-data
rm -rf $DATADIR/BlastDatabases
mkdir $DATADIR/BlastDatabases
cd $DATADIR/BlastDatabases

# REFSEQ DATABASE ======================================
# Download
for i in $(seq 0 9)
do
    srun --exclusive -N1 -n1 wget https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.0${i}.tar.gz &
done
wait

for i in $(seq 10 37)
do
    srun --exclusive -N1 -n1 wget https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.${i}.tar.gz &
done
wait

# Extract files
for i in $(seq 0 9)
do
    srun --exclusive -N1 -n1 tar -xvf refseq_protein.0${i}.tar.gz &
done
wait

for i in $(seq 10 37)
do
    srun --exclusive -N1 -n1 tar -xvf refseq_protein.${i}.tar.gz &
done
wait
rm -f *.tar.gz

# SHEN 2018 332 FUNGI ==================================
rm -rf Shen2018_Pep
# Download and get protein sequences
wget 'https://figshare.com/ndownloader/files/13092791'
mv 13092791 13092791.zip
unzip 13092791.zip
unzip 0_332yeast_genomes/332_genome_annotations.zip
mv pep Shen2018_Pep
rm -rf 0_332yeast_genomes 13092791.zip Candida_albicans_SC5314_A22_current_default_* cds gtf Saccharomyces_cerevisiae_S288C_* statement.txt

# YUE 2017 5 PARA ======================================
rm -rf Yue2017_Pep
mkdir Yue2017_Pep
cd Yue2017_Pep
# Download and extract proteins of 5 Sapa
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_PEP/CBS432.pep.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_PEP/N44.pep.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_PEP/YPS138.pep.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_PEP/UFRJ50816.pep.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_PEP/UWOPS919171.pep.fa.gz
gzip -d *.gz

# BLAST TAXDB ==========================================
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar xzvf taxdb.tar.gz
rm -f taxdb.tar.gz
