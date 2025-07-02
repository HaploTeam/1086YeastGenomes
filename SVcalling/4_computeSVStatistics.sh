#!/bin/bash
#SBATCH -p fast               # partition all / fast
#SBATCH -n 1                      # number of cores
#SBATCH --mem-per-cpu=2000
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR

SCRIPTDIR=/home/vloegler/SVCalling_Sace1000ONT/02-scripts
DATADIR=/home/vloegler/SVCalling_Sace1000ONT/03-data
WORKDIR=/home/vloegler/SVCalling_Sace1000ONT/04-analysis

# Compute Pi and Theta for SVs
rm -rf $WORKDIR/08_SVStatistics
mkdir $WORKDIR/08_SVStatistics
cd $WORKDIR/08_SVStatistics

python $SCRIPTDIR/ComputePi_Parallel.py --vcf $WORKDIR/04_PopulationVCF/SV.1087Samples.vcf.gz --window 10000 --sliding 1000 --threads 16 --output SV.1087Samples.StructuralStatistics.10kWindow.tsv
python $SCRIPTDIR/ComputePi_Parallel.py --vcf $WORKDIR/04_PopulationVCF/SV.1087Samples.vcf.gz --window 10000 --sliding 1000 --threads 16 --output SV.1087Samples.StructuralStatistics.PAV.10kWindow.tsv --types INS DEL
python $SCRIPTDIR/ComputePi_Parallel.py --vcf $WORKDIR/04_PopulationVCF/SV.1087Samples.vcf.gz --window 10000 --sliding 1000 --threads 16 --output SV.1087Samples.StructuralStatistics.CNV.10kWindow.tsv --types CONTR DUP
python $SCRIPTDIR/ComputePi_Parallel.py --vcf $WORKDIR/04_PopulationVCF/SV.1087Samples.vcf.gz --window 10000 --sliding 1000 --threads 16 --output SV.1087Samples.StructuralStatistics.INV.10kWindow.tsv --types INV
python $SCRIPTDIR/ComputePi_Parallel.py --vcf $WORKDIR/04_PopulationVCF/SV.1087Samples.vcf.gz --window 10000 --sliding 1000 --threads 16 --output SV.1087Samples.StructuralStatistics.TRA.10kWindow.tsv --types TRA

