#!/bin/bash
#SBATCH -p fast               # partition all / fast
#SBATCH -n 22                      # number of cores
#SBATCH --mem-per-cpu=200
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR

# Load MiniBUSCO
source /home/vloegler/.bashrc
conda activate miniBuscoEnv

LIBRARY_PATH=BUSCO_SaccharomycetesLibraries
LINEAGE=saccharomycetes

# Collapsed
for ASSEMBLY_PATH in Final_CollapsedAssemblies_20230831/BestAssemblies/*.fasta
do
	OUTPUT_DIR=$(basename $ASSEMBLY_PATH .Final.fasta).miniBusco
	if [ ! -f $OUTPUT_DIR.txt ]; then
		srun --exclusive -N1 -n1 compleasm run -a $ASSEMBLY_PATH -o $OUTPUT_DIR -L $LIBRARY_PATH -l $LINEAGE &
	fi
done
wait

# Phased
for ASSEMBLY_PATH in Final_PhasedAssemblies_20230825/BestAssemblies/*.HP1.Final.fasta
do
	S=$(basename $ASSEMBLY_PATH .HP1.Final.fasta)
	OUTPUT_DIR=$S.Phased.miniBusco
	if [ ! -f $OUTPUT_DIR.txt ]; then
		cat Final_PhasedAssemblies_20230825/BestAssemblies/$S.HP*.Final.fasta > $S.PhasedGenome.fasta
		ASSEMBLY_PATH=$S.PhasedGenome.fasta
		srun --exclusive -N1 -n1 compleasm run -a $ASSEMBLY_PATH -o $OUTPUT_DIR -L $LIBRARY_PATH -l $LINEAGE &
	fi
done
wait

# Phenovar
for ASSEMBLY_PATH in Pheno72Asm/*.fasta
do
	S=$(basename $ASSEMBLY_PATH | cut -d "." -f 1)
	if [ $S = "AAC" ] || [ $S = "AIS" ] || [ $S = "CNT" ]; then
		OUTPUT_DIR=$S.Phased.miniBusco
		HP=$(basename $ASSEMBLY_PATH | cut -d "." -f 2)
		if [ $HP == "HP1" ]; then
			if [ ! -f $OUTPUT_DIR.txt ]; then
				cat Pheno72Asm/$S.HP*.fasta > $S.PhasedGenome.fasta
				ASSEMBLY_PATH=$S.PhasedGenome.fasta
				srun --exclusive -N1 -n1 compleasm run -a $ASSEMBLY_PATH -o $OUTPUT_DIR -L $LIBRARY_PATH -l $LINEAGE &
			fi
		fi
	else
		OUTPUT_DIR=$S.miniBusco
		if [ ! -f $OUTPUT_DIR.txt ]; then
			srun --exclusive -N1 -n1 compleasm run -a $ASSEMBLY_PATH -o $OUTPUT_DIR -L $LIBRARY_PATH -l $LINEAGE &
		fi
	fi
done
wait


for S in */summary.txt
do
	mv $S $(echo $S | cut -d "/" -f 1).txt
done

rm -rf *.miniBusco *.fasta

# Get busco results
for F in *.miniBusco.txt
do
	S=$(grep "^S" $F | cut -d ":" -f 2 | cut -d "%" -f 1) # Single presence of gene
	D=$(grep "^D" $F | cut -d ":" -f 2 | cut -d "%" -f 1) # Duplicated presence of gene
	SCORE=$(echo "${S}+${D}" | bc)
	echo -e $(basename $F .miniBusco.txt)"\t"$SCORE
done > BuscoScores.tsv




