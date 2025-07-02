#!/bin/bash
#MSUB -r MinigraphCactus
#MSUB -n 1
#MSUB -c 4
#MSUB -Q long
#MSUB -T 259200
#MSUB -q milan
#MSUB -m scratch,work,store
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load bcftools
module load r
module load bioconductor
module load glost

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis

cd $WORKDIR/10_PopulationVCF

if false; then
	CMD=Reheader.cmd

	# Reheader: Change header of vcf files because of wrong chromosome ordering

	# Extract header
	bcftools reheader --fai $DATADIR/Sace_S288c_reference_FullMatrixID.fna.fai GraphGenotyping.2874Samples.chromosome1.vcf.gz | bcftools view | head -4000 | grep "^#" > header.txt
	# Extract lines before contigs
	sed -n '/##contig=<ID=chromosome/q;p' header.txt > header1.txt
	# Extract lines after contigs
	sed -n '/##contig=<ID=chromosome/,$p' header.txt | grep -v '##contig=<ID=chromosome' > header3.txt
	# Create ordered contig lines
	cat $DATADIR/Sace_S288c_reference_FullMatrixID.fna.fai | cut -f 1,2 | sed 's/^/##contig=<ID=/g' | sed 's/\t/,length=/g' | sed 's/$/>/g' > header2.txt
	# Reassemble header
	cat header1.txt header2.txt header3.txt > header.txt
	for C in $(seq 1 16)
	do
		echo "bcftools reheader --header header.txt GraphGenotyping.2874Samples.chromosome$C.vcf.gz > GraphGenotyping.2874Samples.chromosome$C.Reheader.vcf.gz"
	done > $CMD

	ccc_mprun glost_launch $CMD
	rm -f $CMD

	# Filter DP2
	CMD=Filtering.cmd
	# Quality filtering
	for C in $(seq 1 16)
	do
		echo "bcftools view -M500 GraphGenotyping.2874Samples.chromosome$C.Reheader.vcf.gz | bcftools +setGT -- -t q -n . -e 'FMT/DP>=2' | bcftools +fill-tags | bcftools view --min-ac 1 -Oz -o GraphGenotyping.2874Samples.chromosome$C.Reheader.DP2.vcf.gz"
	done > $CMD
	ccc_mprun glost_launch $CMD
	rm -f $CMD

	# Trim alt allele
	CMD=Filtering.cmd
	# Quality filtering
	for C in $(seq 1 16)
	do
		echo "bcftools view --trim-alt-alleles GraphGenotyping.2874Samples.chromosome$C.Reheader.DP2.vcf.gz | bcftools +fill-tags | bcftools view --min-ac 1 -Oz -o GraphGenotyping.2874Samples.chromosome$C.Reheader.DP2.TrimAlt.vcf.gz"
	done > $CMD
	ccc_mprun glost_launch $CMD
	rm -f $CMD

	# Atomize
	CMD=Filtering.cmd
	# Quality filtering (Only consider variant with less than 500 alleles to avoid OOM, Drop the FORMAT/GL field because of incompatibilities with trim-alt-alleles)
	for C in $(seq 1 16)
	do
		echo "bcftools view -M500 GraphGenotyping.2874Samples.chromosome$C.Reheader.DP2.TrimAlt.vcf.gz | bcftools norm --atomize --atom-overlap "." --multiallelics +any --threads 16 | bcftools +fill-tags | bcftools annotate -x FORMAT/GL | bcftools view --trim-alt-alleles | bcftools +fill-tags | bcftools view --min-ac 1 -Oz -o GraphGenotyping.2874Samples.chromosome$C.Reheader.DP2.TrimAlt.Atomize.vcf.gz"
	done > $CMD
	ccc_mprun glost_launch $CMD
	rm -f $CMD

	# Merge chromosomes
	bcftools concat GraphGenotyping.2874Samples.chromosome1.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome2.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome3.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome4.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome5.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome6.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome7.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome8.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome9.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome10.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome11.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome12.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome13.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome14.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome15.Reheader.DP2.TrimAlt.Atomize.vcf.gz \
		GraphGenotyping.2874Samples.chromosome16.Reheader.DP2.TrimAlt.Atomize.vcf.gz | \
		bcftools +fill-tags | \
		bcftools annotate --set-id 'rs_%CHROM\_%POS\_%ID' | \
		bcftools view --min-ac 1 -i 'F_MISSING<1' -Oz -o GraphGenotyping.2874Samples.Reheader.DP2.TrimAlt.Atomize.VariantsAnnotated.vcf.gz --threads 16
fi

# Merge chromosomes raw matrix
bcftools concat GraphGenotyping.2874Samples.chromosome1.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome2.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome3.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome4.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome5.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome6.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome7.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome8.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome9.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome10.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome11.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome12.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome13.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome14.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome15.Reheader.vcf.gz \
	GraphGenotyping.2874Samples.chromosome16.Reheader.vcf.gz | \
	bcftools +fill-tags | \
	bcftools annotate --set-id 'rs_%CHROM\_%POS\_%ID' | \
	bcftools view --min-ac 1 -i 'F_MISSING<1' -Oz -o GraphGenotyping.2874Samples.Reheader.VariantsAnnotated.vcf.gz --threads 16
