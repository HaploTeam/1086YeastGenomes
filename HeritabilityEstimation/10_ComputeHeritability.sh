#!/bin/bash
#MSUB -r FilterSNPs
#MSUB -n 128
#MSUB -c 1
#MSUB -T 86400
#MSUB -m scratch,work,store
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load r
module load glost

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/GWAS_Sace1000ONT/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/GWAS_Sace1000ONT/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GWAS_Sace1000ONT/04-analysis

LDAK=$SCRIPTDIR/ldak5.2.linux

if false; then
	# Create grm
	rm -rf $WORKDIR/10_KinshipMatrices
	mkdir $WORKDIR/10_KinshipMatrices
	cd $WORKDIR/10_KinshipMatrices

	# SNPs ================
	# Change genotypes to 1 and 2 to prevent problematic multicharacter genotypes
	cp $WORKDIR/01_Genotypes/full966.SNPs.plink.* .
	cut -f 1-4 full966.SNPs.plink.bim | sed 's/$/\t1\t2/g' > full966.SNPs.plink.bim2
	rm -f full966.SNPs.plink.bim
	mv full966.SNPs.plink.bim2 full966.SNPs.plink.bim
	# Compute variant weights with LDAK-Thin model
	$LDAK --thin thin.SNPs --bfile full966.SNPs.plink --window-prune .98 --window-kb 20
	awk < thin.SNPs.in '{print $1, 1}' > weights.SNPs.thin
	# Compute kinship
	$LDAK --calc-kins-direct LDAK-Thin.SNPs --bfile full966.SNPs.plink --weights weights.SNPs.thin --power -.25

	# INDELs ==============
	# Change genotypes to 1 and 2 to prevent problematic multicharacter genotypes
	cp $WORKDIR/01_Genotypes/full966.INDELs.plink.* .
	cut -f 1-4 full966.INDELs.plink.bim | sed 's/$/\t1\t2/g' > full966.INDELs.plink.bim2
	rm -f full966.INDELs.plink.bim
	mv full966.INDELs.plink.bim2 full966.INDELs.plink.bim
	# Compute variant weights with LDAK-Thin model
	$LDAK --thin thin.INDELs --bfile full966.INDELs.plink --window-prune .98 --window-kb 20
	awk < thin.INDELs.in '{print $1, 1}' > weights.INDELs.thin
	# Compute kinship
	$LDAK --calc-kins-direct LDAK-Thin.INDELs --bfile full966.INDELs.plink --weights weights.INDELs.thin --power -.25

	# SVs ==============
	# Change genotypes to 1 and 2 to prevent problematic multicharacter genotypes
	cp $WORKDIR/01_Genotypes/full966.SVs.plink.* .
	cut -f 1-4 full966.SVs.plink.bim | sed 's/$/\t1\t2/g' > full966.SVs.plink.bim2
	rm -f full966.SVs.plink.bim
	mv full966.SVs.plink.bim2 full966.SVs.plink.bim
	## Compute variant weights with LDAK-Thin model
	#$LDAK --thin thin.SVs --bfile full966.SVs.plink --window-prune .98 --window-kb 20
	#awk < thin.SVs.in '{print $1, 1}' > weights.SVs.thin
	## Compute kinship
	#$LDAK --calc-kins-direct LDAK-Thin.SVs --bfile full966.SVs.plink --weights weights.SVs.thin --power -.25

	# CNVs ==============
	# Change genotypes to 1 and 2 to prevent problematic multicharacter genotypes
	cp $WORKDIR/01_Genotypes/full966.CNVs.plink.* .
	cut -f 1-4 full966.CNVs.plink.bim | sed 's/$/\t1\t2/g' > full966.CNVs.plink.bim2
	rm -f full966.CNVs.plink.bim
	mv full966.CNVs.plink.bim2 full966.CNVs.plink.bim
	## Compute variant weights with LDAK-Thin model
	#$LDAK --thin thin.CNVs --bfile full966.CNVs.plink --window-prune .98 --window-kb 20
	#awk < thin.CNVs.in '{print $1, 1}' > weights.CNVs.thin
	## Compute kinship
	#$LDAK --calc-kins-direct LDAK-Thin.CNVs --bfile full966.CNVs.plink --weights weights.CNVs.thin --power -.25

	# CNVs + INDELs ==============
	# Merge genotypes
	plink --bfile full966.SVs.plink --bmerge full966.CNVs.plink --make-bed --out full966.SVs.CNVs.plink
	# Compute variant weights with LDAK-Thin model
	$LDAK --thin thin.SVs --bfile full966.SVs.CNVs.plink --window-prune .98 --window-kb 20
	awk < thin.SVs.in '{print $1, 1}' > weights.SVs.thin
	# Compute kinship
	$LDAK --calc-kins-direct LDAK-Thin.SVs --bfile full966.SVs.CNVs.plink --weights weights.SVs.thin --power -.25

	# Create file with path to multiple kinship matrices
	(echo "${PWD}/LDAK-Thin.SNPs") > LDAK-Thin.SNPs.txt
	(echo "${PWD}/LDAK-Thin.SNPs"
	echo "${PWD}/LDAK-Thin.INDELs") > LDAK-Thin.SNPs.INDELs.txt
	(echo "${PWD}/LDAK-Thin.SNPs"
	echo "${PWD}/LDAK-Thin.INDELs"
	echo "${PWD}/LDAK-Thin.SVs") > LDAK-Thin.SNPs.INDELs.SVs.txt
	(echo "${PWD}/LDAK-Thin.SNPs"
	echo "${PWD}/LDAK-Thin.INDELs"
	echo "${PWD}/LDAK-Thin.SVs") > LDAK-Thin.SNPs.INDELs.SVs.txt

	# Combine all matrices together and run PCA with plink on it
	plink --bfile full966.SNPs.plink --bmerge full966.INDELs.plink --make-bed --out full966.SNPs.INDELs.plink
	plink --bfile full966.SNPs.INDELs.plink --bmerge full966.SVs.CNVs.plink --make-bed --out full966.SNPs.INDELs.SVs.CNVs.plink
	rm -f full966.SNPs.INDELs.plink.*

	# PCA
	plink --bfile full966.SNPs.INDELs.SVs.CNVs.plink --pca --out full966.SNPs.INDELs.SVs.CNVs.plink

	# Keep the four first PCAs
	cut -d ' ' -f 1-6 full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec > full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.4PC
	# Merge with ploidy
	cut -d ' ' -f 1-2 full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.4PC | tr ' ' '\t' > s.plink
	grep -f s.plink $DATADIR/1086Samples_Ploidy.covar.tsv > 966Samples_Ploidy.covar.tsv
	# Sort files to ensure their order
	sort full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.4PC | tr ' ' '\t' > full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.4PC.sorted
	sort 966Samples_Ploidy.covar.tsv | cut -f 3 > 966Samples_Ploidy.covar.sorted.NoSampleName.tsv
	# Combine covariant files
	paste full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.4PC.sorted 966Samples_Ploidy.covar.sorted.NoSampleName.tsv > full966.4PCCovar.PloidyCovar.plink.tsv

	# Keep the 10 first PCA
	cut -d ' ' -f 1-12 full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec > full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.10PC
	# Merge with ploidy
	cut -d ' ' -f 1-2 full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.10PC | tr ' ' '\t' > s.plink
	grep -f s.plink $DATADIR/1086Samples_Ploidy.covar.tsv > 966Samples_Ploidy.covar.tsv
	# Sort files to ensure their order
	sort full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.10PC | tr ' ' '\t' > full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.10PC.sorted
	sort 966Samples_Ploidy.covar.tsv | cut -f 3 > 966Samples_Ploidy.covar.sorted.NoSampleName.tsv
	# Combine covariant files
	paste full966.SNPs.INDELs.SVs.CNVs.plink.eigenvec.10PC.sorted 966Samples_Ploidy.covar.sorted.NoSampleName.tsv > full966.10PCCovar.PloidyCovar.plink.tsv

	# PCA SNPs
	plink --bfile full966.SNPs.plink --pca --out full966.SNPs.plink

	# Keep the 10 first PCA
	cut -d ' ' -f 1-12 full966.SNPs.plink.eigenvec > full966.SNPs.plink.eigenvec.10PC
	# Merge with ploidy
	cut -d ' ' -f 1-2 full966.SNPs.plink.eigenvec.10PC | tr ' ' '\t' > s.plink
	grep -f s.plink $DATADIR/1086Samples_Ploidy.covar.tsv > 966Samples_Ploidy.covar.tsv
	# Sort files to ensure their order
	sort full966.SNPs.plink.eigenvec.10PC | tr ' ' '\t' > full966.SNPs.plink.eigenvec.10PC.sorted
	sort 966Samples_Ploidy.covar.tsv | cut -f 3 > 966Samples_Ploidy.covar.sorted.NoSampleName.tsv
	# Combine covariant files
	paste full966.SNPs.plink.eigenvec.10PC.sorted 966Samples_Ploidy.covar.sorted.NoSampleName.tsv > full966.10PCCovarSNPs.PloidyCovar.plink.tsv

	# Compute heritability ===============
	rm -rf $WORKDIR/11_Heritability
	mkdir $WORKDIR/11_Heritability
	cd $WORKDIR/11_Heritability

	mkdir $WORKDIR/11_Heritability/NormPheno
	cp $DATADIR/PhenotypesPlink/*.phen $WORKDIR/11_Heritability/NormPheno/
	CMD=NormalizePheno.cmd
	rm -f $CMD

	for PHENO in $WORKDIR/11_Heritability/NormPheno/*.plink.phen
	do
		echo "Rscript $SCRIPTDIR/rank_based_Inverse_Normal_Transformation.PlinkFormat.R $PHENO" >> $CMD
	done
	ccc_mprun glost_launch $CMD
	rm -f $CMD
	rm -f $WORKDIR/11_Heritability/NormPheno/*.plink.phen
fi

cd $WORKDIR/11_Heritability

CMD=EstimateHeritability.cmd
rm -f $CMD

# Pheno with a variant associated with PVal < 1e-20
cat $WORKDIR/02_GWAS_AllVariants/Results/*/*.signif_snps.txt | awk '{if ($5 <= 1E-20) print $0}' | cut -f 13 | sort -u > NonComplexPhenotypes.txt
mkdir K_SNPs K_SNPs.INDELs K_SNPs.INDELs.SVs
mkdir K_SNPs_4PC K_SNPs.INDELs_4PC K_SNPs.INDELs.SVs_4PC
mkdir K_SNPs_10PC K_SNPs.INDELs_10PC K_SNPs.INDELs.SVs_10PC

for PHENO in $WORKDIR/11_Heritability/NormPheno/*.phen
do
	COND=$(basename $PHENO .plink.norm.phen)
	# Look if pheno is complex or not
	if ! grep -Fxq $COND NonComplexPhenotypes.txt; then
		for V in SNPs SNPs.INDELs SNPs.INDELs.SVs
		do
			#echo "$LDAK --reml K_${V}/$COND --pheno $PHENO --covar $DATADIR/1086Samples_Ploidy.covar.tsv --mgrm $WORKDIR/10_KinshipMatrices/LDAK-Thin.$V.txt --constrain YES && rm -f K_${V}/$COND.coeff K_${V}/$COND.cross K_${V}/$COND.indi.* K_${V}/$COND.progress K_${V}/$COND.reg.* K_${V}/$COND.share K_${V}/$COND.vars" >> $CMD
			#echo "$LDAK --reml K_${V}_4PC/$COND --pheno $PHENO --covar $WORKDIR/10_KinshipMatrices/full966.4PCCovar.PloidyCovar.plink.tsv --mgrm $WORKDIR/10_KinshipMatrices/LDAK-Thin.$V.txt --constrain YES && rm -f K_${V}_4PC/$COND.coeff K_${V}_4PC/$COND.cross K_${V}_4PC/$COND.indi.* K_${V}_4PC/$COND.progress K_${V}_4PC/$COND.reg.* K_${V}_4PC/$COND.share K_${V}_4PC/$COND.vars" >> $CMD
			#echo "$LDAK --reml K_${V}_10PC/$COND --pheno $PHENO --covar $WORKDIR/10_KinshipMatrices/full966.10PCCovar.PloidyCovar.plink.tsv --mgrm $WORKDIR/10_KinshipMatrices/LDAK-Thin.$V.txt --constrain YES && rm -f K_${V}_10PC/$COND.coeff K_${V}_10PC/$COND.cross K_${V}_10PC/$COND.indi.* K_${V}_10PC/$COND.progress K_${V}_10PC/$COND.reg.* K_${V}_10PC/$COND.share K_${V}_10PC/$COND.vars" >> $CMD
			echo "$LDAK --reml K_${V}_10PCSNPs/$COND --pheno $PHENO --covar $WORKDIR/10_KinshipMatrices/full966.10PCCovarSNPs.PloidyCovar.plink.tsv --mgrm $WORKDIR/10_KinshipMatrices/LDAK-Thin.$V.txt --constrain YES && rm -f K_${V}_10PC/$COND.coeff K_${V}_10PC/$COND.cross K_${V}_10PC/$COND.indi.* K_${V}_10PC/$COND.progress K_${V}_10PC/$COND.reg.* K_${V}_10PC/$COND.share K_${V}_10PC/$COND.vars" >> $CMD
		done
	fi
done

ccc_mprun glost_launch $CMD
rm -f $CMD

if false; then
	# Gather results: No PC as covariant
	(echo -e 'Pheno\tKinships\tHer_All\tHer_SNPs\tHer_INDELs\tHer_SVs'
	for F in K_SNPs/*.reml
	do
		Pheno=$(basename $F .reml)
		for Kinship in SNPs SNPs.INDELs SNPs.INDELs.SVs
		do
			cat K_${Kinship}/$Pheno.reml | grep ^Her | grep -v Her_Top | sort | cut -d ' ' -f 2 | tr '\n' '\t' | sed 's/\t$/\n/g' | sed "s/^/$Pheno\t$Kinship\t/g"
		done
	done) | gzip > HeritabilityResults.tsv.gz
	# Gather results: 4 PC as covariant
	(echo -e 'Pheno\tKinships\tHer_All\tHer_SNPs\tHer_INDELs\tHer_SVs'
	for F in K_SNPs_4PC/*.reml
	do
		Pheno=$(basename $F .reml)
		for Kinship in SNPs SNPs.INDELs SNPs.INDELs.SVs
		do
			cat K_${Kinship}_4PC/$Pheno.reml | grep ^Her | grep -v Her_Top | sort | cut -d ' ' -f 2 | tr '\n' '\t' | sed 's/\t$/\n/g' | sed "s/^/$Pheno\t$Kinship\t/g"
		done
	done) | gzip > HeritabilityResults.4PCCovar.tsv.gz
	# Gather results: 10 PC as covariant
	(echo -e 'Pheno\tKinships\tHer_All\tHer_SNPs\tHer_INDELs\tHer_SVs'
	for F in K_SNPs_10PC/*.reml
	do
		Pheno=$(basename $F .reml)
		for Kinship in SNPs SNPs.INDELs SNPs.INDELs.SVs
		do
			cat K_${Kinship}_10PC/$Pheno.reml | grep ^Her | grep -v Her_Top | sort | cut -d ' ' -f 2 | tr '\n' '\t' | sed 's/\t$/\n/g' | sed "s/^/$Pheno\t$Kinship\t/g"
		done
	done) | gzip > HeritabilityResults.10PCCovar.tsv.gz
fi





