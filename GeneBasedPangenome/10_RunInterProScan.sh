#!/bin/bash
#MSUB -r CNVcalling
#MSUB -n 1
#MSUB -c 4
#MSUB -T 86400
#MSUB -q milan
#MSUB -m scratch,work,store
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load python
module load java/11


SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/04-analysis

# InterProScan was installed in the $SCRIPTDIR following these instructions: https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html#obtaining-the-core-interproscan-software
# The folder was then converted in a squashfs image because the actual archive contains a folder with more than 50 000 inodes (forbiden by the TGCC), and the folder was deleted
# mksquashfs interproscan-5.65-97.0 interproscan-5.65-97.0.sqsh -no-xattrs

rm -rf $WORKDIR/5_InterProScan
mkdir $WORKDIR/5_InterProScan
cd $WORKDIR/5_InterProScan
mkdir mount

# Retrieve non reference genes of the pangenome
sed 's/\*//g' $DATADIR/Pangenome.NoRedundancy.pep.faa > Pangenome.NoRedundancy.pep.faa

# Run InterProScan
#$SCRIPTDIR/my_interproscan/interproscan-5.65-97.0/interproscan.sh -i Pangenome.NoRedundancy.pep.faa -f tsv -dp -goterms
ccc_mprun -C host-passthrough -E'--ctr-mount src=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/02-scripts/my_interproscan/interproscan-5.65-97.0.sqsh,dst=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/04-analysis/5_InterProScan/mount,type=squashfs' -- mount/interproscan.sh -i Pangenome.NoRedundancy.pep.faa -f tsv -dp -goterms



#ccc_mprun -m work,scratch -p milan -C host-passthrough -E'--ctr-mount src=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/02-scripts/my_interproscan/interproscan-5.65-97.0.sqsh,dst=/ccc/scratch/cont007/fg0006/loeglerv/MappingPangenome/04-analysis/5_InterProScan/mount,type=squashfs' -s
