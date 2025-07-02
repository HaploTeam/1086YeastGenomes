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

module load blast+
module load python
module load minimap2
module load samtools

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis

rm -rf $WORKDIR/12_GraphUniqueSequences
mkdir $WORKDIR/12_GraphUniqueSequences
cd $WORKDIR/12_GraphUniqueSequences

# Get old Minigraph segments
cp $WORKDIR/04_MC_GraphPangenome/SacePangenomeGraph.500Haplotypes.gfa.fa.gz .
gzip -d SacePangenomeGraph.500Haplotypes.gfa.fa.gz

# Filter segments larger than 100 bp
python $CCCWORKDIR/GenomeAssemblyTools/filterContigSize.py -f SacePangenomeGraph.500Haplotypes.gfa.fa -m 0.1 -o SacePangenomeGraph.500Haplotypes.gfa

# Run blast of segments on themselves
blastn -query SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.fasta -subject SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.fasta -perc_identity 95 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -out SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.Self.blastn

# Identify reference segments (100% identity and 100% coverage with the reference genome)
blastn -query SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.fasta -subject $DATADIR/Sace_S288c_reference_FullMatrixID.fna -outfmt "6 qseqid sseqid qcovhsp" -perc_identity 100 -max_target_seqs 1 -max_hsps 1 | awk '{if ($3 == 100) print $1}' > RefSegments.txt

# Remove redundancy and keep unique sequences
python $SCRIPTDIR/getUniqueSequences.py -f SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.fasta -r RefSegments.txt -b SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.Self.blastn

# Extract InRef segments
grep _InRef SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.fasta | tr -d '>' > NonRedundantInRefSegments.txt
python $CCCWORKDIR/GenomeAssemblyTools/extractFragment.py -f SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.fasta -F NonRedundantInRefSegments.txt > SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.InRef.fasta

# Extract OutRef segments
grep _OutRef SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.fasta | tr -d '>' > NonRedundantOutRefSegments.txt
python $CCCWORKDIR/GenomeAssemblyTools/extractFragment.py -f SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.fasta -F NonRedundantOutRefSegments.txt > SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.OutRef.fasta

# Search introgressed segments
# Build the Saccharomyces genomes database (Without Sace)
cat $DATADIR/SaccharomycesGenomes/*.fna > tmp.fna
makeblastdb -dbtype nucl -title SaccharomycesGenomes -out SaccharomycesGenomes -taxid_map $DATADIR/SaccharomycesGenomes/SaccharomycesGenomes.taxid -parse_seqids -in tmp.fna
rm -f tmp.fna

# Align non redundant segments not aligning on Sace ref on the Saccharomyces genomes database
blastn -query SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.OutRef.fasta -db SaccharomycesGenomes -perc_identity 95 -outfmt "6 qseqid sseqid staxids pident qcovs" | awk '{if ($5 > 50) print $1}' | sort -u > NonRedundantOutRefSegments.Introgressed.txt

# Align all outRef segments on Sace ref on the Saccharomyces genomes database
grep _OutRef SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.Annotated.fasta | tr -d '>' > SegmentsOutRef.txt
python $CCCWORKDIR/GenomeAssemblyTools/extractFragment.py -f SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.Annotated.fasta -F SegmentsOutRef.txt > SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.Annotated.OutRef.fasta
# Get introgressed nodes
blastn -query SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.Annotated.OutRef.fasta -db SaccharomycesGenomes -perc_identity 95 -outfmt "6 qseqid sseqid staxids pident qcovs" | awk '{if ($5 > 50) print $1}' | sort -u > SegmentsOutRef.Introgressed.txt

# Get length of reference and non reference segments, and introgressed segments
grep -v '>' SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.InRef.fasta | tr -d '\n' | wc -c
grep -v '>' SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.OutRef.fasta | tr -d '\n' | wc -c
grep -v '>' SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.NoRedundancy.OutRef.Introgressed.fasta | tr -d '\n' | wc -c


# Map redundant segments on reference
minimap2 -ax asm5 $DATADIR/Sace_S288c_reference_FullMatrixID.fna SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.fasta | samtools sort -o SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.MappingOnRef.bam
samtools depth -aa SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.MappingOnRef.bam | gzip > SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.MappingOnRef.bam.depth.tsv.gz
