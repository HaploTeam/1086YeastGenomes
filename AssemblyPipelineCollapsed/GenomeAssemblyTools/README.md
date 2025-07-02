# GenomeAssemblyTools

This repository contains usefull python tools to deal with de novo genome assemblies using MUMmer4 and blast+. 

**Tools.py**
* Contains the main functions and classes used by the other script. 

## Processing fasta files
**extractFragment.py**
* Extract a fragment of a fasta file (either a whole sequence, or a fragment of a sequence (1-based positions))

**maskTelomeres.py**
* Mask X kb at the start and end of each sequence of a fasta file. 

**reorderContigs.py**
* Reorder contigs based on a reference genome. Usefull to have a clean diagonal when plotting against the reference. 

**renameContigs.py**
* Rename contigs based on a reference genome. 

## Genome assembly filtering
**filterContigSize.py**
* Filter out contigs smaller than a specified threshold. 

**getNuclearContigs.py**
* Retrieve contigs corresponding the nuclear genome, based on sequence similarity with a reference genome or a BLAST database containing nuclear chromosomes (several strain avoid reference bias). 

**removeRedundantContigs.py**
* Remove contigs that are entirely covered by other contigs of the same assembly. 

## Genome assembly information
**NbContigsToXCoverage.py**
* Give the number of contig required to cover X% of a reference genome. 

**checkChromosomePresence.py**
* Take a reference genome and indicate how many reference chromosomes are covered at X % (80% default) by the draft assembly. A detailled output (-i) gives the coverage of each chromosome of the reference. 

**findDoubleCentromeres.py**
* Take the centromere sequence and indicate how many contigs contain more than 1 centromere. 

**findMergedChromosomes.py**
* Give the number of contigs covering several reference chromosomes. Can be inacurate because of large translocation. Better use findDoubleCentromere.py to identify merged chromosomes. 

