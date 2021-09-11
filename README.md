## Introduction

tagbam is a C program that accepts a bam file with sequences from amplicon-based sequencing (eg, Haloplex) and a bed with amplicon intervals and creates a new bam file with the name of the best-fitting amplicon in the XN bam tag field. Bam file manipulation is done via the htslib library and overlaps between reads and bed intervals uses Heng Li's  cgranges library.

## Dependencies and Installation

The htslib must be installed. tagbam can be compiled via:

gcc tagbam.c cgranges.c -o tagbam -lhts -lz

## Usage

tagbam [options] <input.bam> <amplicon.bed> <output.bam>
Options:
  -d <int>   maximum sum of the distance bt amplicon and read ends (default=6)
  -v         verbose
	
### Test data

A small test dataset is included, which can be tested via:
	
   tagbam sub.bam amps.bed out.bam 


https://github.com/samtools/htslib
https://github.com/lh3/cgranges/
