#!/bin/bash

sed '/chrM/d;/random/d;/chrUn/d' ./sams/$1.sam > ./filt_sams/$1_filtered.sam # This filters 
# mitochondrial genes + unaligned genes.
ml samtools

samtools view -S -b ./filt_sams/$1_filtered.sam > ./bams/$1.bam
# This converts filtered sam file to bam file.

samtools view -c ./bams/$1.bam # Print count of reads.
