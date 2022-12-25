#!/bin/bash
ml python

inp="bams/"
ext=".bam"
out="macs/"

macs2 callpeak -t "${inp}SRR502327${ext}" -c "${inp}SRR502225_11E6${ext}" -n "STAT1_6h_IFNa" --outdir $out -g hs --bdg -q 0.05 -f BAM
echo "Done for first. Contents are $(ls \macs)"
macs2 callpeak -t "${inp}SRR502329${ext}" -c "${inp}SRR502228_11E6${ext}" -n "STAT1_30m_IFNa" --outdir $out -g hs --bdg -q 0.05 -f BAM
echo "Done for second. Contents are $(ls \macs)"

# -g homo sapiens(hs), --bdp bedgraph output, -q 0.05 q-value -f bAM file input format 
