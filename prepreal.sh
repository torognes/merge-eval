#!/bin/sh

date

THREADS=1
export OMP_NUM_THREADS=$THREADS

# SRR067577, Illumina Genome Analyzer II, 165bp frag, Broad, 101bp read len

GENOME=./Staphylococcus_aureus

XFWD=${GENOME}_1.fastq
XREV=${GENOME}_2.fastq

FWD=${GENOME}_R1.fastq
REV=${GENOME}_R2.fastq

REF=$GENOME.fasta

SAM=$GENOME.sam
LEN=$GENOME.len

if [ ! -e $REV ]; then

    echo
    echo Fixing headers
    echo

    perl ./fixfqheaders.pl $XFWD 1 > $FWD
    perl ./fixfqheaders.pl $XREV 2 > $REV
    
    echo
    date

fi

if [ ! -e $REF.bwt ]; then

    echo
    echo Indexing reference genome
    echo
    
    bwa index $REF

    echo
    date

fi

if [ ! -e $SAM ]; then

    echo
    echo Mapping reads
    echo
    
    bwa mem -t $THREADS $REF $FWD $REV > $SAM
    
    echo
    date

fi

if [ ! -e $LEN ]; then

    echo
    echo Extracting fragment lengths
    echo
    
    perl ./extractlen.pl < $SAM > $LEN
    
    echo
    date

fi

echo
echo Done
