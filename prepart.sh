#!/bin/sh

date

PREFIX=./ecoli
#PREFIX=./Staphylococcus_aureus
#PREFIX=./uneven
#PREFIX=./chlre3
#PREFIX=./staph
#PREFIX=strco
#PREFIX=bacce

REF=$PREFIX.fasta
#REF=./uneven.top25.fasta
#REF=GCF_000203835.1_ASM20383v1_genomic.fna #strco
#REF=GCF_000007825.1_ASM782v1_genomic.fna #bacce

READ_LENGTH=250

COVERAGE=10

FRAGMENT_MEAN=450
FRAGMENT_STDEV=30

SEQUENCER=MSv3
#SEQUENCER=HS25
#SEQUENCER=GA2
#SEQUENCER=HSXn

#RANDOMSEED=893795873

FWD=${PREFIX}_R1.fastq
REV=${PREFIX}_R2.fastq

echo
echo Generate synthetic reads
echo

ART=/projects/rhbioinf/software/art_bin_MountRainier/art_illumina

INDELS="" # default
#INDELS="-ir 0.00009 -ir2 0.00015 -dr 0.00011 -dr2 0.00023" # default
#INDELS="-ir 0.0009  -ir2 0.0015  -dr 0.0011  -dr2 0.0023"  # 10X
#INDELS="-k 0"                                              # no indels

# Generate whole-genome paired-end reads
$ART -sam -p -ss $SEQUENCER -l $READ_LENGTH -f $COVERAGE -m $FRAGMENT_MEAN -s $FRAGMENT_STDEV -i $REF -o $PREFIX -sp -na $INDELS # -rs $RANDOMSEED

# Generate amplicon paired-end reads
#$ART -amp -sam -p -ss $SEQUENCER -l $READ_LENGTH -f $COVERAGE -i $REF -o $PREFIX -sp -na $INDELS # -rs $RANDOMSEED

echo
echo Reformat FASTQ and SAM files, extract lengths
echo

perl ./fixfqheaders.pl ${PREFIX}1.fq 1 > $FWD
perl ./fixfqheaders.pl ${PREFIX}2.fq 2 > $REV
perl ./fixsamheaders.pl < $PREFIX.sam > temp.sam ; mv temp.sam $PREFIX.sam
perl ./extractlen.pl < $PREFIX.sam > $PREFIX.len

rm ${PREFIX}1.fq ${PREFIX}2.fq

date
