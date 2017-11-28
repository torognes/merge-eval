#!/bin/sh

wget -nv ftp://ftp.jgi-psf.org/pub/JGI_data/Chlamy/v3.0/Chlre3.fasta.gz

gunzip Chlre3.fasta.gz

mv Chlre3.fasta chlre3.fasta
