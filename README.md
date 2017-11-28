# merge-eval

Here are a set of scripts to perform evaluation of paired-end read mergers.

Download genome sequences with `download.sh`. Prepare artificial reads with `prepart.sh` or real reads with `prepreal.sh`. Then perform merging and analysis with `merge.sh`. The `merge.sh` script may be called with the name of one tool to analyse only that one, otherwise all tools will be called and analysed. The scripts need to be modified to adjust the names of the genome files, the read lengths etc.

All merging tools need to be installed first. Currently 12 different tools are supported: bbmerge, casper, cope, flash, fastq-join, leehom, pandaseq, pear, stitch, usearch, vsearch and xorro. Two different versions of usearch and vsearch will be run as well as four modes of cope and three variants of pandaseq. Add the paths to the tools in the `merge.sh` file.
