#!/bin/sh

#SBATCH --account=nn9383k
#SBATCH --job-name=merge
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks-per-node=4

#source /cluster/bin/jobsetup

# Note: kmerfreq in cope.1 needs 16G RAM for k=17

module purge
module load flash
module load pandaseq

THREADS=4
export OMP_NUM_THREADS=$THREADS

TIME=/usr/bin/time

BBMERGE=/projects/rhbioinf/torognes/bbmap/bbmerge.sh
CASPER=/projects/rhbioinf/torognes/vsearch-merge/casper/casper_v0.8.2/casper
COPE=/projects/rhbioinf/torognes/cope-1.2.5/cope
FLASH=$(which flash)
FQJOIN=/projects/rhbioinf/software/fastq-join-1.3.1/fastq-join
KMERFREQ=/projects/rhbioinf/torognes/cope-1.2.5/kmerfreq
LEEHOM=/projects/rhbioinf/torognes/leeHom-v.1.1.5/src/leeHomMulti
PANDASEQ=$(which pandaseq)
PEAR=$(which PEAR)
STITCH=/projects/rhbioinf/torognes/stitch/stitch/stitch.py
USEARCH7=$(which usearch)
USEARCH8=$(which usearch8)
VSEARCH=/projects/rhbioinf/torognes/vsearch/bin/vsearch
VSEARCH252=/projects/rhbioinf/torognes/vsearch/bin/vsearch-2.5.2
XORRO=/projects/rhbioinf/torognes/XORRO/xorro

#GENOME=./ecoli
#GENOME=./Staphylococcus_aureus
#GENOME=./uneven
#GENOME=./chlre3
#GENOME=strco
#GENOME=bacce

REF=$GENOME.fasta
#REF=./uneven.top25.fasta
#REF=GCF_000203835.1_ASM20383v1_genomic.fna # strco
#REF=GCF_000007825.1_ASM782v1_genomic.fna # bacce

FWD=${GENOME}_R1.fastq
REV=${GENOME}_R2.fastq
LENFILE=${GENOME}.len

if [ ! "$1" == "" ]; then
    PROGS="$1"
else
    PROGS="bbmerge casper cope.0 cope.1 cope.2 cope.3 flash fqjoin leehom pandaseq.pear pandaseq.rdp_mle pandaseq.simple_bayesian pear stitch usearch7 usearch8 vsearch vsearch252 xorro"
fi

for P in $PROGS; do

    date

    PREFIX=$GENOME.$P
    MERGED=$PREFIX.merged.fastq
    LOG=$PREFIX.log
    ANALYSIS=$PREFIX.analysis.txt

    rm -f $MERGED $LOG

    echo
    echo Running $P
    echo

    case $P in

        bbmerge)
            $TIME $BBMERGE t=$THREADS \
                in1=$FWD in2=$REV out=$MERGED \
                > $LOG 2>&1
            ;;

        casper)
            $TIME $CASPER -t $THREADS \
                $FWD $REV \
                > $LOG 2>&1
            mv casper.fastq $MERGED
            ;;

        cope.0)
            $TIME $COPE -m 0 \
                -a $FWD -b $REV -o $PREFIX \
                -s 33 \
                > $LOG 2>&1
            sed -e "s/_ART.*//g" $PREFIX.connect.fq > $MERGED
            rm $PREFIX.connect.fq $PREFIX.unConnect_*.fq
            ;;

        cope.1|cope.2|cope.3)
            
            K=15
            echo $FWD > $PREFIX.kmerfreq.list
            echo $REV >> $PREFIX.kmerfreq.list

            # Note: RAM-hungry! Adjust with k
            # Do not redirect stderr! Cannot get timing
            $TIME $KMERFREQ -t $THREADS \
                -q 33 -k $K -p $PREFIX.kmerfreq \
                $PREFIX.kmerfreq.list \
                > $LOG

            $TIME $COPE -m 1 \
                -a $FWD -b $REV -o $PREFIX \
                -s 33 -k $K \
                -t $PREFIX.kmerfreq.freq.cz \
                -f $PREFIX.kmerfreq.freq.cz.len \
                >> $LOG 2>&1

            sed -e "s/_ART.*//g" $PREFIX.connect.fq > $MERGED
            rm $PREFIX.connect.fq $PREFIX.unConnect_{1,2}.fq
            rm $PREFIX.kmerfreq.* pair.kmer.list.gz
            ;;

        flash)
            $TIME $FLASH -t $THREADS \
                -q $FWD $REV \
                -M 500 \
                > $LOG 2>&1
            mv out.extendedFrags.fastq $MERGED
            rm out.hist out.histogram
            rm out.notCombined_1.fastq out.notCombined_2.fastq
            ;;
        
        fqjoin)
            $TIME $FQJOIN $FWD $REV \
                -o $PREFIX. \
                > $LOG 2>&1
            mv $PREFIX.join $MERGED
            rm $PREFIX.un{1,2}
            ;;

        leehom)
            $TIME $LEEHOM --ancientdna -t $THREADS \
                -fq1 $FWD -fq2 $REV -fqo $PREFIX \
                > $LOG 2>&1
            gzip -dc $PREFIX.fq.gz > $MERGED
            rm $PREFIX*.fq.gz
            ;;

        pandaseq.pear)
            $TIME $PANDASEQ -A pear -T $THREADS \
                -f $FWD -r $REV -w $MERGED \
                -B -F -g /dev/null \
                > $LOG 2>&1
            ;;

        pandaseq.rdp_mle)
            $TIME $PANDASEQ -A rdp_mle -T $THREADS \
                -f $FWD -r $REV -w $MERGED \
                -B -F -g /dev/null \
                > $LOG 2>&1
            ;;

        pandaseq.simple_bayesian)
            $TIME $PANDASEQ -A simple_bayesian -T $THREADS \
                -f $FWD -r $REV -w $MERGED \
                -B -F -g /dev/null \
                > $LOG 2>&1
            ;;

        pear)
            $TIME $PEAR --threads $THREADS \
                --forward-fastq $FWD --reverse-fastq $REV --output $PREFIX \
                > $LOG 2>&1
            mv $PREFIX.assembled.fastq $MERGED
            rm $PREFIX.unassembled.* $PREFIX.discarded.fastq
            ;;

        stitch)
            $TIME $STITCH \
                -i $FWD -j $REV -o $PREFIX \
                > $LOG 2>&1
            mv $PREFIX-contigs.fastq $MERGED
            rm $PREFIX-nh-s?.fastq
            ;;


        usearch7)
            $TIME $USEARCH7 --threads $THREADS \
                --fastq_mergepairs $FWD --reverse $REV --fastqout $MERGED \
                > $LOG 2>&1
            ;;

        usearch8)
            $TIME $USEARCH8 --threads $THREADS \
                --fastq_mergepairs $FWD --reverse $REV --fastqout $MERGED \
                > $LOG 2>&1
            ;;

        vsearch)
            $TIME $VSEARCH --threads $THREADS \
                --fastq_mergepairs $FWD --reverse $REV --fastqout $MERGED \
                > $LOG 2>&1
            ;;
        
        vsearch252)
            $TIME $VSEARCH252 --threads $THREADS \
                --fastq_mergepairs $FWD --reverse $REV --fastqout $MERGED \
                > $LOG 2>&1
            ;;

        xorro)
            $TIME $XORRO \
                $FWD $REV 1000 10 33 \
                > $MERGED \
                2> $LOG
            ;;

        *)
            echo Unknown program $P
            echo
            ;;

    esac

    echo Timing:
    echo
    grep -A 1 maxresident $LOG
    echo

    echo Analysis of results from $P:
    echo

    ./checkcorrect.pl $LENFILE < $MERGED > $ANALYSIS

    grep "^Pairs with serious errors:" $ANALYSIS
    grep "^Pairs merged, should not have been:" $ANALYSIS
    grep "^Pairs merged, incorrectly (FP):" $ANALYSIS
    grep "^Pairs not merged, incorrectly (FN):" $ANALYSIS
    grep "^Recall (TP/(TP+FN)):" $ANALYSIS
    grep "^Precision (TP/(TP+FP)):" $ANALYSIS
    
    echo

done

echo
echo Done
echo

date
