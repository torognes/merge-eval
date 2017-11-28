#!/usr/bin/perl -w

use strict;
use warnings;

my $maxlen = 1000;
my $minmapqual = 60;

my %pairlength;
my %readlength;

my $count = 0;
my $total = 0;

my $last = "";

while (<>)
{
    next if /^@/;
    my $line = $_;
    my @col = split /\t/;
    my $query = $col[0];
    my $flags = $col[1];
    my $mapqual = $col[4];
    my $cigar = $col[5];
    my $fraglength = $col[8];
    my $seq = $col[9];

    my $readlen = length $seq;

    $total++ if $query ne $last;

    # perfect mapping quality (>= 60)
    # both ends of pair mapped properly
    # primary alignment
    # no QC failure
    # no duplicate
    # no supplementary alignment
    # reasonable insert length (<= 1000)
    # no soft clipping (S)
    # no hard clipping (H)

    if (($mapqual >= $minmapqual)
        && (($flags & 0x003) == 0x003)
        && (($flags & 0xF0C) == 0x000)
        && (abs($fraglength) <= $maxlen)
        && (!($cigar =~ /S/))
        && (!($cigar =~ /H/)))
    {
        my $plen = $pairlength{$query};
        my $rlen = $readlength{$query};
        
        if (defined $plen)
        {
            if ($plen == - $fraglength)
            {
                printf "$query\t%d\t%d\n", abs($fraglength), $readlen + $rlen;
                $count++;
            }
            else
            {
                printf STDERR "Inconsistent fragment size for $query\n";
            }
        }
        else
        {
            $pairlength{$query} = $fraglength;
            $readlength{$query} = $readlen;
        }
    }

    $last = $query;
}

print STDERR "Found fragment lengths for $count pairs out of $total in total.\n";
