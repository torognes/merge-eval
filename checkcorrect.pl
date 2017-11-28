#!/usr/bin/perl -w

use strict;
use warnings;

my $badmax = 10;
my $badmergemax = 10;
my $badmissmax = 10;

# read correct fragment lengths

my ($lengthfilename) = @ARGV;

my %tworeadslen = ();
my %mergelen = ();
my $all = 0;
my $possible = 0;
open F, $lengthfilename;
my %fraglen;
while (<F>)
{
    /^(.*)\t(\d+)\t(\d+)$/;
    my ($query, $length, $maxfraglen) = ($1, $2, $3);
    $fraglen{$query} = $length;
    $tworeadslen{$query} = $maxfraglen;
    $all++;
    $possible++ if $length < $maxfraglen;
    $mergelen{$query} = 0;
}
close F;

my $line = 1;
my $query = "";
my $sequence = "";

my $TP = 0;
my $FP = 0;
my $TN = 0;
my $FN = 0;

while (<STDIN>)
{
    chomp;
    if ($line == 1)
    {
        /^@([^ \t\/]+)/;
        $query = $1;
        if ($query =~ /(.*):$/)
        {
            $query = $1;
        }
    }
    elsif ($line == 2)
    {
        my $length = (length $_);
        $mergelen{$query} = $length;
    }

    $line++;
    if ($line == 5)
    {
        $query = "";
        $line = 1;
    }
}

my $FP_mergeable = 0;
my $FP_nonmergeable = 0;
my $FP_indel = 0;
my $bad = 0;
my $badmerge = 0;
my $badmiss = 0;

my $del = 0;
my $ins = 0;

my $del1 = 0;
my $ins1 = 0;

my $elim = 0;

for (keys %mergelen)
{
    my $maxfraglen = $tworeadslen{$_};

    if (! defined $fraglen{$_})
    {
        $elim++;
    }
    else
    {
        if ($mergelen{$_} > 0)
        {
            if ($mergelen{$_} == $fraglen{$_})
            {
                $TP++;
            }
            else
            {
                $FP++;
                if ($fraglen{$_} < $maxfraglen)
                {
                    
                    $FP_mergeable++;
                    if (abs($mergelen{$_} - $fraglen{$_}) <= 2)
                    {
                        $FP_indel++;
                        if ($mergelen{$_} > $fraglen{$_})
                        {
                            $ins1++;
                        }
                        else
                        {
                            $del1++;
                        }
                    }
                    else
                    {
                        
                        printf "BAD: $_ Merged length=%d, Fragment length=%d Delta=%d\n", $mergelen{$_}, $fraglen{$_}, $mergelen{$_} - $fraglen{$_} if $badmerge < $badmergemax;
                        $badmerge++;
                        
                        if ($mergelen{$_} > $fraglen{$_})
                        {
                            $ins++;
                        }
                        else
                        {
                            $del++;
                        }
                    }
                }
                else
                {
                    printf "BAD: $_ Merged length=%d, should not have been merged\n", $mergelen{$_} if $bad < $badmax;
                    $bad++;
                    $FP_nonmergeable++;
                }
            }
        }
        else
        {
            if ($fraglen{$_} < $maxfraglen)
            {
                printf "BAD: $_ Not merged, Fragment length=%d\n", $fraglen{$_} if $badmiss < $badmissmax;
                $badmiss++;
                $FN++;
            }
            else
            {
                $TN++;
            }
        }
    }
}

# TP, mergeable, merged correctly
# TN, not mergeable, not merged
# FP, incorrectly merged (either mergeable or not)
# FN, mergeable, not merged

my $correct = $TP + $TN;
my $incorrect = $FP + $FN;
my $merged = $TP + $FP;
my $recall = $TP / ($TP + $FN);
my $precision = $TP + $FP == 0 ? 0.0 : $TP / ($TP + $FP);
my $accuracy = ($TP+$TN)/($TP+$FN+$FP+$TN);
my $F1 = 2*$TP/(2*$TP+$FP+$FN);
my $notmerged = $all - $merged;

printf "\n";

printf "Eliminated (merged but not mapped):  %8d\n", $elim;

printf "Total pairs:                         %8d\n", $all;
printf "Pairs merged:                        %8d (%9.5f%%)\n", $merged, 100.0 * $merged / $all;
printf "Pairs not merged:                    %8d (%9.5f%%)\n", $notmerged, 100.0 * $notmerged / $all;
printf "\n";
printf "Total correct:                       %8d (%9.5f%%)\n", $correct, 100.0 * $correct / $all;
printf "Total incorrect:                     %8d (%9.5f%%)\n", $incorrect, 100.0 * $incorrect / $all;
printf "\n";
printf "Pairs merged, correctly (TP):        %8d (%9.5f%%)\n", $TP, 100.0 * $TP / $all;
printf "Pairs merged, incorrectly (FP):      %8d (%9.5f%%)\n", $FP, 100.0 * $FP / $all;
printf "Pairs merged with wrong insert size: %8d (%9.5f%%)\n", $FP_mergeable, 100.0 * $FP_mergeable / $all;
printf "Pairs merged with 1-2bp indel:       %8d (%9.5f%%)\n", $FP_indel, 100.0 * $FP_indel / $all;
printf "Insertions >2bp:                     %8d (%9.5f%%)\n", $ins, 100.0 * $ins / $all;
printf "Deletions >2bp:                      %8d (%9.5f%%)\n", $del, 100.0 * $del / $all;
printf "Insertions 1-2bp:                    %8d (%9.5f%%)\n", $ins1, 100.0 * $ins1 / $all;
printf "Deletions 1-2bp:                     %8d (%9.5f%%)\n", $del1, 100.0 * $del1 / $all;
printf "Pairs merged, should not have been:  %8d (%9.5f%%)\n", $FP_nonmergeable, 100.0 * $FP_nonmergeable / $all;
printf "\n";
printf "Pairs not merged, correctly (TN):    %8d (%9.5f%%)\n", $TN, 100.0 * $TN / $all;
printf "Pairs not merged, incorrectly (FN):  %8d (%9.5f%%)\n", $FN, 100.0 * $FN / $all;
printf "\n";
printf "Pairs with serious errors:           %8d (%9.5f%%)\n", ($FP-$ins1-$del1), 100.0 * ($FP-$ins1-$del1) / $all;
printf "\n";
printf "Recall (TP/(TP+FN)):                 %.6f\n", $TP/($TP+$FN);
printf "Precision (TP/(TP+FP)):              %.6f\n", $precision;
printf "Accuracy:                            %.6f\n", $accuracy;
printf "F1 score:                            %.6f\n", $F1;
