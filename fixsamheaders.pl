#!/usr/bin/perl -w

# reformat sam headers due to the strictness of pandaseq

my $pair = 0;
my $last = "";

while(<>)
{
    next if /^@/;

    if (/^([^\t]*)\t(.*)$/)
    {
        $pair++ if $1 ne $last;
        $last = $1;
        print "ART:1:1:1:1:1:$pair\t$2\n";
    }
}
