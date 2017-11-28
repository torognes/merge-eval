#!/usr/bin/perl -w

# reformat headers due to the strictness of pandaseq

my ($filename, $dir) = @ARGV;

my $line = 1;
my $seq = 1;

open F, $filename;

while(<F>)
{
    if ($line % 4 == 1)
    {
        print "\@ART:1:1:1:1:1:$seq $dir:N:0:\n";
        $seq++;
    }
    else
    {
        print $_;
    }
    $line++;
}

close F;
