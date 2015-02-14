#! /usr/bin/perl -w
#
use strict;
use warnings;

my $genome = shift;
my $outputDir = shift;
if ( not $genome ) { die "Usage: $0 genome <outputDir>\n"; }

my $writeChr = 0;
if ( $outputDir ) {
    if ( not -e $outputDir ) { print STDERR `mkdir $outputDir`; }
    if ( ( -e $outputDir ) and ( -w $outputDir ) ) { $writeChr = 1; }
    else { print STDERR "Cannot open $outputDir for writing! will not output fasta files.\n"; }
}

my $count = 0; my $chr = ""; my $len = 0;  my $first = 1;  my $outFile = "NULL";
open ( GN, $genome );
while ( my $line = <GN> ) {
    $count++;
    #print STDERR $count, "\n" if ( $count % 1000000 == 0 );

    if ( $line =~ /^>/ ) {
        print $chr, "\t", $len, "\n" if ( $len );

        ( $chr ) = ( $line =~ /^>(\S+)\s/ );
        $len = 0;

        if ( $writeChr ) {
            if ( not $first ) {
                close OUT;
                $first = 0;
            }
            $outFile = $outputDir . "/" . $chr . ".fa";
            print STDERR $chr, "\t", $outFile, "\n";
            open ( OUT, ">$outFile" );
            print OUT $line;
        }
    }
    else {  
        $len += ( length($line) - 1 ); 
        if ( $writeChr ) {
            print OUT $line;
        }
    }
}
close OUT;
print $chr, "\t", $len, "\n" if ( $len );
close GN;
