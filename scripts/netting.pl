#! /usr/bin/perl
#
# run lastz jobs for two lists of fasta files
#

use strict;
use warnings;
use Getopt::Std;
use File::Basename;

use lib "/home/qczhang/lib/perllib";
use Locale::Schedule::Simple qw( &remoteCommand );

use vars qw ($opt_h $opt_V $opt_D $opt_d $opt_a $opt_n $opt_o $opt_1 $opt_2 $opt_l );
&getopts('hVDd:a:n:o:1:2:l:');

my $ucsc_tools= "/srv/gs1/software/ucsc_tools/2.7.2/bin/x86_64/";

my $usage = <<_EOH_;
Chain psl alignments into net files

Command:
$0 -d chain_directory -a all_chain_file -n preNet_file -o over_chain_file

# what it is:
 -d   directory of chain files to be netted
 -a   file of all chains
 -n   noClass net file
 -o   final over chain file 

 -1   target size file
 -2   query size file

# example
./netting.pl -d results/hg19_pt4/chain/ -a results/hg19_pt4/hg19ToPanTro4.all.chain -n results/hg19_pt4/hg19ToPanTro4.preNet -o results/hg19_pt4/hg19ToPanTro4.over.chain -1 ~/lib/fasta/hg19.sizes -2 ~/lib/fasta/panTro4.sizes

_EOH_
    ;

&main();

sub main 
{
    my %parameters = &init();

    open ( LOG, ">$parameters{logFile}" );
    print STDERR "collect all chain files...\n";
    print LOG "collect all chain files...\n";
    print LOG "find $parameters{chainDir} -name \"*.chain\" | $ucsc_tools/chainMergeSort -inputList=stdin > $parameters{allChainFile}\n";
    print STDERR `find $parameters{chainDir} -name "*.chain" | $ucsc_tools/chainMergeSort -inputList=stdin > $parameters{allChainFile}`;
    die ( "$parameters{allChainFile} does not exist. Chainning not succesful?\n" ) if ( ( not -e $parameters{allChainFile} ) || ( -z $parameters{allChainFile} ) );

    print STDERR "netting from file $parameters{allChainFile}...\n";
    print LOG "netting from file $parameters{allChainFile}...\n";
    print LOG "$ucsc_tools/chainPreNet $parameters{allChainFile} $parameters{targetSizes} $parameters{querySizes} stdout | $ucsc_tools/chainNet stdin -minSpace=1 $parameters{targetSizes} $parameters{querySizes} stdout /dev/null | $ucsc_tools/netSyntenic stdin $parameters{noClassNetFile}\n";
    print STDERR `$ucsc_tools/chainPreNet $parameters{allChainFile} $parameters{targetSizes} $parameters{querySizes} stdout | $ucsc_tools/chainNet stdin -minSpace=1 $parameters{targetSizes} $parameters{querySizes} stdout /dev/null | $ucsc_tools/netSyntenic stdin $parameters{noClassNetFile}`;
    die ( "$parameters{noClassNetFile} does not exist. Netting not succesful?\n" ) if ( ( not -e $parameters{noClassNetFile} ) || ( -z $parameters{noClassNetFile} ) );

    print STDERR "finalizing all.over chain file...\n";
    print LOG "finalizing all.over chain file...\n";
    print LOG "$ucsc_tools/netChainSubset -verbose=0 $parameters{noClassNetFile} $parameters{allChainFile} stdout | $ucsc_tools/chainStitchId stdin $parameters{overChainFile}\n";
    print STDERR `$ucsc_tools/netChainSubset -verbose=0 $parameters{noClassNetFile} $parameters{allChainFile} stdout | $ucsc_tools/chainStitchId stdin $parameters{overChainFile}`;
    die ( "$parameters{overChainFile} does not exist. Overchainning not succesful?\n" ) if ( ( not -e $parameters{overChainFile} ) || ( -z $parameters{overChainFile} ) );
    close LOG;
}

sub init 
{
    die $usage if ( $opt_h or ( not $opt_d ) or ( not $opt_o ) );

    my %parameters = ();

    $parameters{chainDir} = $opt_d;
    $parameters{overChainFile} = $opt_o;
    $parameters{allChainFile} = $opt_a;
    $parameters{noClassNetFile} = $opt_n;

    $parameters{targetSizes} = $opt_1;
    $parameters{querySizes} = $opt_2;

    if ( defined $opt_l ) { $parameters{logFile} = $opt_l; }
    else { $parameters{logFile} = "netting.log"; }

    return ( %parameters );
}
