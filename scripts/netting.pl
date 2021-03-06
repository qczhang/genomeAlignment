#! /usr/bin/perl
#
# run lastz jobs for two lists of fasta files
#

use strict;
use warnings;
use Getopt::Std;

use vars qw ($opt_h $opt_V $opt_D $opt_1 $opt_2 $opt_d $opt_o $opt_m $opt_p $opt_n $opt_s $opt_c $opt_t $opt_l $opt_a );
&getopts('hVD1:2:d:o:m:p:n:s:c:t:l:a:');

my $usage = <<_EOH_;
Chain psl alignments into net files

Command:
$0 -1 target -2 query -d chain_directory -o output_dir

# what it is:
 -1   target size file
 -2   query size file

 -d   directory of chain files to be netted
 -o   output directory

# more options
 -l   log file

 -m   chainMergeSort binary
 -p   chainPreNet binary
 -n   chainNet binary
 -s   netSyntenic binary
 -c   netChainSubset binary
 -t   chainStitchId binary

 -a   input all chain file

# example
./netting.pl -1 ~/lib/fasta/hg19.sizes -2 ~/lib/fasta/panTro4.sizes -d results/hg19_pt4/chain/ -o results/hg19_pt4/

_EOH_
    ;

&main();

sub main 
{
    my %parameters = &init();

    open ( LOG, ">$parameters{logFile}" );
    my $allChainFile = "$parameters{outDir}/all.chain";
    if ( defined $parameters{allChain} ) { $allChainFile = $parameters{allChain}; }
    else {
        print STDERR "Collect all chain files...\n\t", `date`;
        print LOG "Collect all chain files...\n\t", `date`;
        print LOG "find $parameters{chainDir} -name \"*.chain\" | $parameters{mergeSortBin} -inputList=stdin > $allChainFile\n";
        print STDERR `find $parameters{chainDir} -name "*.chain" | $parameters{mergeSortBin} -inputList=stdin > $allChainFile`;
        die ( "$allChainFile does not exist. Chainning not succesful?\n" ) if ( ( not -e $allChainFile ) || ( -z $allChainFile ) );
    }

    my $noClassNetFile = "$parameters{outDir}/all.preNet";
    print STDERR "\nNetting from file $allChainFile...\n\t", `date`;
    print LOG "\nNetting from file $allChainFile...\n\t", `date`;
    print LOG "$parameters{preNetBin} $allChainFile $parameters{targetSizes} $parameters{querySizes} stdout | $parameters{netBin} stdin -minSpace=1 $parameters{targetSizes} $parameters{querySizes} stdout /dev/null | $parameters{netSynBin} stdin $noClassNetFile\n";
    print LOG `$parameters{preNetBin} $allChainFile $parameters{targetSizes} $parameters{querySizes} stdout | $parameters{netBin} stdin -minSpace=1 $parameters{targetSizes} $parameters{querySizes} stdout /dev/null | $parameters{netSynBin} stdin $noClassNetFile`;
    die ( "$noClassNetFile does not exist. Netting not succesful?\n" ) if ( ( not -e $noClassNetFile ) || ( -z $noClassNetFile ) );

    my $overChainFile = "$parameters{outDir}/over.chain";
    print STDERR "\nFinalizing over chain file...\n\t", `date`;
    print LOG "\nFinalizing over chain file...\n\t", `date`;
    print LOG "$parameters{netChainBin} -verbose=0 $noClassNetFile $allChainFile stdout | $parameters{stitchBin} stdin $overChainFile\n";
    print LOG `$parameters{netChainBin} -verbose=0 $noClassNetFile $allChainFile stdout | $parameters{stitchBin} stdin $overChainFile`;
    die ( "$overChainFile does not exist. Overchainning not succesful?\n" ) if ( ( not -e $overChainFile ) || ( -z $overChainFile ) );
    close LOG;
}

sub init 
{
    die $usage if ( $opt_h or ( ( not $opt_d ) and ( not $opt_a ) ) or ( not $opt_o ) );

    my %parameters = ();

    $parameters{targetSizes} = $opt_1;
    $parameters{querySizes} = $opt_2;
    $parameters{outDir} = $opt_o;
    if ( defined $opt_d ) { $parameters{chainDir} = $opt_d; }
    if ( defined $opt_a ) { $parameters{allChain} = $opt_a; }

    if ( defined $opt_l ) { $parameters{logFile} = $opt_l; }
    else { $parameters{logFile} = "$parameters{outDir}/netting.log"; }

    if ( defined $opt_m ) { $parameters{mergeSortBin} = $opt_m; }
    else { $parameters{mergeSortBin} = "chainMergeSort";  } ## check availability of binary chainMergeSort
    if ( defined $opt_p ) { $parameters{preNetBin} = $opt_p; }
    else { $parameters{preNetBin} = "chainPreNet";  } ## check availability of binary chainPreNet
    if ( defined $opt_n ) { $parameters{netBin} = $opt_n; }
    else { $parameters{netBin} = "chainNet";  } ## check availability of binary chainNet
    if ( defined $opt_s ) { $parameters{netSynBin} = $opt_s; }
    else { $parameters{netSynBin} = "netSyntenic";  } ## check availability of binary netSyntenic
    if ( defined $opt_c ) { $parameters{netChainBin} = $opt_c; }
    else { $parameters{netChainBin} = "netChainSubset";  } ## check availability of binary netChainSubset
    if ( defined $opt_t ) { $parameters{stitchBin} = $opt_t; }
    else { $parameters{stitchBin} = "chainStitchId";  } ## check availability of binary chainStitchId

    return ( %parameters );
}
