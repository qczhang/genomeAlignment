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

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_d $opt_1 $opt_2 $opt_o $opt_p $opt_l );
&getopts('hVDi:1:2:d:o:p:l:');

my $axtChainBin = "/srv/gs1/software/ucsc_tools/2.7.2/bin/x86_64/axtChain";
my $chainAntiRepeatBin = "/srv/gs1/software/ucsc_tools/2.7.2/bin/x86_64/chainAntiRepeat";

my $usage = <<_EOH_;
Chain psl alignments into net files

Command:
$0 -i target_list -1 target_2bit -2 query_2bit -o output_directory

# what it is:
 -i   target list
 -1   target 2bit file
 -2   query 2bit file
 -o   output directory of chain file

# more options
 -p   options in running axtChain
 -d   directory of target alignment psl files

# example
 ./chainning.pl -i ~/lib/fasta/hg19.sizes -1 ~/lib/fasta/hg19.2bit -2 ~/lib/fasta/panTro4.2bit -o results/hg19_pt4/ -d results/hg19_pt4/psl/

_EOH_
    ;

&main();

sub main 
{
    my ( %parameters ) = &init();

    #my $ref_targetList = readFastaList ( $parameters{targetList}, skip => "_" );
    my $ref_targetList = readFastaList ( $parameters{targetList} );

    open ( LOG, ">$parameters{logFile}" );
    foreach my $targetName ( keys ( %{$ref_targetList} ) ) {
        my $inputLav = $parameters{pslDir} . "/" . $targetName . ".psl";
        my $outputChain = $parameters{outputDir} . "/chain/" . $targetName . ".chain";
        my $chainningCmd = "cat $inputLav";
        $chainningCmd .= " | $axtChainBin $parameters{axtChainOptions} stdin $parameters{target2bit} $parameters{query2bit} stdout";
        $chainningCmd .= " | $chainAntiRepeatBin $parameters{target2bit} $parameters{query2bit} stdin $outputChain\n";

        if ( not $parameters{remote} ) {
            print LOG "\n$chainningCmd\n\tTime: ", `date`;
            print LOG `$chainningCmd`;
            if ( not $? ) { print LOG "Job successful, check $outputChain for output chain. \n"; }
            else {
                print STDERR "Warning! Job not successful, $outputChain not finished.\n";
                print LOG "Warning! Job not successful, $outputChain not finished.\n";
            }
        }
        else {
            $chainningCmd .= "\ntouch $outputChain.done\n";
            my $chainningScriptFile = $parameters{outputDir} . "/chainning_shell/" . $targetName . ".sh";
            my $chainningJobSubmitted = $parameters{outputDir} . "/chainning_shell/" . $targetName . ".sh.submitted";
            my $chainningJobName = "$targetName.chainning";  my $chainningJobMemory = "4G";  my $chainningJobTime = "3:00:00";
            remoteCommand( $chainningJobName, $chainningJobMemory, $chainningJobTime, $chainningScriptFile, $chainningCmd, $chainningJobSubmitted );
        }
    }
}

sub init 
{
    die $usage if ( $opt_h or ( not $opt_i ) or ( not $opt_1 ) or ( not $opt_2 ) or ( not $opt_o ) );

    my %parameters = ();

    $parameters{targetList} = $opt_i;
    $parameters{outputDir} = $opt_o;
    mkdir ( $parameters{outputDir} ) if ( not -e $parameters{outputDir} );
    mkdir ( "$parameters{outputDir}/chain" ) if ( not -e "$parameters{outputDir}/chain" );
    $parameters{target2bit} = $opt_1;
    $parameters{query2bit} = $opt_2;

    if ( defined $opt_l ) { $parameters{logFile} = $opt_l; }
    else { $parameters{logFile} = "chainning.log"; }
    if ( defined $opt_p ) { $parameters{axtChainOptions} = $opt_p; }
    else { $parameters{axtChainOptions} = "-psl -verbose-0 -minScore=5000 -linearGap=medium"; }
    if ( defined $opt_d ) { $parameters{pslDir} = $opt_d; }
    $parameters{remote} = 0;
    #mkdir ( "$parameters{outputDir}/chainning_shell" ) if ( not -e "$parameters{outputDir}/chainning_shell" );

    return ( %parameters );
}

sub readFastaList {
    my $listFile = shift;
    my %parameters = @_;

    my %fasta = ();
    open ( LIST, $listFile ) or die "Cannot open $listFile for reading!\n";
    while ( my $line = <LIST> ) {
        next if ( $line =~ /^#/ );
        chomp $line;

        my @fields = split ( /\t/, $line );
        if ( defined $parameters{skip} ) { next if ( $fields[0] =~ /$parameters{skip}/ ); }

        $fasta{$fields[0]} = 1;
    }
    close LIST;

    return \%fasta;
}
