#! /usr/bin/perl
#
# run lastz jobs for two lists of fasta files
#

use strict;
use warnings;
use Getopt::Std;

use vars qw ($opt_h $opt_V $opt_D $opt_1 $opt_2 $opt_o $opt_p $opt_r $opt_l $opt_z $opt_s );
&getopts('hVD1:2:o:p:l:z:s:r');

my $usage = <<_EOH_;
Run lastz for two genomes...

Command:
$0 -1 fasta_list_1 -2 fasta_list_2 -o output_directory

# what it is:
 -1   input list 1 of fasta files to be aligned
 -2   input list 2 of fasta files to be aligned
 -o   output directory of LastZ output .lav alignment files

# more options
 -z   lastz binary
 -s   lav to psl binary
 -p   options in running Lastz

 -l   log file
 -r   remote running (not yet)

# example
 ./doLastz.pl -1 ~/lib/fasta/hg19.file -2 ~/lib/fasta/panTro4.file -o results/hg19_pt4/

_EOH_
;

&main();

sub main 
{
    my ( %parameters ) = &init();

    my $ref_fastaList1 = readFastaList ( $parameters{fastaList1} );
    my $ref_fastaList2 = readFastaList ( $parameters{fastaList2} );

    open ( LOG, ">$parameters{logFile}" );
    print STDERR "Running lastz jobs using options: \n\t$parameters{lastzOptions}\n";
    print LOG "Running lastz jobs using options: \n\t$parameters{lastzOptions}\n";

    my $countTotalPercent = scalar ( keys %{$ref_fastaList1} );
    my $count = 0; my $lastPerc = 0;
    foreach my $fastaName1 ( sort { $a cmp $b } ( keys ( %{$ref_fastaList1} ) ) ) {
        next if ( $fastaName1 eq "genome" );

        my $fasta1 = $ref_fastaList1->{$fastaName1}; 
        if ( not -e $fasta1 ) { print STDERR "Warning! fasta1 $fasta1 does not exist!\n"; next; }
        my $fasta2 = $ref_fastaList2->{genome}; 
        if ( not -e $fasta2 ) { print STDERR "Warning! fasta2 $fasta2 does not exist!\n"; next; }

        $count++;
        if ( int ( 100*$count/$countTotalPercent ) > $lastPerc ) {
            print STDERR "Job progress: ", int ( 100*$count/$countTotalPercent), "\%...Time: ", `date`;
            $lastPerc = int( 100*$count/$countTotalPercent );
        }

        my $outputLav = $parameters{outputDir} . "/lav/" . $fastaName1 . ".lav";
        my $outputPsl = $parameters{outputDir} . "/psl/" . $fastaName1 . ".psl";
        my $lastzCmd = "$parameters{lastzBin} $fasta1 $fasta2 $parameters{lastzOptions} > $outputLav";
        if ( not $parameters{remote} ) {
            print LOG "\n$lastzCmd\n\tTime: ", `date`;
            die;
            print LOG `$lastzCmd`;
            if ( not $? ) { 
                print LOG "Job successful, check $outputLav for output alignment. \n"; 
                print LOG "$parameters{lavToPslBin} $outputLav $outputPsl\n\tTime: ", `date`; 
                print LOG `$parameters{lavToPslBin} $outputLav $outputPsl`; 
                if ( not $? ) { print LOG "Job successful, check $outputPsl for output alignment. \n"; }
                else { 
                    print STDERR "Warning! Job not successful, $outputPsl not finished.\n";
                    print LOG "Warning! Job not successful, $outputPsl not finished.\n";
                }
            }
            else { 
                print STDERR "Warning! Job not successful, $outputLav not finished.\n";
                print LOG "Warning! Job not successful, $outputLav not finished.\n";
            }
        }
        else {
            $lastzCmd .= "\ntouch $outputLav.done\n\n";
            $lastzCmd .= "\n$parameters{lavToPslBin} $outputLav $outputPsl";
            $lastzCmd .= "\ntouch $outputLav.done\n";
            my $lastzJobSubmitted = $parameters{outputDir} . "/lastz_shell/" . $fastaName1 . ".sh.submitted";
            my $lastzJobName = "$fastaName1.lastz";  my $lastzJobMemory = "4G";  my $lastzJobTime = "3:00:00";
            my $lastzScriptFile = $parameters{outputDir} . "/lastz_shell/" . $fastaName1 . ".sh";
            #remoteCommand( $lastzJobName, $lastzJobMemory, $lastzJobTime, $lastzScriptFile, $lastzCmd, $lastzJobSubmitted );
        }
    }
}

sub init 
{
    die $usage if ( $opt_h or ( not $opt_1 ) or ( not $opt_2 ) or ( not $opt_o ) );

    my %parameters = ();
    $parameters{fastaList1} = $opt_1;
    $parameters{fastaList2} = $opt_2;
    $parameters{outputDir} = $opt_o;
    mkdir ( $parameters{outputDir} ) if ( not -e $parameters{outputDir} );
    mkdir ( "$parameters{outputDir}/lav" ) if ( not -e "$parameters{outputDir}/lav" );
    mkdir ( "$parameters{outputDir}/psl" ) if ( not -e "$parameters{outputDir}/psl" );

    if ( defined $opt_z ) { $parameters{lastzBin} = $opt_z; }
    else { $parameters{lastzBin} = "lastz";  } ## check availability of binary lastz
    if ( defined $opt_s ) { $parameters{lavToPslBin} = $opt_s; }
    else { $parameters{lastzBin} = "lavToPsl";  } ## check availability of binary lavToPsl
    if ( defined $opt_p ) { $parameters{lastzOptions} = $opt_p; }
    else { $parameters{lastzOptions} = ""; }

    if ( defined $opt_l ) { $parameters{logFile} = $opt_l; }
    else { $parameters{logFile} = "$parameters{outputDir}/doLastz.log"; }

    if ( defined $opt_r ) { 
        $parameters{remote} = 1; 
        mkdir ( "$parameters{outputDir}/lastz_shell" ) if ( not -e "$parameters{outputDir}/lastz_shell" );

        $parameters{remote} = 0; 
        print STDERR `echo "not yet implemented. still running in local" > $parameters{outputDir}/lastz_shell/note`;
    }

    return ( %parameters );
}

sub readFastaList {
    my $listFile = shift;
    my %parameters = @_;

    my %fasta = ();
    my @fastaFile = ();
    open ( LIST, $listFile ) or die "Cannot open $listFile for reading!\n";
    while ( my $line = <LIST> ) {
        next if ( $line =~ /^#/ );
        #if ( defined $parameters{skip} ) { next if ( $line =~ /$parameters{skip}/ ); }
        chomp $line;

        my ( $name, $file, @fields ) = split ( /\t/, $line );
        #$name = "chr" . $name if ( $name !~ /chr/ );
        $fasta{$name} = $file;
        push @fastaFile, $file;
    }
    close LIST;

    if ( not defined $fasta{genome} ) {
        my $genomeFile = "$parameters{outputDir}/genome_" . int ( rand ( 9000000 ) + 1000000 ) . ".fa";
        foreach my $scf ( @fastaFile ) { print STDERR `cat $scf >> $genomeFile`; }
        $fasta{genome} = $genomeFile;
    }

    return \%fasta;
}

## obsolete
sub _estimateTotal
{
    my $ref_array1 = shift;
    my $ref_array2 = shift;

    return ( scalar ( keys ( %{$ref_array1} ) ) * scalar ( keys ( %{$ref_array2} ) ) );
}

