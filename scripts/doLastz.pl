#! /usr/bin/perl
#
# run lastz jobs for two lists of fasta files
#

use strict;
use warnings;
use Getopt::Std;

use vars qw ($opt_h $opt_V $opt_D $opt_1 $opt_2 $opt_o $opt_p $opt_r $opt_l $opt_z $opt_s $opt_b );
&getopts('hVD1:2:o:p:l:z:s:b:r');

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
 -b   number of lastz jobs in a batch 

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

    my $count = 0; my $lastPerc = 0; my $countTotalPercent = scalar ( keys %{$ref_fastaList1} );
    my $batch = 0; if ( ( $countTotalPercent > 1000 ) and ( not defined $parameters{batch} ) ) { $parameters{batch} = 1000; }

    my $lastzCmd = "";
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
        if ( not $parameters{remote} ) {
            $lastzCmd = "$parameters{lastzBin} $fasta1 $fasta2 $parameters{lastzOptions} > $outputLav";
            print LOG "\n$lastzCmd\n\tTime: ", `date`;
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
            $lastzCmd .= "$parameters{lastzBin} $fasta1 $fasta2 $parameters{lastzOptions} > $outputLav";
            $lastzCmd .= "\nstatus=\$\?";
            $lastzCmd .= "\nif [ \$status -eq 0 ] \n    then touch $outputLav.done \nfi";
            $lastzCmd .= "\n$parameters{lavToPslBin} $outputLav $outputPsl";
            $lastzCmd .= "\nstatus=\$\?";
            $lastzCmd .= "\nif [ \$status -eq 0 ] \n    then touch $outputPsl.done \nfi\n\n";

            my $jobLabel = $fastaName1;
            if ( defined $parameters{batch} ) {
                if  ( $count % $parameters{batch} != 0 ) { next; }
                else {
                    $batch++;
                    $jobLabel = "jobBatch." . $batch;
                }
            }

            my $lastzJobName = "$jobLabel.lastz";
            my $lastzScriptFile = $parameters{outputDir} . "/lastz_shell/" . $jobLabel . ".sh";
            my $lastzJobSubmitted = $parameters{outputDir} . "/lastz_shell/" . $jobLabel . ".sh.submitted";
            my $lastzJobMemory = "10G";  my $lastzJobTime = "168:00:00";
            remoteCommand( $lastzJobName, $lastzJobMemory, $lastzJobTime, $lastzScriptFile, $lastzCmd, $lastzJobSubmitted );
            $lastzCmd = "";
        }
    }

    if ( ( defined $parameters{batch} ) and ( $count % $parameters{batch} != 0 ) ) { 
        $batch++;
        my $jobLabel = "jobBatch." . $batch;
        my $lastzJobName = "$jobLabel.lastz";
        my $lastzScriptFile = $parameters{outputDir} . "/lastz_shell/" . $jobLabel . ".sh";
        my $lastzJobSubmitted = $parameters{outputDir} . "/lastz_shell/" . $jobLabel . ".sh.submitted";
        my $lastzJobMemory = "10G";  my $lastzJobTime = "168:00:00";
        remoteCommand( $lastzJobName, $lastzJobMemory, $lastzJobTime, $lastzScriptFile, $lastzCmd, $lastzJobSubmitted );
    }

    1;
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
        if ( defined $opt_b ) { $parameters{batch} = $opt_b }

        #    $parameters{remote} = 0; 
        #    print STDERR `echo "not yet implemented. still running in local" > $parameters{outputDir}/lastz_shell/note`;
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

sub remoteCommand {
    my ( $jobName, $jobMemory, $jobTime, $scriptFile, $cmd, $jobSubmittedTag, $largeJob ) = @_;

    open ( SH, ">$scriptFile" ) or die "Cannot open script file $scriptFile for writing.\n";
    print SH "#! /bin/bash\n\n";
    print SH "source ~/.bash_profile\n\n";
    print SH $cmd, "\n\n";

    my $isLargeJob = ""; if ( ( defined $largeJob ) and $largeJob ) { $isLargeJob = " -P large_mem" }

    if ( defined $jobSubmittedTag ) {
        if ( $jobSubmittedTag ) {
            print STDERR `qsub -cwd -N $jobName $isLargeJob -l h_vmem=$jobMemory -l h_rt=$jobTime -m ea -w e $scriptFile`;
            #print STDERR `qsub -cwd -N $jobName -l h_vmem=$jobMemory -l h_rt=$jobTime -m ea -M qczhang\@stanford.edu -w e $scriptFile`;
            ## there must be a simple way to check the submit status and the return of the job id
            print STDERR `touch $jobSubmittedTag`;
        }
        else {
            print STDERR "please submit the job using the following command:\n";
            print STDERR "qsub -cwd -N $jobName $isLargeJob -l h_vmem=$jobMemory -l h_rt=$jobTime -m ea -w e $scriptFile\n";
        }
    }
    else { print STDERR `qsub -cwd -N $jobName $isLargeJob -l h_vmem=$jobMemory -l h_rt=$jobTime -m ea -w e $scriptFile`; }

    1;
}


## obsolete
sub _estimateTotal
{
    my $ref_array1 = shift;
    my $ref_array2 = shift;

    return ( scalar ( keys ( %{$ref_array1} ) ) * scalar ( keys ( %{$ref_array2} ) ) );
}

