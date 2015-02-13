#! /usr/bin/perl
# wrapper of icSHAPE pipeline
# copy right qiangfeng.zhang@gmail.com
# history: 0.01 
#   date: 01/06/2015

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

my $_debug = 1;
my %config = ();

use vars qw ($opt_h $opt_V $opt_D $opt_1 $opt_2 $opt_o $opt_c $opt_a $opt_b $opt_A $opt_B $opt_s $opt_t );
&getopts('hVD1:2:o:c:a:b:A:B:s:t:');

my $usage = <<_EOH_;
## --------------------------------------

Command:
$0 -1 genome1_file -2 genome2_file -o output_directory

# what it is:
 -1     a list of either the whole genome 1 or all of its chromosome/scaffold file locations
 -2     a list of either the whole genome 2 or all of its chromosome/scaffold file locations

 -a     fasta of genome1
 -b     fasta of genome2

 -A     2bit of genome1
 -B     2bit of genome2

 -s     sizes of genome1
 -t     sizes of genome2

 -o     output directory

# more options:
 -c     configuration file
        set parameters for lastz, chaining and others

_EOH_
;

&main ();

sub main
{
    &init ();
    if ( $_debug ) { foreach my $key ( keys %config ) { print $key, "\t", $config{$key}, "\n"; }  }

    &checkForGenomeFasta ();
    &checkForGenome2Bit ();
    &checkForChrAndSize ();
    if ( $_debug ) { foreach my $key ( keys %config ) { print $key, "\t", $config{$key}, "\n"; }  }
    exit;

    #  set parameters and do lastz

    #  if successful, set parameters and continue to chainning
    #  if successful, set parameters and continue to netting
    #
    1;
}

sub init
{
    die $usage if ( $opt_h or 
        ( ( ( not defined $opt_1 ) or ( not defined $opt_2 ) ) and ( ( not defined $opt_a ) or ( not defined $opt_b ) ) and ( ( not defined $opt_A ) or ( not defined $opt_B ) ) ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    ## looking for default config
    my $configFile = "";
    if ( defined $opt_c ) { $configFile = $opt_c; }
    else {
        $configFile = ".config";
        if ( not -e $configFile ) { if ( defined $ENV{"ICSHAPE"} ) { $configFile = $ENV{"ICSHAPE"} . "/.config"; } }
    }

    &config_pipeline ( $configFile );

    $config{list1} = $opt_1 if ( defined $opt_1 );
    $config{list2} = $opt_2 if ( defined $opt_2 );
    $config{genomeFasta1} = $opt_a if ( defined $opt_a );
    $config{genomeFasta2} = $opt_b if ( defined $opt_b );
    $config{genomeTwobit1} = $opt_A if ( defined $opt_A );
    $config{genomeTwobit2} = $opt_B if ( defined $opt_B );
    $config{genomeSize1} = $opt_A if ( defined $opt_s );
    $config{genomeSize2} = $opt_B if ( defined $opt_t );

    if ( defined $opt_o ) { $config{outDir} = $opt_o; }
    else { $config{outDir} = "out$$"; }

    if ( $config{outDir} =~ /^~/ ) { $config{outDir} = $ENV{HOME} . substr ( $config{outDir}, 1 ); }
    elsif ( $config{outDir} !~ /^\// ) {
        my $pwd = `pwd`; chomp $pwd;
        $config{outDir} = $pwd . "/" . $config{outDir};
    }
    if ( not -e $config{outDir} ) {
        print STDERR `mkdir $config{outDir}`;
        die "Cannot output to $config{outDir}!\n" if ( not -e $config{outDir} );
    }
    if ( not -e "$config{outDir}/data" ) {
        print STDERR `mkdir "$config{outDir}/data"`;
        die "Cannot output to $config{outDir}/data!\n" if ( not -e "$config{outDir}/data" );
    }

    1;
}

sub config_pipeline
{
    my $configFile = shift;

    open ( CONFIG, $configFile ) or ( die "Cannot configure pipeline to run!\n" );
    while ( my $line = <CONFIG> ) {
        next if ( $line =~ /^#/ );
        chomp $line;

        my ( $key, $value ) = ( $line =~ /^(\S+)\s+(.+)$/ );
        $value =~ s/^\s+//; $value =~ s/\s+$//;
        $value =~ s/^"//; $value =~ s/"$//;
        $config{$key} = $value;
    }
    close CONFIG;

    1;
}

sub checkForGenomeFasta
{
    if ( defined $config{list1} ) { 
        if ( not defined $config{genomeFasta1} ) { 
            my $tmptFile = $config{outDir} . "/data/genome1.fa";
            $config{genomeFasta1} = mergeFastaList ( $config{list1}, $tmptFile ); 
        } 
    }
    if ( defined $config{list2} ) { 
        if ( not defined $config{genomeFasta2} ) { 
            my $tmptFile = $config{outDir} . "/data/genome2.fa";
            $config{genomeFasta2} = mergeFastaList ( $config{list2}, $tmptFile ); 
        } 
    }

    if ( defined $config{genomeTwobit1} ) { 
        if ( not defined $config{genomeFasta1} ) { 
            $config{genomeFasta1} = $config{outDir} . "/data/genome1.fa";
            print STDERR `$config{UCSC_TOOLS}/twoBitToFa $config{genomeTwobit1} $config{genomeFasta1}`; 
        } 
    }
    if ( defined $config{genomeTwobit2} ) { 
        if ( not defined $config{genomeFasta2} ) { 
            $config{genomeFasta2} = $config{outDir} . "/data/genome2.fa";
            print STDERR `$config{UCSC_TOOLS}/twoBitToFa $config{genomeTwobit2} $config{genomeFasta2}`; 
        } 
    }

    1;
}

sub checkForChrAndSize
{
    my $tmpSize = $config{outDir} . "/data/genome1.tmp.sizes";
    if ( not defined $config{list1} ) { 
        $config{list1} = $config{outDir} . "/data/genome1.lst";
        print STDERR `echo "genome	$config{genomeFasta1}" > $config{list1}`;
        print STDERR `$config{SCRIPTS}/getChr.pl $config{genomeFasta1} $config{outDir}/data/genome1 1> $tmpSize 2>> $config{list1}`; 
    } 
    if ( not defined $config{genomeSize1} ) {
        if ( -e $tmpSize ) { $config{genomeSize1} = $tmpSize; }
        else {
            print STDERR `$config{SCRIPTS}/getChr.pl $config{genomeFasta1} 1> $tmpSize`; 
            if ( -e $tmpSize ) { $config{genomeSize1} = $tmpSize; }
        }
    }

    $tmpSize = $config{outDir} . "/data/genome2.tmp.sizes";
    if ( not defined $config{list2} ) { 
        $config{list2} = $config{outDir} . "/data/genome2.lst";
        print STDERR `echo "genome	$config{genomeFasta2}" > $config{list2}`;
        print STDERR `$config{SCRIPTS}/getChr.pl $config{genomeFasta2} $config{outDir}/data/genome2 1> $tmpSize 2>> $config{list2}`; 
    } 
    if ( not defined $config{genomeSize2} ) {
        if ( -e $tmpSize ) { $config{genomeSize2} = $tmpSize; }
        else {
            print STDERR `$config{SCRIPTS}/getChr.pl $config{genomeFasta2} 1> $tmpSize`; 
            if ( -e $tmpSize ) { $config{genomeSize2} = $tmpSize; }
        }
    }

    1;
}

sub checkForGenome2bit
{
    if ( not defined $config{genomeTwoBit1} ) { 
        $config{genomeTwoBit1} = $config{outDir} . "/data/genome1.2bit";
        print STDERR `$config{UCSC_TOOLS}/faToTwoBit $config{genomeFasta1} $config{genomeTwoBit1}`;
    } 
    if ( not defined $config{genomeTwoBit2} ) { 
        $config{genomeTwoBit2} = $config{outDir} . "/data/genome2.2bit";
        print STDERR `$config{UCSC_TOOLS}/faToTwoBit $config{genomeFasta2} $config{genomeTwoBit2}`;
    } 

    1;
}

sub mergeFastaList 
{
    my $listFile = shift;
    my $merged = shift;

    my @fastaFile = ();
    my $foundGenome = 0;
    open ( LIST, $listFile ) or die "Cannot open $listFile for reading!\n";
    my $genomeFile = "";
    while ( my $line = <LIST> ) {
        next if ( $line =~ /^#/ );
        chomp $line;

        my ( $name, $file, @fields ) = split ( /\t/, $line );
        if ( $name eq "genome" ) { $foundGenome = 1; $genomeFile = $file;  last; }
        push @fastaFile, $file;
    }
    close LIST;

    if ( not $foundGenome ) {
        $genomeFile = $merged;
        foreach my $scf ( @fastaFile ) { print STDERR `cat $scf >> $genomeFile`; }
        open ( LS, ">>$listFile" ) or die "Cannot add information to $listFile\n";
        print LS "genome\t", $genomeFile, "\n";
        close LS;
    }

    return $genomeFile;
}


