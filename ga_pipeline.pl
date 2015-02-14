#! /usr/bin/perl
# wrapper of genome alignment pipeline
# copy right qiangfeng.zhang@gmail.com
# history: 0.01 
#   date: 01/13/2015

use strict;
use warnings;
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
    if ( $_debug ) { print STDERR "\nInitialized..\n"; foreach my $key ( sort { $a cmp $b } ( keys %config ) ) { print STDERR $key, "\t", $config{$key}, "\n"; }  }

    &checkForGenomeFasta ();
    &checkForGenome2Bit ();
    &checkForChrAndSize ();
    if ( $_debug ) { print STDERR "\nFiles ready..\n"; foreach my $key ( sort { $a cmp $b } ( keys %config ) ) { print STDERR $key, "\t", $config{$key}, "\n"; }  }

    print STDERR "\n$config{SCRIPTS}/doLastz.pl -1 $config{list1} -2 $config{list2} -o $config{outDir} -z $config{LASTZ}/lastz -m $config{LASTZMATRIX} -p $config{LASTZMOREOPT} -s $config{UCSC_TOOLS}/lavToPsl\n";
    print STDERR `$config{SCRIPTS}/doLastz.pl -1 $config{list1} -2 $config{list2} -o $config{outDir} -z $config{LASTZ}/lastz -m $config{LASTZMATRIX} -p $config{LASTZMOREOPT} -s $config{UCSC_TOOLS}/lavToPsl\n`;

    if ( not $? ) {
        print STDERR "\n$config{SCRIPTS}/chainning.pl -i $config{genomeSize1} -1 $config{genomeTwoBit1} -2 $config{genomeTwoBit2} -d $config{outDir}/psl -o $config{outDir} -a $config{UCSC_TOOLS}/axtChain -p $config{AXTCHAINOPT} -t $config{UCSC_TOOLS}/chainAntiRepeat\n";
        print STDERR `$config{SCRIPTS}/chainning.pl -i $config{genomeSize1} -1 $config{genomeTwoBit1} -2 $config{genomeTwoBit2} -d $config{outDir}/psl -o $config{outDir} -a $config{UCSC_TOOLS}/axtChain -p $config{AXTCHAINOPT} -t $config{UCSC_TOOLS}/chainAntiRepeat`;
    }
    else { die "Error in finishing lastz jobs\n"; }

    if ( not $? ) {
        print STDERR "\n$config{SCRIPTS}/netting.pl -1 $config{genomeSize1} -2 $config{genomeSize2} -d $config{outDir}/chain -o $config{outDir} -m $config{UCSC_TOOLS}/chainMergeSort -p $config{UCSC_TOOLS}/chainPreNet -n $config{UCSC_TOOLS}/chainNet -s $config{UCSC_TOOLS}/netSyntenic -c $config{UCSC_TOOLS}/netChainSubset -t $config{UCSC_TOOLS}/chainStitchId\n";
        print STDERR `$config{SCRIPTS}/netting.pl -1 $config{genomeSize1} -2 $config{genomeSize2} -d $config{outDir}/chain -o $config{outDir} -m $config{UCSC_TOOLS}/chainMergeSort -p $config{UCSC_TOOLS}/chainPreNet -n $config{UCSC_TOOLS}/chainNet -s $config{UCSC_TOOLS}/netSyntenic -c $config{UCSC_TOOLS}/netChainSubset -t $config{UCSC_TOOLS}/chainStitchId`;
    }
    else { die "Error in chainning\n"; }

    if ( not $? ) { 
        print STDERR "genome alignment successful! Check for over.chain file in $config{outDir}\n\tTime: ", `date`; 
        print STDERR "Now swap the two species in the over.chain file to get reverse.over.chain:\n";
        print STDERR `$config{UCSC_TOOLS}/chainSwap $config{outDir}/over.chain $config{outDir}/reverse.over.chain`;

    }
    else { die "Error in netting\n"; }

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
        if ( not -e $configFile ) { if ( defined $ENV{"GENOMEALIGNMENT"} ) { $configFile = $ENV{"GENOMEALIGNMENT"} . "/.config"; } }
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

    open ( CONFIG, $configFile ) or ( die "Error in openning configuration file $configFile\n\t...cannot configure pipeline to run!\n" );
    print STDERR "Configure pipeline from $configFile\n";
    while ( my $line = <CONFIG> ) {
        next if ( $line =~ /^#/ );
        chomp $line;

        my ( $key, $value ) = ( $line =~ /^(\S+)\s+(.+)$/ );
        $value =~ s/^\s+//; $value =~ s/\s+$//;
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

sub checkForGenome2Bit
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


