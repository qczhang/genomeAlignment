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

use vars qw ($opt_h $opt_V $opt_D $opt_1 $opt_2 $opt_o $opt_c $opt_a $opt_b $opt_A $opt_B $opt_s $opt_t $opt_j );
&getopts('hVD1:2:o:c:a:b:A:B:s:t:j:');

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
        
 -j     job [all(default)|prepare|lastz|chainning|netting]

_EOH_
;

&main ();

sub main
{
    &init ();
    if ( $_debug ) { print STDERR "\nInitialized..\n";  print STDERR &printConfigure(), "\n"; }

    &checkForGenomeFasta ( job => $config{job} );
    &checkForGenome2Bit ( job => $config{job} );
    &checkForChrAndSize ( job => $config{job} );
    if ( $_debug ) { print STDERR "Files are ready.\n\n"; print STDERR &printConfigure(), "\n"; }

    if ( defined $config{job}{lastz} ) {
        print STDERR "\n$config{SCRIPTS}/doLastz.pl -1 $config{list1} -2 $config{list2} -o $config{outDir} -z $config{LASTZ}/lastz -p $config{LASTZOPT} -s $config{UCSC_TOOLS}/lavToPsl\n";
        print STDERR `$config{SCRIPTS}/doLastz.pl -1 $config{list1} -2 $config{list2} -o $config{outDir} -z $config{LASTZ}/lastz -p $config{LASTZOPT} -s $config{UCSC_TOOLS}/lavToPsl`;
    }

    if ( not $? ) {
        if ( defined $config{job}{chainning} ) {
            print STDERR "\n$config{SCRIPTS}/chainning.pl -i $config{genomeSize1} -1 $config{genomeTwoBit1} -2 $config{genomeTwoBit2} -d $config{outDir}/psl -o $config{outDir} -a $config{UCSC_TOOLS}/axtChain -p $config{AXTCHAINOPT} -t $config{UCSC_TOOLS}/chainAntiRepeat\n";
            print STDERR `$config{SCRIPTS}/chainning.pl -i $config{genomeSize1} -1 $config{genomeTwoBit1} -2 $config{genomeTwoBit2} -d $config{outDir}/psl -o $config{outDir} -a $config{UCSC_TOOLS}/axtChain -p $config{AXTCHAINOPT} -t $config{UCSC_TOOLS}/chainAntiRepeat`;
        }
    }
    else { die "Error in finishing lastz jobs\n"; }

    if ( not $? ) {
        if ( defined $config{job}{netting} ) {
            print STDERR "\n$config{SCRIPTS}/netting.pl -1 $config{genomeSize1} -2 $config{genomeSize2} -d $config{outDir}/chain -o $config{outDir} -m $config{UCSC_TOOLS}/chainMergeSort -p $config{UCSC_TOOLS}/chainPreNet -n $config{UCSC_TOOLS}/chainNet -s $config{UCSC_TOOLS}/netSyntenic -c $config{UCSC_TOOLS}/netChainSubset -t $config{UCSC_TOOLS}/chainStitchId\n";
            print STDERR `$config{SCRIPTS}/netting.pl -1 $config{genomeSize1} -2 $config{genomeSize2} -d $config{outDir}/chain -o $config{outDir} -m $config{UCSC_TOOLS}/chainMergeSort -p $config{UCSC_TOOLS}/chainPreNet -n $config{UCSC_TOOLS}/chainNet -s $config{UCSC_TOOLS}/netSyntenic -c $config{UCSC_TOOLS}/netChainSubset -t $config{UCSC_TOOLS}/chainStitchId`;
        }
    }
    else { die "Error in chainning\n"; }

    if ( not $? ) { 
        if ( defined $config{job}{swaping} ) {
            print STDERR "netting successful! Check for over.chain file in $config{outDir}\n\tTime: ", `date`; 
            print STDERR "Now swap the two species in the over.chain file to get reverse.over.chain:\n";
            print STDERR `$config{UCSC_TOOLS}/chainSwap $config{outDir}/over.chain $config{outDir}/reverse.over.chain`;
            print STDERR "and reverse.over.chain file in $config{outDir} for reverse chain\n\tTime: ", `date` if ( not $? ); 
        }
    }
    else { die "Error in netting\n"; }

    if ( not $? ) { 
        print STDERR "genome alignment successful! Check for over.chain file in $config{outDir}\n\tTime: ", `date`; 
        print STDERR "and reverse.over.chain file in $config{outDir} for reverse chain\n\tTime: ", `date` if ( not $? ); 
    }

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

    print STDERR "Initializing...\n\t", `date`;
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

    my %allJobs = ( "prepare" => 1, "lastz" => 1, "chainning" => 1, "netting" => 1, "swaping" => 1 );
    if ( defined $opt_j ) {
        my @jobs = split ( /,/, $opt_j ); 
        foreach my $job ( @jobs ) {
            if ( defined $allJobs{$job} ) { $config{job}{$job} = 1; }
            elsif ( $job eq "all" ) {
                foreach my $aj ( keys %allJobs ) {
                    $config{job}{$aj} = 1;
                    last;
                }
            }
        }
    }
    else { foreach my $aj ( keys %allJobs ) { $config{job}{$aj} = 1; } }
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
    my %parameters = @_;

    print STDERR "Check for genome fasta...\n\t", `date`;
    if ( defined $config{list1} ) { 
        if ( not defined $config{genomeFasta1} ) { 
            my ( $foundGenome, @fastaFile ) = searchInList ( $config{list1} );
            if ( $foundGenome ) { $config{genomeFasta1} = $foundGenome; }
            else {
                my $tmptListFile = $config{outDir} . "/data/genome1.lst";
                my $genomeFile = $config{outDir} . "/data/genome1.fa";
                if ( defined $parameters{job}->{prepare} ) {
                    $config{list1} = updateList ( $config{list1}, $tmptListFile, $genomeFile );
                    $config{genomeFasta1} = mergeFastaList ( $genomeFile, @fastaFile ); 
                }
                else {
                    $config{list1} = $tmptListFile;
                    $config{genomeFasta1} = $genomeFile;
                }
            }
        } 
    }
    if ( defined $config{list2} ) { 
        if ( not defined $config{genomeFasta2} ) { 
            my ( $foundGenome, @fastaFile ) = searchInList ( $config{list2} );
            if ( $foundGenome ) { $config{genomeFasta2} = $foundGenome; }
            else {
                my $tmptListFile = $config{outDir} . "/data/genome2.lst";
                my $genomeFile = $config{outDir} . "/data/genome2.fa";
                if ( defined $parameters{job}->{prepare} ) {
                    $config{list2} = updateList ( $config{list2}, $tmptListFile, $genomeFile );
                    $config{genomeFasta2} = mergeFastaList ( $genomeFile, @fastaFile ); 
                }
                else {
                    $config{list2} = $tmptListFile;
                    $config{genomeFasta2} = $genomeFile;
                }
            }
        } 
    }

    if ( defined $config{genomeTwobit1} ) { 
        if ( not defined $config{genomeFasta1} ) { 
            $config{genomeFasta1} = $config{outDir} . "/data/genome1.fa";
            if ( defined $parameters{job}->{prepare} ) { print STDERR `$config{UCSC_TOOLS}/twoBitToFa $config{genomeTwobit1} $config{genomeFasta1}`; }
        } 
    }
    if ( defined $config{genomeTwobit2} ) { 
        if ( not defined $config{genomeFasta2} ) { 
            $config{genomeFasta2} = $config{outDir} . "/data/genome2.fa";
            if ( defined $parameters{job}->{prepare} ) { print STDERR `$config{UCSC_TOOLS}/twoBitToFa $config{genomeTwobit2} $config{genomeFasta2}`; }
        } 
    }

    1;
}

sub checkForGenome2Bit
{
    my %parameters = @_;

    print STDERR "Check for genome 2bit...\n\t", `date`;
    if ( not defined $config{genomeTwoBit1} ) { 
        $config{genomeTwoBit1} = $config{outDir} . "/data/genome1.2bit";
        if ( defined $parameters{job}->{prepare} ) {  print STDERR `$config{UCSC_TOOLS}/faToTwoBit $config{genomeFasta1} $config{genomeTwoBit1}`; }
    } 
    if ( not defined $config{genomeTwoBit2} ) { 
        $config{genomeTwoBit2} = $config{outDir} . "/data/genome2.2bit";
        if ( defined $parameters{job}->{prepare} ) {  print STDERR `$config{UCSC_TOOLS}/faToTwoBit $config{genomeFasta2} $config{genomeTwoBit2}`; }
    } 

    1;
}

sub checkForChrAndSize
{
    my %parameters = @_;

    print STDERR "Check for chromosome and sizes...\n\t", `date`;
    my $tmpSize = $config{outDir} . "/data/genome1.tmp.sizes";
    if ( not defined $config{list1} ) { 
        $config{list1} = $config{outDir} . "/data/genome1.lst";
        if ( defined $parameters{job}->{prepare} ) {
            print STDERR `echo "genome	$config{genomeFasta1}" > $config{list1}`;
            print STDERR `$config{SCRIPTS}/getChr.pl $config{genomeFasta1} $config{outDir}/data/genome1 1> $tmpSize 2>> $config{list1}`; 
        }
    } 
    if ( not defined $config{genomeSize1} ) {
        if ( -e $tmpSize ) { $config{genomeSize1} = $tmpSize; }
        else {
            if ( defined $parameters{job}->{prepare} ) {
                print STDERR `$config{SCRIPTS}/getChr.pl $config{genomeFasta1} 1> $tmpSize`; 
                if ( -e $tmpSize ) { $config{genomeSize1} = $tmpSize; }
            }
        }
    }

    $tmpSize = $config{outDir} . "/data/genome2.tmp.sizes";
    if ( not defined $config{list2} ) { 
        $config{list2} = $config{outDir} . "/data/genome2.lst";
        if ( defined $parameters{job}->{prepare} ) {
            print STDERR `echo "genome	$config{genomeFasta2}" > $config{list2}`;
            print STDERR `$config{SCRIPTS}/getChr.pl $config{genomeFasta2} $config{outDir}/data/genome2 1> $tmpSize 2>> $config{list2}`; 
        }
    } 
    if ( not defined $config{genomeSize2} ) {
        if ( -e $tmpSize ) { $config{genomeSize2} = $tmpSize; }
        else {
            if ( defined $parameters{job}->{prepare} ) {
                print STDERR `$config{SCRIPTS}/getChr.pl $config{genomeFasta2} 1> $tmpSize`; 
                if ( -e $tmpSize ) { $config{genomeSize2} = $tmpSize; }
            }
        }
    }

    1;
}

sub searchInList 
{
    my $list = shift;

    print STDERR "Search for genome fasta file in input list\n\t", `date`;
    open ( LIST, $list ) or die "Cannot open $list for reading!\n";
    my $foundGenome = 0;
    my @fastaFile = ();
    while ( my $line = <LIST> ) {
        next if ( $line =~ /^#/ ); chomp $line;
        my ( $name, $file, @fields ) = split ( /\t/, $line );
        if ( $name eq "genome" ) { $foundGenome = $file;  last; }
        push @fastaFile, $file;
    }
    close LIST;

    return ( $foundGenome, @fastaFile );
}

sub mergeFastaList 
{
    my $merged = shift;
    my @toMerge = @_;

    print STDERR "Merge chromosome fasta to get genome fasta file\n\t", `date`;
    foreach my $file ( @toMerge ) { print STDERR `cat $file >> $merged`; }
    return $merged;
}

sub updateList
{
    my $oldListFile = shift;
    my $newListFile = shift;
    my $genomeFile = shift;

    print STDERR "Update fasta list\n\t", `date`;
    print STDERR `/bin/cp $oldListFile $newListFile`;
    open ( LS, ">>$newListFile" ) or die "Cannot add information to $newListFile\n";
    print LS "genome\t", $genomeFile, "\n";
    close LS;

    return $newListFile;
}

sub printConfigure
{
    my $conf = "";

    foreach my $key ( sort { $a cmp $b } ( keys %config ) ) { 
        if ( $key ne "job" ) { $conf .= $key . "\t" . $config{$key} . "\n"; }
        else {
            my $jobs = "";
            foreach my $j ( "prepare", "lastz", "chainning", "netting", "swaping" ) { if ( defined $config{job}{$j} ) { $jobs .= "," . $j; } }
            $jobs =~ s/^,//;
            $conf .= $key . "\t" . $jobs . "\n";
        }
    }  

    return $conf;
}

