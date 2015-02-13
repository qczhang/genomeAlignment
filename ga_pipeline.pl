#! /usr/bin/perl
# wrapper of icSHAPE pipeline
# copy right qiangfeng.zhang@gmail.com
# history: 0.01 
#   date: 01/06/2015

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

my $_debug = 0;
my %config = ();

use vars qw ($opt_h $opt_V $opt_D $opt_1 $opt_2 $opt_o $opt_c );
&getopts('hVD1:2:o:c:');

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

    ## check the input format, check whether genome 2bit/fa files are available, if not generate them
    #  check whether chromsome files are available, if not generate them
    #  generate sizes files

    #  set parameters and do lastz

    #  if successful, set parameters and continue to chainning
    #  if successful, set parameters and continue to netting
    #
    1;
}

sub init
{
    die $usage if ( $opt_h || ( ( ( not defined $opt_1 ) || ( not defined $opt_2 ) ) and ( ( not defined $opt_a ) || ( not defined $opt_b ) ) and ( ( not defined $opt_A ) || ( not defined $opt_B ) ) ) );
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
    $config{genome2bit} = $opt_A if ( defined $opt_A );
    $config{genome2bit} = $opt_B if ( defined $opt_B );

    if ( defined $opt_o ) { $config{outDir} = $opt_o; }
    else { $config{outDir} = "out$$"; }
    if ( not -e $config{outDir} ) {
        print STDERR `mkdir $config{outDir}`;
        die "Cannot output to $config{outDir}!\n" if ( not -e $config{outDir} );
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
