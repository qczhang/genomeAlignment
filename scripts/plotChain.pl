#! /usr/bin/perl
#
# generate gnuplot scripts from chain files
#
use strict;
use warnings;

use Getopt::Std;

my @preDefinedChromosome = ( "chr1", "chr2", "chr2a", "chr2b", "chr2L", "chr2R", "chr3", "chr3L", "chr3R", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chrX", "chrY", "chrM" );

use vars qw ($opt_h $opt_V $opt_D $opt_c $opt_o $opt_a $opt_b $opt_s $opt_x $opt_y );
&getopts('hVDc:o:a:b:s:x:y:');

my $usage = <<_EOH_;

Use gnuplot to draw dotplot with chain file as input. 

Command:
$0 -c chain file <options>

# options:
 -o   output file (default=name.png)
 -a   plot chromosome ID 1
 -b   plot chromosome ID 2
 -s   score cutoff
 -x   chromosome size file 1
 -y   chromosome size file 2

# more options:
 -h   print help
 -V   Verbose
 -D   Debug

_EOH_
    ;

&main();

sub main 
{
    my ( $inFile, $outFile, $plot1, $plot2, $score_cutoff, $chromosomeSizeFile1, $chromosomeSizeFile2 ) = &init();
    print "Parameters are:\n\t", join "\n\t", $inFile, $outFile, $plot1, $plot2, $score_cutoff, $chromosomeSizeFile1, $chromosomeSizeFile2, "\n" if ( $opt_D );   

    my %tchr_size = (); %tchr_size = &getChromosomeSize ( $chromosomeSizeFile1 ) if ( $chromosomeSizeFile1 ne "NULL" );
    my %qchr_size = (); %qchr_size = &getChromosomeSize ( $chromosomeSizeFile2 ) if ( $chromosomeSizeFile2 ne "NULL" );

    my ( $ref_plot_coordinates, $ref_tchr_inData, $ref_qchr_inData ) 
        = &readChainFile ( $inFile, $plot1, $plot2, $score_cutoff, $chromosomeSizeFile1, $chromosomeSizeFile2, \%tchr_size, \%qchr_size );

    foreach my $chr ( keys %$ref_tchr_inData ) {   $$ref_tchr_inData{$chr} = $tchr_size{$chr};  }
    foreach my $chr ( keys %$ref_qchr_inData ) {   $$ref_qchr_inData{$chr} = $qchr_size{$chr};  }

    my @tchrlist = &defineChromosomes ( $ref_tchr_inData, 
        lengthcutoff => 10000,
        skip => "hap:random"
    );
    my @qchrlist = &defineChromosomes ( $ref_qchr_inData, 
        lengthcutoff => 10000,
        skip => "hap:random"
    );

    my $tchrTotalLen = 0; my $qchrTotalLen = 0;
    foreach my $tchr ( @tchrlist ) {   $tchrTotalLen += $tchr_size{$tchr};   }
    foreach my $qchr ( @qchrlist ) {   $qchrTotalLen += $qchr_size{$qchr};   }

    my  ( $plotChrGridFile, $plotLabelFile, $plotSynFile ) 
        = &genPlotFile ( $ref_plot_coordinates, \@tchrlist, \@qchrlist, \%tchr_size, \%qchr_size, $tchrTotalLen, $qchrTotalLen );
    my $gnuFile = &genGNUfile ( $outFile, $plotChrGridFile, $plotLabelFile, $plotSynFile, $tchrTotalLen, $qchrTotalLen );
    print STDERR `gnuplot $gnuFile`;
}



sub init 
{
    die $usage if $opt_h;

    my $inFile=$opt_c or die $usage;
    my $outFile;
    if ( $inFile=~/^(.+)\.chain$/ )   {
        if ( defined $opt_a && defined $opt_b ) {  $outFile=$1 ."." . $opt_a . "vs" . $opt_b . ".png";  }
        else {  $outFile=$1 .".png";  }
    }
    $outFile=$opt_o if $opt_o;

    my $plot1 = ( defined $opt_a ) ? $opt_a : "NULL";;
    my $plot2 = ( defined $opt_b ) ? $opt_b : "NULL";;
    my $score_cutoff = ( defined $opt_s ) ? $opt_s : 0;
    my $chromosomeSizeFile1 = ( defined $opt_x ) ? $opt_x : "NULL";
    my $chromosomeSizeFile2 = ( defined $opt_y ) ? $opt_y : "NULL";

    return ( $inFile, $outFile, $plot1, $plot2, $score_cutoff, $chromosomeSizeFile1, $chromosomeSizeFile2 );
}


sub getChromosomeSize 
{
    my $chromosomeSizeFile = shift;

    print STDERR "Read chromsome size file: $chromosomeSizeFile\n";
    if ( not -e $chromosomeSizeFile )
        {  die "$chromosomeSizeFile not exists";  }

    my @chr_size;
    open ( CS, $chromosomeSizeFile );
    while ( my $line = <CS> )  {
        next if ( $line =~ /^#/ );
        chomp $line;
        next if ( not $line );

        $line =~ s/\s/\t/g;
        $line =~ s/\t+/\t/g;

        my ( $chr, $size ) = split ( /\t/, $line );
        push @chr_size, $chr, $size;
    }

    print STDERR scalar(@chr_size), " chromosomes read.\n";
    return @chr_size;
}


sub updateChromosomeSize
{
    my ( $chr, $end, $ref_chr_size ) = @_;

    if ( ( not defined $$ref_chr_size{$chr} ) || ( $$ref_chr_size{$chr} < $end ) )
        {   $$ref_chr_size{$chr} = $end;  }

    return 1;
}


sub defineChromosomes 
{
    my $ref_chrmap = shift;
    my %parameters = @_;

    print STDERR "Sort chromsomese\n";
    foreach ( keys %$ref_chrmap ) {
        if ( defined $parameters{lengthcutoff} ) {
            delete $$ref_chrmap{$_} if ( $$ref_chrmap{$_} < $parameters{lengthcutoff} ); 
        }
        if ( defined $parameters{skip} ) {
            my @skipWords = split ( /:/, $parameters{skip} );
            foreach my $word ( @skipWords ) {
                if ( $_ =~ /$word/ ) {
                    delete $$ref_chrmap{$_}; 
                }
            }
        }
    }
        
    my @sortedChrList = ();
    if ( defined $parameters{predefined} ) {
        foreach ( @preDefinedChromosome ) {   push @sortedChrList, $_ if ( defined $$ref_chrmap{$_} );   }
    }
    else {
        foreach ( @preDefinedChromosome ) { 
            if ( defined $$ref_chrmap{$_} ) {   push @sortedChrList, $_;  delete $$ref_chrmap{$_};   }
        }
        foreach ( sort {$a cmp $b} ( keys %$ref_chrmap ) ) {   push @sortedChrList, $_;   }
    }
    print STDERR "Chromosomes sorted in order:\n";
    print STDERR "\t@sortedChrList\n";

    return @sortedChrList;
}

sub readChainFile 
{
    my ( $inFile, $plot1, $plot2, $score_cutoff, $chromosomeSizeFile1, $chromosomeSizeFile2, $ref_tchr_size, $ref_qchr_size ) = @_;

    print STDERR "Read chain file: $inFile\n";
    my %plot_coordinates = ();
    my %tchr_inData = (); my %qchr_inData = ();

    my %plot1_name = ();
    foreach ( split ( /:/, $plot1 ) ) { $plot1_name{$_} = 1; }
    my %plot2_name = ();
    foreach ( split ( /:/, $plot2 ) ) { $plot2_name{$_} = 2; }

    my ( $isChain,$score,$tName,$tSize,$tStrand,$tStart,$tEnd,$qName,$qSize,$qStrand,$qStart,$qEnd,$id );
    my $lineCount = 0; my $chainCount = 0;
    my $selected = 0; my $color = 0;
    open ( IN, $inFile ) or die "Can't open $inFile to read.";
    while ( <IN> ) {
        chomp;
        next if ( $_ =~ /^#/ ); next if ( not $_ );

        $lineCount++;
        print STDERR "\t$lineCount lines read\n" if ( $lineCount % 1000000 == 0 );

        if ( $_ =~ /^chain/ ) {
            ( $isChain,$score,$tName,$tSize,$tStrand,$tStart,$tEnd,$qName,$qSize,$qStrand,$qStart,$qEnd,$id ) = split ( /\s/, $_ );
            next if ( $score < $score_cutoff );

            $chainCount++;
            print STDERR "\t$chainCount chains read\n" if ( $chainCount % 10000 == 0 );
            #if ( ( $plot1 eq "NULL" || $plot1 eq $tName ) and ( $plot2 eq "NULL" || $plot2 eq $qName ) ) {
            if ( ( $plot1 eq "NULL" || ( defined $plot1_name{$tName} ) ) and ( $plot2 eq "NULL" || ( defined $plot2_name{$qName} ) ) ) {
                &updateChromosomeSize ( $tName, $tEnd, $ref_tchr_size ) if ( $chromosomeSizeFile1 eq "NULL" );
                &updateChromosomeSize ( $qName, $qEnd, $ref_qchr_size ) if ( $chromosomeSizeFile2 eq "NULL" );

                $selected = 1;
                $color = int ( rand ( 1000000 ) ); 
            }
            else {   $selected = 0;  }

            if ( $selected ) {
                $tchr_inData{$tName} = 1;
                $qchr_inData{$qName} = 1;

                if ( not defined $plot_coordinates{$tName}{$qName} ) {
                    $plot_coordinates{$tName}{$qName} = "";  
                    $color = 1;  
                }
            }
        }
        elsif ( $selected ) {
            my ( $aLen, $dt, $dq ) = split ( /\t/, $_ );
            next if ( not $aLen );
            $tEnd = $tStart + $aLen;
            $qEnd = $qStart + $aLen;

            if ( $tStrand eq $qStrand ) {
                $plot_coordinates{$tName}{$qName} .= "$tStart $qStart $color\n";
                $plot_coordinates{$tName}{$qName} .= "$tEnd $qEnd $color\n";
            } 
            else {
                $plot_coordinates{$tName}{$qName} .= "$tStart $qEnd $color\n";
                $plot_coordinates{$tName}{$qName} .= "$tEnd $qStart $color\n";
            }

            if ( ( defined $dt ) && ( defined $dq ) ) {
                $tStart = $tEnd + $dt;
                $qStart = $qEnd + $dq;
            }
        }
    }
    close IN;

    print STDERR "In total $chainCount chains and $lineCount lines read\n";
    return ( \%plot_coordinates, \%tchr_inData, \%qchr_inData );
}

sub genPlotFile 
{
    ## write plot date file
    my ( $ref_plot_coordinates, $ref_tchrlist, $ref_qchrlist, $ref_tchr_size, $ref_qchr_size, $tchrTotalLen, $qchrTotalLen ) = @_;

    my $plotChrGridFile="plot.grid.$$"; my $plotChrLabelFile="plot.label.$$";
    print STDERR "Generate plot date file: $plotChrGridFile, and $plotChrLabelFile\n";
    if ( $opt_D ) {
        foreach my $key ( keys %$ref_tchr_size ) {    print STDERR "\t$key, $$ref_tchr_size{$key}\n";   }
        foreach my $key ( keys %$ref_qchr_size ) {    print STDERR "\t$key, $$ref_qchr_size{$key}\n";   }
    }

    open ( PLOTG,">$plotChrGridFile" ) or die "Can't open $plotChrGridFile to write.";
    open ( PLOTL,">$plotChrLabelFile" ) or die "Can't open $plotChrLabelFile to write.";
    my $blockStart = 0;
    foreach my $tchr ( @$ref_tchrlist ) {
        print $tchr, "\n" if ( $opt_D );
        if ( $blockStart && ( $blockStart < $tchrTotalLen ) ) {
            print PLOTG "$blockStart 0 10000\n" if ( $blockStart );
            print PLOTG "$blockStart $qchrTotalLen 10000\n\n" if ( $blockStart );
        }

        my $blockCenter = $blockStart + $$ref_tchr_size{$tchr}/2;
        print PLOTL "$tchr $blockCenter 10\n";
        $blockStart += $$ref_tchr_size{$tchr};
    }

    $blockStart = 0;
    foreach my $qchr ( @$ref_qchrlist ) {
        print $qchr, "\n" if ( $opt_D );
        if ( $blockStart && ( $blockStart < $qchrTotalLen ) ) {
            print PLOTG "0 $blockStart 10000\n" if ( $blockStart );
            print PLOTG "$tchrTotalLen $blockStart 10000\n\n" if ( $blockStart );
        }

        my $blockCenter = $blockStart + $$ref_qchr_size{$qchr}/2;
        print PLOTL "$qchr 10 $blockCenter\n";

        $blockStart += $$ref_qchr_size{$qchr};
    }
    close PLOTG;
    close PLOTL;

    my $plotSynFile="plot.syn.$$";
    open ( PLOTS,">$plotSynFile" ) or die "Can't open $plotSynFile to write.";
    my $plotLineCount = 0;
    my $blockStart1 = 0; my $blockStart2 = 0;
    foreach my $tchr ( @$ref_tchrlist ) {
        if ( defined $$ref_plot_coordinates{$tchr} ) {
            print $tchr, "\n" if ( $opt_V );

            $blockStart2 = 0;
            foreach my $qchr ( @$ref_qchrlist ) {
                if ( defined $$ref_plot_coordinates{$tchr}{$qchr} ) {
                    print $qchr, "\t" if ( $opt_V );

                    my @lines = split ( /\n/, $$ref_plot_coordinates{$tchr}{$qchr} );
                    foreach my $line ( @lines ) {
                        my ( $pos1, $pos2, $color ) = split ( /\s/, $line );
                        print PLOTS $pos1+$blockStart1, " ", $pos2+$blockStart2, " ", $color, "\n";

                        $plotLineCount++;
                        print PLOTS "\n" if ( $plotLineCount % 2 == 0 );
                    }
                }

                $blockStart2 += $$ref_qchr_size{$qchr};
            }
        }

        $blockStart1 += $$ref_tchr_size{$tchr};
        print "\n" if ( $opt_V );
    }
    close PLOTS;

    return ( $plotChrGridFile, $plotChrLabelFile, $plotSynFile );
}


sub genGNUfile 
{
    my ( $outFile, $plotChrGridFile, $plotChrLabelFile, $plotSynFile, $max1, $max2 ) = @_;

    my $gnuFile="gnu.$$";
    print STDERR "Generate gnuplot script file: $gnuFile\n";
    my $width = int(sqrt($max1/100000000)*100)*30;
    my $height = int(sqrt($max2/100000000)*100)*30;
    $max1 = sprintf("%12.10e", $max1);
    $max2 = sprintf("%12.10e", $max2);

    open(GNU,">$gnuFile") or die "Can't open $gnuFile to write.";
    print GNU "set terminal png size $width,$height\n";
    print GNU "set output \"$outFile\"\n";
    print GNU "set xrange [0:$max1]\n";
    print GNU "set yrange [0:$max2]\n";
    print GNU "unset colorbox\n";
    print GNU "set palette model HSV\n";
    print GNU "set palette rgb 3,2,2\n";
    print GNU "set style line 1 lt 3 lc rgb \"grey\" lw 1\n";
    print GNU "set style line 2 lt 1 palette lw 2\n";
    print GNU "plot \"$plotChrGridFile\" u 1:2 with lines ls 1 notitle, \\\n";
    print GNU "\t\"$plotChrLabelFile\" u 2:3:1 with labels point offset character 0,character 1 tc rgb \"blue\", \\\n";
    print GNU "\t\"$plotSynFile\" u 1:2:3 with lines ls 2 notitle\n";
    close GNU;

    return $gnuFile;
}

