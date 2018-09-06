#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Which;


################################################################################
# Author : Isabelle GUIGON
# Date : 2014-11-27
# This script converts a BAM file into a BED file for use by miRkwood.
#
# Dependencies : samtools, perl module File::Which
# Make sure to have it installed in your PATH.
# For Ubuntu/Debian distributions `sudo apt-get install samtools` is enough.
################################################################################


########## Variables
my $input_file = '';
my $bed_file   = '';
my $min_length = 18;
my $max_length = 25;
my $time = time();
my $sorted_bam_file = "/tmp/mirkwood_bam2bed_${time}_sorted";
my $sorted_sam_file = "/tmp/mirkwood_bam2bed_${time}_sorted.sam";

my $counts;
my $help;
my $help_message = <<'EOF';
mirkwood-bam2bed.pl
----------
Script to convert a BAM/SAM file into a BED file for use by miRkwood.

Usage : ./mirkwood-bam2bed.pl -in <input BAM/SAM file> -bed <output BED file>

Options :
    -in   : input BAM or SAM
    -bed  : output BED
    -min  : keep only reads with length >= min (default 18)
    -max  : keep only reads with length <= max (default 25)
    -help : display this message and quit

Dependencies : samtools
Make sure to have it installed in your PATH. For Ubuntu/Debian distributions `sudo apt-get install samtools` is enough.

EOF


########## Check that the samtools are installed
if ( ! which( 'samtools' ) ){
    die "ERROR: samtools are missing. Please install them in your PATH.\n";
}


########## Get options
GetOptions ('in=s'  => \$input_file,
            'bed=s' => \$bed_file,
            'min=s' => \$min_length,
            'max=s' => \$max_length,
            'help'  => \$help);


########## Validate options
if ( $help ){
    print $help_message;
    exit;
}

if ( ( ! -r $input_file ) || ( $bed_file eq '' ) ){
    die "Missing parameter!\n" . $help_message;
}


if ( $min_length >= $max_length ){
    die "--min should be strictly lower than --max\n" . $help_message;
}


########## Create sorted SAM file
if ( $input_file =~ /([^\/\\]+)[.]sam/ ){
    system("sort -k 3,3 -k 4,4n $input_file | grep -v \"^@\" > $sorted_sam_file");
}
elsif ( $input_file =~ /([^\/\\]+)[.](bam|dat)/ ){
    ##### Sort BAM file
    system("samtools sort $input_file $sorted_bam_file");

    ##### Convert sorted BAM into SAM and filter out unmapped reads
    system("samtools view -F 4 $sorted_bam_file.bam > $sorted_sam_file");
    unlink $sorted_bam_file . '.bam';
}
else{
    die "Non correct input file. We accept BAM and SAM as input formats.\n$help_message";
}



########## Read the SAM file to store counts into a hash
open(my $SAM, '<', $sorted_sam_file) or die "ERROR while reading SAM file. Program will end prematurely.\n";

while ( <$SAM> ){

    chomp;

    my @line = split ( /\t/smx );

    if ( $line[1] eq '0x4' || $line[1] eq '4' ){
        next;
    }
    if ( length($line[9]) < $min_length || length($line[9]) > $max_length ){
        next;
    }

    my $chromosome = $line[2];
    my $start = $line[3] - 1;
    my $sequence = $line[9];
    my $strand = '+';
    if ( $line[1] eq '16' or $line[1] eq '0x10' ){
        $strand = '-';
    }

    if ( ! exists( $counts->{$chromosome}{$start}{$sequence}{$strand} ) ){
        $counts->{$chromosome}{$start}{$sequence}{$strand} = 0;
    }
    $counts->{$chromosome}{$start}{$sequence}{$strand}++;

}

close $SAM;
unlink $sorted_sam_file;


########## Browse hash tables and print data in BED file
open(my $BED, '>', $bed_file) or die "ERROR while creating $bed_file. Program will end prematurely.\n";
foreach my $chromosome ( sort ( keys%{$counts} ) ){
    foreach my $start ( sort {$a <=> $b} keys%{ $counts->{$chromosome} } ){
        foreach my $sequence ( sort (keys%{ $counts->{$chromosome}{$start} } ) ){
            foreach my $strand ( sort (keys%{ $counts->{$chromosome}{$start}{$sequence} } ) ){
                my $end = $start + length($sequence);
                print $BED "$chromosome\t";
                print $BED "$start\t";
                print $BED "$end\t";
                print $BED "$sequence\t";
                print $BED "$counts->{$chromosome}{$start}{$sequence}{$strand}\t";
                print $BED "$strand\n";
            }
        }
    }
}
close $BED;

my $total_time = time() - $time;
my $day  = int( $total_time / 86_400 );
my $hour = int( ($total_time % 86_400 ) / 3_600 );
my $min  = int( ( ($total_time % 86_400 ) % 3_600 ) / 60 );
my $sec  = int( ( ($total_time % 86_400 ) % 3_600 ) % 60 );

print "Done in $day day $hour h $min min $sec sec.\n";


