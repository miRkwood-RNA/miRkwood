#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


################################################################################
# Author : Isabelle GUIGON
# Date : 2014-11-27
# This script converts a BAM file into a BED file for use by miRkwood.
#
# Dependancies : samtools
# Make sure to have it installed in your PATH.
# For Ubuntu/Debian distributions `sudo apt-get install samtools` is enough.
################################################################################


########## Variables
my $bam_file = '';
my $sam_file;
my $bed_file;
my $output_directory;
my $tmp_file;
my $counts;
my $names;
my $sequence;
my $id;
my $chromosome;
my $start;
my $end;
my $help;
my $help_message = "mirkwood-bam2bed.pl
----------
Script to convert a BAM into a BED file for use by miRkwood.

Usage : ./mirkwood-bam2bed.pl -bam <input BAM file> -out <output directory> 

Dependancies : samtools
Make sure to have it installed in your PATH. For Ubuntu/Debian distributions `sudo apt-get install samtools` is enough.\n";


########## Get options
GetOptions ('bam=s' => \$bam_file,
            'out=s' => \$output_directory,
	        'help'  => \$help);


########## Validate options
if ( $help or ! -r $bam_file or ! $output_directory ){
    print $help_message;
    exit;
}

if (! -d $output_directory){
	mkdir $output_directory, 0777;
}


########## Convert BAM into SAM and filtering out unmapped reads
if ( $bam_file =~ /([^\/\\]+)[.]bam/ ){
    $sam_file = "$output_directory/$1.sam";
    $bed_file = "$output_directory/$1.bed";
    $tmp_file = "$output_directory/$1_sorted";
}
else{
    die "Non correct input file.\n$help_message";
}

system("samtools sort $bam_file $tmp_file");
system("samtools view -h -F 4 $tmp_file.bam > $sam_file");


########## Read the SAM file to store counts into a hash
open(my $SAM, '<', $sam_file) or die "ERROR while reading SAM file. Program will end prematurely.\n";

while ( <$SAM> ){

    if ( ! /^@/ ){

        chomp;

        my @line = split ( /\t/smx );

        $id = $line[0];
        $chromosome = $line[2];
        $start = $line[3] - 1;
        $sequence = $line[9];

        if ( ! exists( $counts->{$chromosome}{$start}{$sequence} ) ){
            $counts->{$chromosome}{$start}{$sequence}{'count'} = 0;
        }
        $counts->{$chromosome}{$start}{$sequence}{'count'}++;
        $counts->{$chromosome}{$start}{$sequence}{'strand'} = '+';
        if ( $line[1] eq '16' or $line[1] eq '0x10' ){
            $counts->{$chromosome}{$start}{$sequence}{'strand'} = '-';
        }

        if ( ! exists( $names->{$sequence} ) ){
            $names->{$sequence} = $id;
        }

    }

}

close $SAM;


########## Browse hash tables and print data in BED file
open(my $BED, '>', $bed_file) or die "ERROR while creating $bed_file. Program will end prematurely.\n";

foreach $chromosome ( sort (keys%$counts) ){

    foreach $start ( sort {$a <=> $b} keys%{$counts->{$chromosome}} ){

        foreach $sequence (sort (keys%{$counts->{$chromosome}{$start}}) ){
            $end = $start + length($sequence);

            print $BED "$chromosome\t";
            print $BED "$start\t";
            print $BED "$end\t";
            print $BED "$names->{$sequence}\t";
            print $BED "$counts->{$chromosome}{$start}{$sequence}{'count'}\t";
            print $BED "$counts->{$chromosome}{$start}{$sequence}{'strand'}\n";
        }

    }

}

close $BED;


unlink "$tmp_file.bam";
unlink $sam_file;
