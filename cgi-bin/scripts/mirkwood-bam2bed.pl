#!/usr/bin/perl -w
use strict;
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
my $bam_file = "";
my $sam_file;
my $bed_file;
my $counts;
my $help;
my $help_message = "mirkwood-bam2bed.pl
----------
Script to convert a BAM into a BED file for use by miRkwood.

Usage : ./mirkwood-bam2bed.pl -bam <input BAM file> -bed <output BED file> 

Dependancies : samtools
Make sure to have it installed in your PATH. For Ubuntu/Debian distributions `sudo apt-get install samtools` is enough.\n";


########## Get options
GetOptions ('bam=s' => \$bam_file,
            'bed=s' => \$bed_file,
	        'help'  => \$help);
        
        
########## Validate options
if ( $help or ! -r $bam_file or ! $bed_file ){
    print $help_message;
    exit;
}


########## Convert BAM into SAM and filtering out unmapped reads
if ( $bam_file =~ /(.*)\.bam/ ){
    $sam_file = "$1.sam";
}
else{
    die "Non correct input file.\n$help_message";
}

system("samtools view -h -F 4 $bam_file > $sam_file");


########## Read the SAM file to store counts into a hash
open(SAM, $sam_file) or die "ERROR while reading SAM file. Program will end prematurely.\n";

while ( <SAM> ){

    if ( ! /^@/ ){
        
        chomp;
        
        my @line = split ("\t");
        if ( ! exists( $counts->{$line[9]} ) ){
            $counts->{$line[9]}{"id"} = $line[0];
            $counts->{$line[9]}{"count"} = 0;
        }
        $counts->{$line[9]}{"count"}++;
        
    }
    
}

close SAM;


########## Read the SAM file to convert each line into BED
open(SAM, $sam_file) or die "ERROR while reading $sam_file. Program will end prematurely.\n";
open(BED, ">$bed_file") or die "ERROR while creating $bed_file. Program will end prematurely.\n";

while ( <SAM> ){
    
    if ( ! /^@/ ){
        
        chomp;   
        
        my @line = split ("\t");
        
        if ( $counts->{$line[9]}{"id"} eq $line[0] ){
            
            my $chrom  = $line[2];
            my $start  = $line[3] -1;
            my $end    = $start + length($line[9]);
            my $name   = $line[0];
            my $score  = $counts->{$line[9]}{'count'};
            my $strand = "+";
            if ( $line[1] eq "16" or $line[1] eq "0x10" ){
                $strand = "-";
            }
            
            print BED "$chrom\t";
            print BED "$start\t";
            print BED "$end\t";
            print BED "$name\t";
            print BED "$score\t";
            print BED "$strand\n";
            
        }
        
    } 
    
}

close SAM;
close BED;

