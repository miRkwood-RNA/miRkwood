#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;


################################################################################
# Author : Isabelle GUIGON
# Date : 2014-11-27
# This script filters a GFF file to keep only lines corresponding to the
# allowed features, which are CDS, tRNA, rRNA and snoRNA.
################################################################################


########## Variables
my $gff_file = "";
my $output_file = "";
my $help;
my $help_message = "filter_gff.pl
----------
Script to filter a GFF file according to the features.
Only CDS, tRNA, rRNA and snoRNA will be kept.

Usage : filter_gff.pl -gff <input GFF file> -out <output file>
";

my %allowed_features = ( "CDS"    => 1,
                         "tRNA"   => 1,
                         "rRNA"   => 1,
                         "snoRNA" => 1 );


########## Get options
GetOptions ('gff=s' => \$gff_file,
            'out=s' => \$output_file,
	        'help'  => \$help);
        
        
########## Validate options
if ( $help or ! -r $gff_file or ! -r $output_file ){
    print $help_message;
    exit;
}


########## Filter gff
open(GFF, $gff_file) or die "ERROR while reading $gff_file : $!";
open(OUT, ">$output_file") or die "ERROR while creating $output_file : $!";

while ( <GFF> ){
    if ( /^#/ ){
        print OUT $_;
    }
    else{
        my @fields = split("\t");
        if ( exists($allowed_features{ $fields[2] }) ){
            print OUT $_;
        }
    }
}

close GFF;
close OUT;
