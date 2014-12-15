#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;


################################################################################
# Author : Isabelle GUIGON
# Date : 2014-11-27
# This script filters a GFF file to keep only lines corresponding to the
# features given in parameter.
################################################################################


########## Variables
my $gff_file = '';
my $output_file = '';
my $features;
my %allowed_features;
my $help;
my $help_message = "filter_gff.pl
----------
Script to filter a GFF file according to the features given in parameter.
Only the features given in parameter will be kept.

Usage : filter_gff.pl -gff <input GFF file> -out <output file> -features <CDS,tRNA,rRNA...>

Features to keep must be separated by commas.
";


########## Get options
GetOptions ('gff=s' => \$gff_file,
            'features=s' => \$features,
            'out=s' => \$output_file,
	        'help'  => \$help);


########## Validate options
if ( $help or ( ! -r $gff_file ) or $output_file eq '' ){
    print $help_message;
    exit;
}

foreach ( split( /,/xms, $features) ) {
    $allowed_features{$_} = 1;
}

print STDERR 'Output file will contain only the following features:';
foreach (keys %allowed_features){
    print STDERR " $_";
}
print STDERR "\n";


########## Filter gff
open(my $GFF, '<', $gff_file) or die "ERROR while reading $gff_file : $!";
open(my $OUT, '>', $output_file) or die "ERROR while creating $output_file : $!";

while ( <$GFF> ){
    if ( /^#/ ){
        print $OUT $_;
    }
    else{
        my @fields = split( /\t/xms );
        if ( exists($allowed_features{ $fields[2] }) ){
            print $OUT $_;
        }
    }
}

close $GFF;
close $OUT;
