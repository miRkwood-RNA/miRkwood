#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Utils');
}
require_ok('PipelineMiRNA::Utils');

my $fastaFile1 = input_file('fasta1.fa');
file_exists_ok($fastaFile1);

open( my $INPUT_FH, '<', $fastaFile1 ) or die "Fail to open: $!";
ok( my %tab = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH),
    'Can call parse_multi_fasta()' );
close $INPUT_FH;

my %expected = (
    '>fasta11' => 'AATGTGCCAATCCCAATGTTAACCAAAAACTAAAAAAGTGAAACGAACATTGTC',
    '>fasta12' =>
'ACTGAGATCGCAACTAATTTATTTATTCGCTCGTATAATGTATACATTAGATAGAGGCCTAGCCTCTTAGTCGAAAAGCCC',
);
is_deeply( \%tab, \%expected, 'FASTA parsing with parse_multi_fasta is ok' );

my $fastaFile2 = input_file('fasta_long_header.fa');
file_exists_ok($fastaFile2);
open( my $INPUT_FH2, '<', $fastaFile2 ) or die "Fail to open: $!";
ok( my %tab2 = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH2),
    'Can call parse_multi_fasta()' );
close $INPUT_FH2;
my %expected2 =
  ( '>fasta1' => 'AATGTGCCAATCCCAATGTTAACCAAAAACTAAAAAAGTGAAACGAACATTGTC', );
is_deeply( \%tab2, \%expected2,
    'Parsing FASTA with long header with parse_multi_fasta ok' );

my $fastaFile3 = input_file('fasta_with_pipes.fa');
file_exists_ok($fastaFile3);
open( my $INPUT_FH3, '<', $fastaFile3 ) or die "Fail to open: $!";
ok( my %tab3 = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH3),
    'Can call parse_multi_fasta()' );
close $INPUT_FH3;
my %expected3 = ( '>gi-425626932-gb-JX648278.1-' =>
      'AATGTGCCAATCCCAATGTTAACCAAAAACTAAAAAAGTGAAACGAACATTGTC', );

is_deeply( \%tab3, \%expected3,
    'Parsing FASTA with pipes using parse_multi_fasta ok' );
