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

my $fastaFile4 = input_file('fasta2.fa');
file_exists_ok($fastaFile4);
open( my $INPUT_FH4, '<', $fastaFile4 ) or die "Fail to open: $!";
ok( my %tab4 = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH4),
    'Can call parse_multi_fasta()' );
close $INPUT_FH4;
my %expected4 = ( '>contig15750' =>
                  'aatgagtaagataaattgctaattaaatgcgacgagaggttcatacatgaagagaagagtgctcttattatgtagccaaggatgaattgcctaatgacagctcaagtcgtttaaaaaacgactctttgttggtttattaggcgttcatttcttgactgacttaatcggctttttttcatcatgttagatcttctcaacttgttacgagcatatcgttcaatattttcatagtcttcttgtaatatgactttgtcaagtcatttcatatagctacttatgtgtagctattattgtcataattattatatagattatatacttaaagagagacttgtaagggatttaagatgtttagataatcatgtaacattcttgtcaagttatgatcaagcattat',
                  '>contig15916' =>
                  'aaaaaacctcacatacagcccccgtatctctctctctctataattgataggctattttcttctctctctagaaatgagcttacatggcatgcagatccattgcttatttataggtatagatacagcagatatatattatttattcatatatgtgtatcgaggtatcggaagaagaaattttcattgttacggcggttttctgattcgcttggtgcaggtcgggaacggcttggccgacggtttcatatttgtctccactgtgtgaaacctcgtagcttgagtactgtcctgccttgcatcaactgaatctgaaccgatgtaaatgatctgtgaccggtgtaggagaattggatgaatattgttggagat'
                  );
is_deeply( \%tab4, \%expected4,
           'Parsing FASTA with pipes using parse_multi_fasta ok' );

##################################################################
diag('Testing find_matching_count()');
is( PipelineMiRNA::Utils::find_matching_count('()..'), 1, 'one - matching left' );
is( PipelineMiRNA::Utils::find_matching_count('..()'), 1, 'one - matching right' );
is( PipelineMiRNA::Utils::find_matching_count('.().'), 1, 'one - matching middle' );
is( PipelineMiRNA::Utils::find_matching_count('(..)'), 1, 'one - matching across' );
is( PipelineMiRNA::Utils::find_matching_count('(.().)'), 2, 'two - matching across' );

##################################################################
diag('Testing make_loop()');
use Data::Dumper;
ok( my @res5 = PipelineMiRNA::Utils::make_loop('123'),
    'Can call make_loop() with 3 elements');
#print Dumper(@res5);
is_deeply( \@res5, [ [], ['1'], ['2'], ['3'], [] ], 'make_loop() with 3 elements ok' );

ok( my @res6 = PipelineMiRNA::Utils::make_loop('1234'),
    'Can call make_loop() with 4 elements');
#print Dumper(@res6);
is_deeply( \@res6, [ ['1'], [' ', '2'], [' '], [' ', '3'], ['4'] ], 'make_loop() with 4 elements ok' );

ok( my @res7 = PipelineMiRNA::Utils::make_loop('12345'),
    'Can call make_loop() with 5 elements');
is_deeply( \@res7, [ ['1'], [' ', '2'], [' ', '3'], [' ', '4'], ['5'] ], 'make_loop() with 5 elements ok' );

my $sequence1 = "cuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucgg";
my $vienna1   = "((.((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))";

ok( my $result8 = PipelineMiRNA::Utils::make_ASCII_viz($sequence1, $vienna1),
    'Can call make_ASCII_viz()');
my $expected8 = "  u    u  -  c      -     u         --    uca-        
cu auua gu ag caagga ugaau gccuaauga  cagc    agucguuua
|| |||| || || |||||| ||||| |||||||||  ||||    ||||||||a
gg uaau ca uc guucuu acuug cggauuauu  guug    ucagcaaaa
  c    u  g  a      u     -         ug    uuuc        
";
is( $result8, $expected8, 'make_ASCII_viz returns a correct hairpin');