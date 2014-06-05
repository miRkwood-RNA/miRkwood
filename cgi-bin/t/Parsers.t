#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Parsers');
}
require_ok('PipelineMiRNA::Parsers');

# Vienna line

my $input1    = '(((((.(((.(((((.))....))).)))...))))) (-16.93)';
my @expected1 = [ '(((((.(((.(((((.))....))).)))...)))))', '-16.93' ];
my @result1   = PipelineMiRNA::Parsers::parse_Vienna_line($input1);
is_deeply( \@result1, @expected1, 'Standard line correctly parsed' );

my $input2    = '(((((.(((.(((((.))....))).)))...))))) ( -8.58)';
my @expected2 = [ '(((((.(((.(((((.))....))).)))...)))))', '-8.58' ];
my @result2   = PipelineMiRNA::Parsers::parse_Vienna_line($input2);
is_deeply( \@result2, @expected2,
    'Line with whitespace in energy correctly parsed' );

# P_value

my $p_value_file = input_file('Parsers.pvalue.txt');
file_exists_ok($p_value_file);
my $expected3 = '0.125000';
my $result3   = PipelineMiRNA::Parsers::parse_pvalue($p_value_file);
is( $result3, $expected3, 'p_value file is correctly parsed' );

# MFEI

my $mfei_file = input_file('Parsers.outMFEI.txt');
file_exists_ok($mfei_file);
my @expected4 = [ '-1.00166666666667', '-60.1', '-37.0987654320988' ];
my @result4 = PipelineMiRNA::Parsers::parse_mfei($mfei_file);
is_deeply( \@result4, @expected4, 'MFEI file is correctly parsed' );

# Vienna

my $vienna_file1 = input_file('Parsers.outViennaTraited.txt');
my $vienna_file2 = input_file('Parsers.outViennaTraited2.txt');
file_exists_ok($vienna_file1);
file_exists_ok($vienna_file2);
my @expected5 = ['gagagguucauacaugaagagaagagugcucuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucggcuuuuuuucaucauguuagaucuucuc',
                 '((((((.((..((((((.((((((((.(((...((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))))))))))).))))))..)).))))))' ];
my @result5   = PipelineMiRNA::Parsers::parse_vienna($vienna_file1);
is_deeply( \@result5, @expected5, 'Vienna file "all in line" is correctly parsed' );

my @result6   = PipelineMiRNA::Parsers::parse_vienna($vienna_file2);
is_deeply( \@result6, @expected5, 'Vienna file "with breakspace" is correctly parsed' );

# Selfcontain

my $selfcontain_file = input_file('Parsers.selfContain.txt');
file_exists_ok($selfcontain_file);
my $expected7 = '0';
my $result7   = PipelineMiRNA::Parsers::parse_selfcontain($selfcontain_file);
is( $result7, $expected7, 'Selfcontain file is correctly parsed' );

# RNAFold

my @expected8 = ('contig15750__34-195',
               'gagagguucauacaugaagagaagagugcucuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucggcuuuuuuucaucauguuagaucuucuc',
               '((((((.((..((((((.((((((((.(((...((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))))))))))).))))))..)).))))))',
               '-60.10',
               );

my $rnafold_file = input_file('candidate1', 'outRNAFold_stemloop.txt');
ok( my @result8 = PipelineMiRNA::Parsers::parse_RNAfold_output($rnafold_file),
    'Can call parse_RNAfold_output()' );
is_deeply (\@result8 , \@expected8, 'Converting RNAfold output in ViennaTraited ok');

my $alternative_candidates_file = input_file('Parsers.alternativeCandidates.txt');
file_exists_ok($alternative_candidates_file);
ok( my %alternatives = PipelineMiRNA::Parsers::parse_alternative_candidates_file($alternative_candidates_file),
    'Can callparse_alternative_candidates_file');

my %expected6;
$expected6{'contig15750__64-168'} = {
  'mfei' => '-0.879487179487179',
  'sequence' => 'cuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucgg',
  'structure' => '((.((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))'
};
$expected6{'contig15750__34-195'} = {
  'mfei' => '-1.00166666666667',
  'sequence' => 'gagagguucauacaugaagagaagagugcucuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucggcuuuuuuucaucauguuagaucuucuc',
  'structure' => '((((((.((..((((((.((((((((.(((...((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))))))))))).))))))..)).))))))'
};

is_deeply( \%alternatives, \%expected6, 'Alternative candidate file is correctly parsed' );

my $blastout_file = input_file('Parsers.blastx.out');
file_exists_ok($blastout_file);
ok(
    my $index_blast_output_result =
      PipelineMiRNA::Parsers::index_blast_output($blastout_file),
    'Can call index_blast_output'
);
my $index_blast_output_expected = { 'arabidopsis_filtered' => 1 };
is_deeply( $index_blast_output_result, $index_blast_output_expected,
    'BLAST output file correctly parsed' );


my $tRNAscanSE_file = input_file('Parsers.tRNAscanSE.out');
ok(
    my @parse_tRNAscanSE_output =
      PipelineMiRNA::Parsers::parse_tRNAscanSE_output($tRNAscanSE_file),
    'Can call parse_tRNAscanSE_output'
);
my @parse_tRNAscanSE_expected_1 = ('C28G1', 9738, 9809);
my @parse_tRNAscanSE_expected_2 = ('C28G1', 20346, 20417);
my @parse_tRNAscanSE_expected_3 = ('CELF22B7', 12619, 12738);
my @parse_tRNAscanSE_expected = (\@parse_tRNAscanSE_expected_1, \@parse_tRNAscanSE_expected_2, \@parse_tRNAscanSE_expected_3);
is_deeply( \@parse_tRNAscanSE_output, \@parse_tRNAscanSE_expected,
    'tRNAscanSE output file correctly parsed' );
