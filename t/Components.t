#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Components');
}
require_ok('PipelineMiRNA::Components');

my $input1 = '- -
T T
G g
T T
A a
G g
C c
C c
A a
A a
G g
G g
A a
C -
- -
A T
G g
A a
C a
T T
T T
G g
C c
C c
- -';

my $expected1 = '-TGTAGCCAAGGAC-AGACTTGCC-
 ||||||||||||   || ||||| 
-TGTAGCCAAGGA--TGAATTGCC-
';
my $output1 = PipelineMiRNA::Components::parse_exonerate_alignment($input1);
is( $output1, $expected1, 'Parsing Exonerate alignment ok' );

my $exonerate_output_file = input_file('exonerate_output.yaml');
file_exists_ok($exonerate_output_file);

#my %exonerate_expected = ('A' => 'B');
my %exonerate_expected = (
    '37-58' => [
        {
            'alignment' => '-TGTAGCCAAGGACAACTTGCC-
 ||||||||||||  | ||||| 
-TGTAGCCAAGGATGAATTGCC-
',
            'name'  => 'aly-miR169a*',
            'begin_target' => '37',
            'end_target' => '58',
            'begin_query' => '21',
            'end_query' => '0',
            'score' => '-3',
            'seq'   => 'TGTAGCCAAGGACAACTTGCC'
        },
        {
            'alignment' => '-TG--AGCCAAGGATGACTTGCCG--
 ||  |||||||||||| |||||   
-TGT-AGCCAAGGATGAATTGCC---
',
            'name'  => 'aly-miR169d',
            'begin_target' => '37',
            'end_target' => '58',
            'begin_query' => '0',
            'end_query' => '21',
            'score' => '-3',
            'seq'   => 'TGAGCCAAGGATGACTTGCCG'
        }
    ],
    '39-60' => [
        {
            'alignment' => '-CAGCCAAGGATGACTTGCCGA-
  |||||||||||| ||||| | 
-TAGCCAAGGATGAATTGCCTA-
',
            'name'  => 'aly-miR169a',
            'begin_target' => '39',
            'end_target' => '60',
            'begin_query' => '0',
            'end_query' => '21',
            'score' => '-3',
            'seq'   => 'CAGCCAAGGATGACTTGCCGA'
        }
    ]
);

my %output = PipelineMiRNA::Components::parse_custom_exonerate_output(
    $exonerate_output_file);
is_deeply( \%output, \%exonerate_expected, 'Parsing Exonerate output ok' );

my $input_vienna = "
>contig15916__234-365
ccgacgguuucauauuugucuccacugugugaaaccucguagcuugaguacuguccugccuugcaucaacugaaucugaaccgauguaaaugaucugugaccgguguaggagaauuggaugaauauuguugg
(((((((((((((.....(((((((((.((.(....((((.((..((......))..)).(((((((...............)))))))))))....).))))))...))))).....))))).)))))))) (-30.70)
";
my $expected = 'contig15916__234-365'. "\t".
               'ccgacgguuucauauuugucuccacugugugaaaccucguagcuugaguacuguccugccuugcaucaacugaaucugaaccgauguaaaugaucugugaccgguguaggagaauuggaugaauauuguugg' . "\t".
               '(((((((((((((.....(((((((((.((.(....((((.((..((......))..)).(((((((...............)))))))))))....).))))))...))))).....))))).))))))))' . "\n";

open my ($input_vienna_fh), '<', \$input_vienna;
ok( my $result = PipelineMiRNA::Components::parse_RNAfold_output($input_vienna_fh),
    'Can call parse_RNAfold_output()' );
close $input_vienna_fh;
is ($result , $expected, 'Converting RNAfold output in ViennaTraited ok');
