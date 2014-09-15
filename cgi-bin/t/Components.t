#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::Components');
}
require_ok('miRkwood::Components');

my @alignment_values = (
['- - none
G g match
A - gap
- - none
U U match
C c match
C c match
C U match
G g match
C c match
C c match
U U match
U U match
G g match
C c match
A a match
U U match
C c match
A a match
A a match
C c match
U U match
G g match
A a match
A a match
U U match
- - none',
'g-uccugccuugcaucaacugaau
| ||| ||||||||||||||||||
gaucccgccuugcaucaacugaau
'],
['- - none
C c match
C c match
C U match
G g match
C c match
C c match
U U match
U U match
G g match
C c match
A a match
U U match
C c match
A a match
A a match
C c match
U U match
G g match
A a match
A a match
U U match
- - none',
'ccugccuugcaucaacugaau
|| ||||||||||||||||||
cccgccuugcaucaacugaau
'],
['- - none
A U match
G c match
G g match
C c match
U U match
U U match
G g match
G g match
U U match
G g match
C c match
A a match
G g match
C g match
U U match
C c match
G g match
G g match
G g match
A a match
A a match
- - none
','ucgcuuggugcaggucgggaa
  ||||||||||| |||||||
aggcuuggugcagcucgggaa
'],
['- - none
U U match
G g match
U U match
A a match
G g match
C c match
C c match
A a match
A a match
G g match
G g match
A a match
C U match
A g match
A a match
C a match
U U match
U U match
G g match
C c match
C c match
- - none
','uguagccaaggaugaauugcc
||||||||||||  | |||||
uguagccaaggacaacuugcc
'],
['- - none
U U match
G g match
- U gap
- - none
A a match
G g match
C c match
C c match
A a match
A a match
G g match
G g match
A a match
U U match
G g match
A a match
C a match
U U match
U U match
G g match
C c match
C c match
G - gap
- - none
- - none
','uguagccaaggaugaauugcc-
|| |||||||||||| ||||| 
ug-agccaaggaugacuugccg
']
);

foreach (@alignment_values){
    my ($input, $expected) = @{$_};
    my $output = miRkwood::Components::parse_exonerate_alignment($input);
    is( $output, $expected, 'Parsing Exonerate alignment ok' );
}


my $exonerate_output_file = input_file('Components.parse_custom_exonerate_output.in');
file_exists_ok($exonerate_output_file);

#my %exonerate_expected = ('A' => 'B');
my %exonerate_expected = (
    '38-58' => [
        {
            'alignment' => 'uguagccaaggaugaauugcc
||||||||||||  | |||||
uguagccaaggacaacuugcc
',
            'name'  => 'aly-miR169a*',
            'begin_target' => '38',
            'end_target' => '58',
            'strand_target' => '+',
            'begin_query' => '21',
            'end_query' => '1',
            'def_query' => 'MIMAT0017486 Arabidopsis lyrata miR169a*:[revcomp]',
            'score' => '-3',
            'seq' => 'TGTAGCCAAGGACAACTTGCC'
        },
        {
            'alignment' => 'uguagccaaggaugaauugcc-
|| |||||||||||| ||||| 
ug-agccaaggaugacuugccg
',
            'name'  => 'aly-miR169d',
            'begin_target' => '38',
            'end_target' => '58',
            'strand_target' => '+',
            'begin_query' => '1',
            'end_query' => '21',
            'def_query' => 'MIMAT0017491 Arabidopsis lyrata miR169d',
            'score' => '-3',
            'seq'   => 'TGAGCCAAGGATGACTTGCCG'
        }
    ],
    '40-60' => [
        {
            'alignment' => 'uagccaaggaugaauugccua
 |||||||||||| ||||| |
cagccaaggaugacuugccga
',
            'name'  => 'aly-miR169a',
            'begin_target' => '40',
            'end_target' => '60',
            'strand_target' => '+',
            'begin_query' => '1',
            'end_query' => '21',
            'def_query' => 'MIMAT0017485 Arabidopsis lyrata miR169a',
            'score' => '-3',
            'seq'   => 'CAGCCAAGGATGACTTGCCGA'
        }
    ]
);

    my %output = miRkwood::Components::parse_custom_exonerate_output(
        $exonerate_output_file);

is_deeply( \%output, \%exonerate_expected, 'Parsing Exonerate output ok' );

my $CT_file_input = input_file('workspace', '1', '1', 'outB2ct_stemloop.ct' );
my $output_file = 'output_file.txt';
ok (miRkwood::Components::mask_CT_file($CT_file_input, $output_file),
    'can call mask_CT_file()' );
my $mask_CT_file_output = slurp_file($output_file);
my $mask_CT_file_expected = slurp_file(input_file('workspace', '1', '1', 'seqWithN.txt' ));
is( $mask_CT_file_output, $mask_CT_file_expected,
    'mask_CT_file ok' );

my $candidate_ct_file = input_file('workspace', '1', '1', 'outB2ct_optimal.ct' );
