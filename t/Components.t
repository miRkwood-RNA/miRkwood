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
    my $output = PipelineMiRNA::Components::parse_exonerate_alignment($input);
    is( $output, $expected, 'Parsing Exonerate alignment ok' );
}


my $exonerate_output_file = input_file('exonerate_output.yaml');
file_exists_ok($exonerate_output_file);

#my %exonerate_expected = ('A' => 'B');
my %exonerate_expected = (
    '37-58' => [
        {
            'alignment' => 'uguagccaaggaugaauugcc
||||||||||||  | |||||
uguagccaaggacaacuugcc
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
            'alignment' => 'uguagccaaggaugaauugcc-
|| |||||||||||| ||||| 
ug-agccaaggaugacuugccg
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
            'alignment' => 'uagccaaggaugaauugccua
 |||||||||||| ||||| |
cagccaaggaugacuugccga
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
