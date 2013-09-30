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
T T match
C c match
C c match
C T match
G g match
C c match
C c match
T T match
T T match
G g match
C c match
A a match
T T match
C c match
A a match
A a match
C c match
T T match
G g match
A a match
A a match
T T match
- - none',
'GATCCCGCCTTGCATCAACTGAAT
| ||| ||||||||||||||||||
g-TccTgccTTgcaTcaacTgaaT
'],
['- - none
C c match
C c match
C T match
G g match
C c match
C c match
T T match
T T match
G g match
C c match
A a match
T T match
C c match
A a match
A a match
C c match
T T match
G g match
A a match
A a match
T T match
- - none',
'CCCGCCTTGCATCAACTGAAT
|| ||||||||||||||||||
ccTgccTTgcaTcaacTgaaT
'],
['- - none
A T match
G c match
G g match
C c match
T T match
T T match
G g match
G g match
T T match
G g match
C c match
A a match
G g match
C g match
T T match
C c match
G g match
G g match
G g match
A a match
A a match
- - none
','AGGCTTGGTGCAGCTCGGGAA
  ||||||||||| |||||||
TcgcTTggTgcaggTcgggaa
'],
['- - none
T T match
G g match
T T match
A a match
G g match
C c match
C c match
A a match
A a match
G g match
G g match
A a match
C T match
A g match
A a match
C a match
T T match
T T match
G g match
C c match
C c match
- - none
','TGTAGCCAAGGACAACTTGCC
||||||||||||  | |||||
TgTagccaaggaTgaaTTgcc
'],
['- - none
T T match
G g match
- T gap
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
T T match
G g match
A a match
C a match
T T match
T T match
G g match
C c match
C c match
G - gap
- - none
- - none
','TG-AGCCAAGGATGACTTGCCG
|| |||||||||||| ||||| 
TgTagccaaggaTgaaTTgcc-
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
            'alignment' => 'TGTAGCCAAGGACAACTTGCC
||||||||||||  | |||||
TgTagccaaggaTgaaTTgcc
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
            'alignment' => 'TG-AGCCAAGGATGACTTGCCG
|| |||||||||||| ||||| 
TgTagccaaggaTgaaTTgcc-
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
            'alignment' => 'CAGCCAAGGATGACTTGCCGA
 |||||||||||| ||||| |
TagccaaggaTgaaTTgccTa
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
