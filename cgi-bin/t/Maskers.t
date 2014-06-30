#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::Maskers');
}
require_ok('miRkwood::Maskers');

## mask_sequence #

ok(
    my $mask_sequence_output =
      miRkwood::Maskers::mask_sequence( '123456789', 3, 6),
    'Can call mask_sequence'
);
my $mask_sequence_expected = '123NNNN89';
is ($mask_sequence_output, $mask_sequence_expected,
    'mask_sequence returs the correct value' );

## mask_sequences #

my %masker_11 = ( 'start' => 5, 'end' => 10 );
my %masker_12 = ( 'start' => 17, 'end' => 19 );
my @masker_1  = ( \%masker_11, \%masker_12 );
my %masker = ( 'A' => \@masker_1 );

my @input_sequence1 = ('A', '12345678901234567890');
my @input_sequence2 = ('B', '12345678901234567890');
my @input_sequences = ( \@input_sequence1, \@input_sequence2 );
ok (
    my @result_sequences =
      miRkwood::Maskers::mask_sequences(\%masker, @input_sequences),
      'Can call mask_sequences'
);

my @expected_sequence1 = ('A', '12345NNNNNN234567NNN');
my @expected_sequence2 = ('B', '12345678901234567890');
my @expected_sequences = ( \@expected_sequence1, \@expected_sequence2 );

is_deeply (\@result_sequences, \@expected_sequences,
    'mask_sequence returs the correct value' );
