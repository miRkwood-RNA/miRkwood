#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::PosterioriTests');
}
require_ok('miRkwood::PosterioriTests');

## merge_alignments

my %dummy_alignments = ();
push @{ $dummy_alignments{'0-20'} }, ( 'A', 'B' );
push @{ $dummy_alignments{'0-21'} }, ( 'C', 'D', 'E' );
push @{ $dummy_alignments{'4-30'} }, ( 'F', 'G' );
push @{ $dummy_alignments{'4-31'} },  ('H');
push @{ $dummy_alignments{'10-40'} }, ('I');

ok(
    my %merged_alignments =
      miRkwood::PosterioriTests->merge_alignments( \%dummy_alignments ),
    'Can call merge_alignments'
);

my %expected_merged = ();
push @{ $expected_merged{'0-21'} }, ( 'A', 'B', 'C', 'D', 'E' );
push @{ $expected_merged{'4-30'} }, ( 'F', 'G', 'H' );
push @{ $expected_merged{'10-40'} }, ('I');
is_deeply( \%merged_alignments, \%expected_merged, 'merge_alignments ok' );

my %dummy_alignments2 = ();
push @{ $dummy_alignments2{'37-58'} }, ( 'A', 'B' );
push @{ $dummy_alignments2{'38-58'} }, ( 'C', 'D', 'E' );
push @{ $dummy_alignments2{'39-58'} }, ('F');
push @{ $dummy_alignments2{'39-59'} }, ( 'G', 'H' );
push @{ $dummy_alignments2{'39-60'} }, ('I');
push @{ $dummy_alignments2{'40-58'} }, ( 'J', 'K' );
ok(
    my %merged_alignments2 =
      miRkwood::PosterioriTests->merge_alignments( \%dummy_alignments2 ),
    'Can call merge_alignments'
);

my %expected_merged2 = ();
push @{ $expected_merged2{'38-58'} },
  ( 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K' );
is_deeply( \%merged_alignments2, \%expected_merged2, 'merge_alignments ok' );

my %dummy_alignments3 = ();
push @{ $dummy_alignments3{'1-21'} }, ( 'A', 'B', 'C' );
push @{ $dummy_alignments3{'2-21'} }, ('D');
push @{ $dummy_alignments3{'2-23'} }, ( 'E', 'F' );
push @{ $dummy_alignments3{'2-24'} }, ('G');
push @{ $dummy_alignments3{'3-22'} }, ( 'H', 'I' );
push @{ $dummy_alignments3{'3-23'} }, ( 'J', 'K' );
push @{ $dummy_alignments3{'3-24'} }, ( 'L', 'M' );

ok(
    my %merged_alignments3 =
      miRkwood::PosterioriTests->merge_alignments( \%dummy_alignments3 ),
    'Can call merge_alignments'
);
my %expected_merged3 = ();
push @{ $expected_merged3{'1-21'} },
  ( 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M' );
is_deeply( \%merged_alignments3, \%expected_merged3,
    'merge_alignments ok on edge case' );


## mask_CT_file ##

my $CT_file_input = input_file('workspace', '1', '1', 'outB2ct_stemloop.ct' );
my $output_file = 'output_file.txt';
ok (miRkwood::PosterioriTests->mask_CT_file($CT_file_input, $output_file),
    'can call mask_CT_file()' );
my $mask_CT_file_output = slurp_file($output_file);
my $mask_CT_file_expected = slurp_file(input_file('workspace', '1', '1', 'seqWithN.txt' ));
is( $mask_CT_file_output, $mask_CT_file_expected,
    'mask_CT_file ok' );

my $candidate_ct_file = input_file('workspace', '1', '1', 'outB2ct_optimal.ct' );
    