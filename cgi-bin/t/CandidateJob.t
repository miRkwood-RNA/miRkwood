#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use YAML;
use File::Temp;
use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::CandidateJob');
}
require_ok('miRkwood::CandidateJob');

my $candidate_dir = input_file('workspace', '1', '1');

my @args = ($candidate_dir, 'my_seq', 'my_id');
my $real_candidate_job = new_ok( 'miRkwood::CandidateJob' => \@args );

my @funcs = qw(get_directory);
can_ok( $real_candidate_job, @funcs );

is( $real_candidate_job->get_directory(), $candidate_dir,
    'get_directory() returns the correct value');

my $tmp_dir = File::Temp::tempdir();
my $proto_candidate_file = input_file('CandidateJob.protocandidate.in.yml');
my $proto_candidate = YAML::LoadFile($proto_candidate_file);

my @args2 = ($tmp_dir, 'my_seq', 'my_id', $proto_candidate, undef);
my $dummy_candidate_job = new_ok( 'miRkwood::CandidateJob' => \@args2 );

my $expected_write_RNAFold_stemloop_output = slurp_file(input_file('CandidateJob.write_RNAFold_stemloop_output.out'));

my $output_write_RNAFold_stemloop_output_file = $dummy_candidate_job->write_RNAFold_stemloop_output();
file_exists_ok($output_write_RNAFold_stemloop_output_file);
my $output_write_RNAFold_stemloop_output = slurp_file($output_write_RNAFold_stemloop_output_file);
is( $expected_write_RNAFold_stemloop_output, $output_write_RNAFold_stemloop_output,
    'write_RNAFold_stemloop_output returns the expected value' );


## merge_alignments

my %dummy_alignments = ();
push @{ $dummy_alignments{'0-20'} }, ( 'A', 'B' );
push @{ $dummy_alignments{'0-21'} }, ( 'C', 'D', 'E' );
push @{ $dummy_alignments{'4-30'} }, ( 'F', 'G' );
push @{ $dummy_alignments{'4-31'} },  ('H');
push @{ $dummy_alignments{'10-40'} }, ('I');

ok(
    my %merged_alignments =
      miRkwood::CandidateJob->merge_alignments( \%dummy_alignments ),
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
      miRkwood::CandidateJob->merge_alignments( \%dummy_alignments2 ),
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
      miRkwood::CandidateJob->merge_alignments( \%dummy_alignments3 ),
    'Can call merge_alignments'
);
my %expected_merged3 = ();
push @{ $expected_merged3{'1-21'} },
  ( 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M' );
is_deeply( \%merged_alignments3, \%expected_merged3,
    'merge_alignments ok on edge case' );
