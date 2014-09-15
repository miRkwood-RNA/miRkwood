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
