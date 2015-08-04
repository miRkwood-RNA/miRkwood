#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;
use File::Temp;
use File::Spec;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::CandidateHandler');
}
require_ok('miRkwood::CandidateHandler');


dies_ok {
    my $candidate =
      miRkwood::CandidateHandler->retrieve_candidate_information( 'a', 'b', 'c' );
} 'retrieve_candidate_information with wrong parameters dies as expected';

my $candidate_input = input_file('');

ok(
    my $candidate = miRkwood::CandidateHandler->retrieve_candidate_information(
        $candidate_input, '1-1'
    ),
    'Can call retrieve_candidate_information'
);
isa_ok($candidate, 'miRkwood::Candidate');

ok( my $candidate_filename = miRkwood::CandidateHandler->make_candidate_filename('id'),
    'Can call make_candidate_filename');
is( $candidate_filename, 'id.yml',
    'make_candidate_filename returns the expected value');

ok( my $candidate_filepath = miRkwood::CandidateHandler->get_candidate_filepath('a', 'b'),
    'Can call get_candidate_filepath');
is( $candidate_filepath, 'a/YML/b.yml',
    'get_candidate_filepath returns the expected value');

my $tmp_dir = File::Temp::tempdir();

ok(  miRkwood::CandidateHandler->serialize_candidate_information($tmp_dir, $candidate),
     'can call serialize_candidate_information');
my $tmp_file = File::Spec->catfile($tmp_dir, '1-1.yml');
file_exists_ok($tmp_file,
               "Serialized file exists (in $tmp_file)");
