#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::Candidate');
}
require_ok('miRkwood::Candidate');

## Set up ##

# Necessary as the headers are fetched in the config
use miRkwood;
my $config_file = input_file('run_options.cfg');
miRkwood->CONFIG_FILE($config_file);

#new_from_serialized
my $candidate_file = input_file('candidates', '1-1.yml');
file_exists_ok($candidate_file);

my $empty_candidate = miRkwood::Candidate->new();
isa_ok($empty_candidate, 'miRkwood::Candidate');

my $candidate = miRkwood::Candidate->new_from_serialized($candidate_file);

isa_ok($candidate, 'miRkwood::Candidate');

can_ok($candidate, qw[get_identifier has_mirdup_validation get_absolute_image
                      get_relative_image get_name get_shortened_sequence_name
                      get_shortened_name candidateAsVienna candidateAsFasta
                      candidate_as_gff alternativeCandidatesAsVienna
                      make_alignments_HTML candidate_as_pseudoXML]);

is($candidate->get_identifier(), '1-1',
   'get_identifier returns the expected value');

is($candidate->has_mirdup_validation(), 1,
   'has_mirdup_validation returns the expected value');

is($candidate->get_absolute_image(),
   '/bio1/www/html/mirkwood/results/jobDec181721032013/1/1/image.png',
   'get_absolute_image returns the expected value');

is($candidate->get_relative_image(),
    '/mirkwood/results/jobDec181721032013/1/1/image.png',
    'get_relative_image returns the expected value'
);

is($candidate->get_name(), 'contig15750__34-195',
   'get_name returns the expected value');

is($candidate->get_shortened_sequence_name(), 'contig15750',
   'get_shortened_sequence_name returns the expected value');

is($candidate->get_shortened_name(), 'contig15750__34-195',
   'get_shortened_name returns the expected value');

ok( my $result4 = $candidate->candidateAsVienna(),
    'can call candidateAsVienna()' );
my $expected_file4 = input_file('Candidate.candidateAsVienna.out');
file_exists_ok($expected_file4);
my $expected4 = slurp_file($expected_file4);
is( $result4, $expected4, 'candidateAsVienna returns the expected value' );

ok( my $result5 = $candidate->candidateAsFasta(),
    'can call candidateAsFasta()' );
my $expected_file5 = input_file('Candidate.candidateAsFasta.out');
file_exists_ok($expected_file5);
my $expected5 = slurp_file($expected_file5);
is( $result5, $expected5, 'candidateAsFasta returns the expected value' );

ok( my $result_gff = $candidate->candidate_as_gff(),
    'can call candidate_as_gff()' );
my $expected_file_gff = input_file('Candidate.candidate_as_gff.out');
file_exists_ok($expected_file_gff);
my $expected_gff = slurp_file($expected_file_gff);
is( $result_gff, $expected_gff, 'candidate_as_gff returns the expected value' );

ok( my $result_pseudoXML = $candidate->candidate_as_pseudoXML() . "\n",
    'can call candidate_as_pseudoXML()' );
my $expected_file_pseudoXML = input_file('Candidate.candidate_as_pseudoXML.out');
file_exists_ok($expected_file_pseudoXML);
my $expected_pseudoXML = slurp_file($expected_file_pseudoXML);
is( $result_pseudoXML, $expected_pseudoXML, 'candidate_as_pseudoXML returns the expected value' );

ok(
    my $result6 =
      $candidate->alternativeCandidatesAsVienna(),
    'can call alternativeCandidatesAsVienna()'
);
my $expected_file6 = input_file('Candidate.alternativeCandidatesAsVienna.out');
file_exists_ok($expected_file6);
my $expected6 = slurp_file($expected_file6);
is( $result6, $expected6,
    'alternativeCandidatesAsVienna returns the expected value' );

ok(
    my $result7 = $candidate->make_alignments_HTML(),
    'can call make_alignments_HTML()'
);
my $expected_file7 = input_file('Candidate.make_alignments_HTML.out');
file_exists_ok($expected_file7);
my $expected7 = slurp_file($expected_file7);
#is( $result7, $expected7,
#    'make_alignments_HTML returns the expected value' );

#print $result7;
