#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Candidate');
}
require_ok('PipelineMiRNA::Candidate');

my $candidate_dir = input_file('candidate1');

ok(
    my %pseudo_candidate =
      PipelineMiRNA::Candidate->parse_candidate_information($candidate_dir),
    'Can call parse_candidate_information()'
);

my $seq    = 'gagagguucauacaugaagagaagagugcucuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucggcuuuuuuucaucauguuagaucuucuc';
my $struct = '((((((.((..((((((.((((((((.(((...((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))))))))))).))))))..)).))))))';
my %expected = ('name' => 'contig15750',
                'position' => '34-195',
                'p_value' => '0.125000',
                'amfe'    => '-37.10',
                'mfei' => '-1.00',
                'mfe' => '-60.1',
                'DNASequence' => $seq,
                'Vienna' => $struct,
                'Vienna_optimal' => $struct,
                'alignment' => 2,
                'alignment_existence' => 1,
                'quality' => '3',
                'strand'  => '+',
);
delete $pseudo_candidate{'alignments'};
delete $pseudo_candidate{'mirdup_validation'};
is_deeply( \%pseudo_candidate, \%expected, 'parse candidate information ok' );


dies_ok {
    my %candidate =
      PipelineMiRNA::Candidate->retrieve_candidate_information( 'a', 'b', 'c' );
} 'retrieve_candidate_information with wrong parameters dies as expected';

my $candidate_input = input_file('');

ok(
    my %candidate = PipelineMiRNA::Candidate->retrieve_candidate_information(
        $candidate_input, '1-1'
    ),
    'Can call retrieve_candidate_information'
);

ok( my $result1 = PipelineMiRNA::Candidate->compute_quality( \%candidate ),
    'can call compute_quality()' );
is( $result1, 3, 'compute_quality returns the expected value' );

ok( my $has_mirdup_validation_result = PipelineMiRNA::Candidate->has_mirdup_validation( \%candidate ),
    'can call has_mirdup_validation()' );
is( $has_mirdup_validation_result, 1, '$has_mirdup_validation returns the expected value' );

ok( my $result2 = PipelineMiRNA::Candidate->get_absolute_image( \%candidate ),
    'can call get_absolute_image()' );
is(
    $result2,
    '/var/www/arn/results/jobDec181721032013/1/1/image.png',
    'get_absolute_image returns the expected value'
);

ok( my $result3 = PipelineMiRNA::Candidate->get_relative_image( \%candidate ),
    'can call get_relative_image()' );
is(
    $result3,
    '/arn/results/jobDec181721032013/1/1/image.png',
    'get_relative_image returns the expected value'
);

ok( my $result4 = PipelineMiRNA::Candidate->candidateAsVienna( \%candidate ),
    'can call candidateAsVienna()' );
my $expected_file4 = input_file('candidateAsVienna.out');
file_exists_ok($expected_file4);
my $expected4 = slurp_file($expected_file4);
is( $result4, $expected4, 'candidateAsVienna returns the expected value' );

ok( my $result5 = PipelineMiRNA::Candidate->candidateAsFasta( \%candidate ),
    'can call candidateAsFasta()' );
my $expected_file5 = input_file('candidateAsFasta.out');
file_exists_ok($expected_file5);
my $expected5 = slurp_file($expected_file5);
is( $result5, $expected5, 'candidateAsFasta returns the expected value' );

ok( my $result_gff = PipelineMiRNA::Candidate->candidate_as_gff( \%candidate ),
    'can call candidate_as_gff()' );
my $expected_file_gff = input_file('candidate_as_gff.out');
file_exists_ok($expected_file_gff);
my $expected_gff = slurp_file($expected_file_gff);
is( $result_gff, $expected_gff, 'candidate_as_gff returns the expected value' );


ok(
    my $result6 =
      PipelineMiRNA::Candidate->alternativeCandidatesAsVienna( \%candidate ),
    'can call alternativeCandidatesAsVienna()'
);
my $expected_file6 = input_file('alternativeCandidatesAsVienna.out');
file_exists_ok($expected_file6);
my $expected6 = slurp_file($expected_file6);
is( $result6, $expected6,
    'alternativeCandidatesAsVienna returns the expected value' );

ok(
    my $result7 = PipelineMiRNA::Candidate->make_alignments_HTML(
        \%candidate, 'a', 'b', 'c'
    ),
    'can call make_alignments_HTML()'
);
my $expected_file7 = input_file('make_alignments_HTML.out');
file_exists_ok($expected_file7);
my $expected7 = slurp_file($expected_file7);
#is( $result7, $expected7,
#    'make_alignments_HTML returns the expected value' );

#print $result7;
