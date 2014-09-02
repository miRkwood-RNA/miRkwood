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

use miRkwood;

my $candidate_dir = input_file('candidate1');

ok(
    my $candidate_object_from_dir =
      miRkwood::CandidateHandler->make_candidate_from_directory($candidate_dir),
    'Can call make_candidate_from_directory()'
);
isa_ok($candidate_object_from_dir, 'miRkwood::Candidate');

ok(
    my $pseudo_candidate_contents =
      miRkwood::CandidateHandler->parse_candidate_information($candidate_dir),
    'Can call parse_candidate_information()'
);

my $seq    = 'gagagguucauacaugaagagaagagugcucuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucggcuuuuuuucaucauguuagaucuucuc';
my $struct = '((((((.((..((((((.((((((((.(((...((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))))))))))).))))))..)).))))))';
my %expected = ('name' => 'contig15750',
#                'position' => '34-195',
                'start_position' => '34',
                'end_position' => '195',
                'shuffles' => '0.125000',
                'amfe'    => '-37.10',
                'mfei' => '-1.00',
                'mfe' => '-60.1',
                'DNASequence' => $seq,
                'Vienna' => $struct,
                'Vienna_optimal' => $struct,
                'alignment_existence' => 1,
                'strand'  => '+',
);
delete $pseudo_candidate_contents->{'alignments'};
delete $pseudo_candidate_contents->{'mirdup_validation'};
is_deeply( $pseudo_candidate_contents, \%expected, 'parse candidate information ok' );


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
is( $candidate_filepath, 'a/candidates/b.yml',
    'get_candidate_filepath returns the expected value');

my $tmp_dir = File::Temp::tempdir();

ok(  miRkwood::CandidateHandler->serialize_candidate_information($tmp_dir, $candidate),
     'can call serialize_candidate_information');
my $tmp_file = File::Spec->catfile($tmp_dir, '1-1.yml');
file_exists_ok($tmp_file,
               "Serialized file exists (in $tmp_file)");

my $cfg_file = input_file('run_options.cfg');
miRkwood->CONFIG_FILE($cfg_file);

my $dummy_dir = input_file();
ok( my $result3 = miRkwood::CandidateHandler->get_candidate_information_from_run($dummy_dir, '1', '1'),
    'can call get_candidate_information_from_run()');
isa_ok($result3, 'miRkwood::Candidate');
