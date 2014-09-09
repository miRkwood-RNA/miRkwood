#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::CandidateJob');
}
require_ok('miRkwood::CandidateJob');

my $candidate_dir = input_file('workspace', '1', '1');

my @args = ($candidate_dir);
my $candidate_job = new_ok( 'miRkwood::CandidateJob' => \@args );

my @funcs = qw(get_directory);
can_ok( $candidate_job, @funcs );

use miRkwood;


ok(
    my $candidate_object_from_dir =
      $candidate_job->make_candidate_from_directory(),
    'Can call make_candidate_from_directory()'
);
isa_ok($candidate_object_from_dir, 'miRkwood::Candidate');

ok(
    my $pseudo_candidate_contents =
      $candidate_job->parse_candidate_information($candidate_dir),
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

my $cfg_file = input_file('run_options.cfg');
miRkwood->CONFIG_FILE($cfg_file);

my $dummy_dir = input_file();
ok( my $result3 = $candidate_job->get_candidate_information_from_run(),
    'can call get_candidate_information_from_run()');
isa_ok($result3, 'miRkwood::Candidate');
