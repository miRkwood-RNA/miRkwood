#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Candidate');
}
require_ok('PipelineMiRNA::Candidate');

my $candidate_dir = input_file('contig15750__34-195');

ok( my %candidate = PipelineMiRNA::Candidate->parse_candidate_information('/dummy/', $candidate_dir),
    'Can call parse_candidate_information()' );
$candidate{'alignment'} = 'None'; # Hack to get around the full path thing


my $seq    = 'gagagguucauacaugaagagaagagugcucuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucggcuuuuuuucaucauguuagaucuucuc';
my $struct = '((((((.((..((((((.((((((((.(((...((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))))))))))).))))))..)).))))))';
my %expected = ('p_value' => '0.125000',
                'amfe'    => '-37.0987654320988',
                'mfei' => '-1.00166666666667',
                'mfe' => '-60.1',
                'DNASequence' => $seq,
                'Vienna' => $struct,
                'Vienna_optimal' => $struct,
                'image' => '/dummy/image.png',
                'alignment' => 'None',
                'quality' => '2',
);
is_deeply( \%candidate, \%expected, 'retrieve candidate information ok' );

my %dummy_alignments = ();
push @{$dummy_alignments{'0-20'}}, ('A', 'B');
push @{$dummy_alignments{'0-21'}}, ('C', 'D', 'E');
push @{$dummy_alignments{'4-30'}}, ('F', 'G');
push @{$dummy_alignments{'4-31'}}, ('H');
push @{$dummy_alignments{'10-40'}}, ('I');

ok( my %merged_alignments = PipelineMiRNA::Candidate->merge_alignments(\%dummy_alignments),
    'Can call merge_alignments') ;

my %expected_merged = ();
push @{$expected_merged{'0-21'}}, ('A', 'B', 'C', 'D', 'E');
push @{$expected_merged{'4-30'}}, ('F', 'G', 'H');
push @{$expected_merged{'10-40'}}, ('I');
is_deeply( \%merged_alignments, \%expected_merged, 'merge_alignments ok' );


my %dummy_alignments2 = ();
push @{$dummy_alignments2{'37-58'}}, ('A', 'B');
push @{$dummy_alignments2{'38-58'}}, ('C', 'D', 'E');
push @{$dummy_alignments2{'39-58'}}, ('F');
push @{$dummy_alignments2{'39-59'}}, ('G', 'H');
push @{$dummy_alignments2{'39-60'}}, ('I');
push @{$dummy_alignments2{'40-58'}}, ('J', 'K');
ok( my %merged_alignments2 = PipelineMiRNA::Candidate->merge_alignments(\%dummy_alignments2),
    'Can call merge_alignments') ;

my %expected_merged2 = ();
push @{$expected_merged2{'38-58'}}, ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K');
is_deeply( \%merged_alignments2, \%expected_merged2, 'merge_alignments ok' );
