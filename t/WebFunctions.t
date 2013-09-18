#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::WebFunctions');
}
require_ok('PipelineMiRNA::WebFunctions');

my $candidate_dir = input_file('contig15750__34-195');

ok( my %results = PipelineMiRNA::WebFunctions->actual_retrieve_candidate_information('/dummy/', $candidate_dir),
    'Can call actual_retrieve_candidate_information()' );
$results{'alignment'} = 'None'; # Hack to get around the full path thing


my $seq    = 'gagagguucauacaugaagagaagagugcucuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucggcuuuuuuucaucauguuagaucuucuc';
my $struct = '((((((.((..((((((.((((((((.(((...((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))))))))))).))))))..)).))))))';
my %expected = ('p_value' => '0.125000',
                'amfe'    => '-37.0987654320988',
                'mfei' => '-1.00166666666667',
                'mfe' => '-60.1',
                'DNASequence' => $seq,
                'Vienna' => $struct,
                'image' => '/dummy/image.png',
                'alignment' => 'None',
                'quality' => '2',
);
is_deeply( \%results, \%expected, 'retrieve candidate information ok' );
