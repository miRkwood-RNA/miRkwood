#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::Components');
}
require_ok('miRkwood::Components');

my $CT_file_input = input_file('workspace', '1', '1', 'outB2ct_stemloop.ct' );
my $output_file = 'output_file.txt';
ok (miRkwood::Components::mask_CT_file($CT_file_input, $output_file),
    'can call mask_CT_file()' );
my $mask_CT_file_output = slurp_file($output_file);
my $mask_CT_file_expected = slurp_file(input_file('workspace', '1', '1', 'seqWithN.txt' ));
is( $mask_CT_file_output, $mask_CT_file_expected,
    'mask_CT_file ok' );

my $candidate_ct_file = input_file('workspace', '1', '1', 'outB2ct_optimal.ct' );
