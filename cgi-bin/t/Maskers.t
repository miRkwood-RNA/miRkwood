#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Maskers');
}
require_ok('PipelineMiRNA::Maskers');

## mask_sequence #

ok(
    my $mask_sequence_output =
      PipelineMiRNA::Maskers::mask_sequence( '123456789', 3, 6),
    'Can call mask_sequence'
);
my $mask_sequence_expected = '123NNNN89';
is ($mask_sequence_output, $mask_sequence_expected,
    'mask_sequence returs the correct value' );
