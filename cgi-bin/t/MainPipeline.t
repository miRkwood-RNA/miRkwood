#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;
use Data::Dumper;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::MainPipeline');
}
require_ok('PipelineMiRNA::MainPipeline');


my @is_to_merge_values = (
    [ [ 2, 7,  1, 9  ],  1],
    [ [ 2, 12, 1, 9  ],  1],
    [ [ 5, 12,  1, 9 ],  ''],
);

foreach my $couple (@is_to_merge_values) {
    my @input = @{ @{$couple}[0] };
    my $expected = @{$couple}[1];
    my $is_to_merge_res =
      PipelineMiRNA::MainPipeline::is_to_merge(@input);
    is( $is_to_merge_res, $expected, "is_to_merge (@input) --> $expected ok" );
}
