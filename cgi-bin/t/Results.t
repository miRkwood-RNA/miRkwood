#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Results');
}
require_ok('PipelineMiRNA::Results');

my $candidates_dir = input_file('candidates');

ok(
    my %results = PipelineMiRNA::Results->deserialize_results($candidates_dir),
    'Can call deserialize_results'
);

my @keys       = keys %results;
my $identifier = $keys[0];
is( $identifier, '1-1', 'deserialize_results correctly deserialized data' );

ok( my $has_candidates = PipelineMiRNA::Results->has_candidates( \%results ),
    'can call has_candidates' );
ok( $has_candidates, 'has_candidates ok' );
