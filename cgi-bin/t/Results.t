#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::Results');
}
require_ok('miRkwood::Results');

my $candidates_dir = input_file('YML');

ok(
    my %results = miRkwood::Results->deserialize_results($candidates_dir),
    'Can call deserialize_results'
);

my @keys       = keys %results;
my $identifier = $keys[0];
is( $identifier, '1-1', 'deserialize_results correctly deserialized data' );

ok( my $has_candidates = miRkwood::Results->has_candidates( \%results ),
    'can call has_candidates' );
ok( $has_candidates, 'has_candidates ok' );

ok(
    my $number_of_results_output =
      miRkwood::Results->number_of_results( \%results ),
    'can call number_of_results'
);

is( $number_of_results_output, 1,
    'number_of_results return the correct value' );
