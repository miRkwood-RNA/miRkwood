#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::Paths');
}
require_ok('miRkwood::Paths');

ok( my $result1 = miRkwood::Paths->get_job_config_path('a/b'),
    'can call get_job_config_path()');
is( $result1, 'a/b/run_options.cfg',
    'get_job_config_path returns expected value');

ok( my $result2 = miRkwood::Paths->get_candidate_paths('a/b', 'c', 'd'),
    'can call get_candidate_paths()');
is( $result2, 'a/b/workspace/c/d',
    'get_candidate_paths returns expected value');
