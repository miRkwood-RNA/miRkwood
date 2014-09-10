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
