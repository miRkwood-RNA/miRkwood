#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use File::Temp;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::SequenceJob');
}
require_ok('miRkwood::SequenceJob');

my $tmp_dir = File::Temp::tempdir();
my @args = ($tmp_dir, 'sample', 'ACGATGCTGAGCTAGCGTAGCTAAT', '-');
my $sequence_job = new_ok( 'miRkwood::SequenceJob' => \@args );

my @funcs = qw(get_strand
               is_opposite_strand get_sequence_length);
can_ok( $sequence_job, @funcs );

is($sequence_job->get_sequence_length(), 25,
   'get_sequence_length returns the correct value');