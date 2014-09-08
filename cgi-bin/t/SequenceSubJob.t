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
    use_ok('miRkwood::SequenceSubJob');
}
require_ok('miRkwood::SequenceSubJob');

my $tmp_dir = File::Temp::tempdir();
my $sequence = 'ACGATGCTGAGCTAGCGTAGCTAAT';
my @args = ($tmp_dir, 'sample', $sequence , '-');
my $sequence_job = new_ok( 'miRkwood::SequenceSubJob' => \@args );

my @funcs = qw(get_strand get_sequence get_directory
               is_opposite_strand get_sequence_length);
can_ok( $sequence_job, @funcs );

is($sequence_job->get_sequence_length(), 25,
   'get_sequence_length returns the correct value');

is($sequence_job->get_sequence(), $sequence,
   'get_sequence returns the correct value');

is($sequence_job->get_directory(), $tmp_dir,
   'get_directory returns the correct value');
