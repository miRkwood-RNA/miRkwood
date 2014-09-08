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
    use_ok('miRkwood::Pipeline');
}
require_ok('miRkwood::Pipeline');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

my $sequences_file = input_file('Pipeline.input_sequences.fa');
file_exists_ok($sequences_file);

my $pipeline = miRkwood::Pipeline->new($tmp_dir);

isa_ok($pipeline, 'miRkwood::Pipeline');

