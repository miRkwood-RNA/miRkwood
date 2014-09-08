#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::SequenceJob');
}
require_ok('miRkwood::SequenceJob');

my @args = ();
my $sequence_job = new_ok( 'miRkwood::SequenceJob' => \@args );

#my @funcs = qw();
#can_ok( $sequence_job, @funcs );
