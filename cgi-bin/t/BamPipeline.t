#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::BamPipeline');
}
require_ok('miRkwood::BamPipeline');

