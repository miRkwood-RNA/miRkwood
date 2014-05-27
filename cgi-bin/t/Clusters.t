#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Clusters');
}
require_ok('PipelineMiRNA::Clusters');

my $bamfile = input_file('Clusters.reads-Athaliana_167-ChrC.bam');
file_exists_ok($bamfile);

my $genome_file = input_file('Clusters.Athaliana_167-ChrC.fa');
file_exists_ok($genome_file);
