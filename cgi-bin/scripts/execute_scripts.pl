#!/usr/bin/perl -w

use warnings;
use strict;

use FindBin;

BEGIN {
    use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
    use miRkwood;
    use miRkwood::FastaPipeline;
}


## Code ##
my ( $job_dir ) = @ARGV;

my $pipeline = miRkwood::FastaPipeline->new($job_dir);
$pipeline->run_pipeline();
