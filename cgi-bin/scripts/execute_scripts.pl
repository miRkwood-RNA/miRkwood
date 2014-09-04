#!/usr/bin/perl -w

use warnings;
use strict;

use FindBin;

BEGIN {
    use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
    use miRkwood;
    use miRkwood::Pipeline;
}


## Code ##
my ( $job_dir ) = @ARGV;

my $pipeline = miRkwood::Pipeline->new($job_dir);
$pipeline->run_pipeline();
