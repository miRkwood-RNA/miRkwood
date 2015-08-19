#!/usr/bin/perl -w

use warnings;
use strict;

use FindBin;

BEGIN {
    use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
    use miRkwood;
}


## Code ##
#~ my ( $job_dir ) = @ARGV;
my $mode = $ARGV[0];
my $job_dir = $ARGV[1];
my $pipeline;

if ( $mode eq 'smallRNAseq' ){
    use miRkwood::BEDPipeline;
    my $localBED = $ARGV[2];
    my $genome   = $ARGV[3];
    $pipeline = miRkwood::BEDPipeline->new($job_dir, $localBED, $genome);
}
else{
    use miRkwood::FastaPipeline;
    $pipeline = miRkwood::FastaPipeline->new($job_dir);
}
$pipeline->run_pipeline();
