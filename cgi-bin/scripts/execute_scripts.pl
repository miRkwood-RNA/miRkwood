#!/usr/bin/perl -w

use warnings;
use strict;

use FindBin;

BEGIN {
    use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
    use miRkwood;
    use miRkwood::MainPipeline;
}


## Code ##
my ( $idirJob ) = @ARGV;

miRkwood::MainPipeline::fasta_pipeline( $idirJob );
