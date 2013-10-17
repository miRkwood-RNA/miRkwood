#!/usr/bin/perl -w

use warnings;
use strict;

use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory
use PipelineMiRNA::MainPipeline;


## Code ##
my ( $ifilterCDS, $imfei, $irandfold,  $ialign, $idirJob, $iplant ) = @ARGV;

PipelineMiRNA::MainPipeline::main_entry( $ifilterCDS, $imfei, $irandfold, $ialign, $idirJob, $iplant );
