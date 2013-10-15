#!/usr/bin/perl -w

use warnings;
use strict;

use CGI::Carp qw(fatalsToBrowser);

use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory
use PipelineMiRNA::MainPipeline;


## Code ##
my ( $icheck, $imfei, $irandfold,  $ialign, $idirJob, $iplant ) = @ARGV;

PipelineMiRNA::MainPipeline::main_entry( $icheck, $imfei, $irandfold, $ialign, $idirJob, $iplant );
