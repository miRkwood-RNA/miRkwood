#!/usr/bin/perl -w
use strict;
use warnings;

use Class::Struct;
use File::Basename qw(dirname);
use File::Spec;
use CGI::Carp qw(fatalsToBrowser);
use CGI;
my $cgi = CGI->new;
use Cwd qw( abs_path );
use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory
use PipelineMiRNA::WebFunctions;

my $local_dir = dirname( abs_path($0) );
my $rootdir = File::Spec->catdir( $local_dir, '..' );

my $id_job      = $cgi->param('run_id');    # récupération id job
my $dirJob_name = 'job' . $id_job;
my $dirJob =
  abs_path( File::Spec->catdir( $rootdir, 'results', $dirJob_name ) );

my $csv_file = File::Spec->catfile( $dirJob, 'result.csv' );
open( my $CSV, '<', $csv_file ) or die "Error when opening $csv_file: $!";
my @csv = <$CSV>;
close $CSV or die "Cannot close $csv_file: $!";
print <<"DATA" or die "Error when printing content: $!";
Content-type: text/csv
Content-disposition: attachment;filename=Results.csv

@csv
DATA
