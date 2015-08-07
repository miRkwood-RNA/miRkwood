#!/usr/bin/perl
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::WebTemplate;
use miRkwood::FileUtils;

my $cgi = CGI->new();

my $id_job = $cgi->param('jobId');
my $type   = $cgi->param('type');

my $jobDir = miRkwood::Results->jobId_to_jobPath($id_job);

my $bed_file = miRkwood::Paths::get_bed_file ( $jobDir, $type );

my ( $filename, $contents, $disposition );

( $filename, $contents, $disposition ) = exportFile( $bed_file );

print <<"DATA" or miRkwood::WebTemplate::web_die("Error when printing content: $!");
Content-type: text/plain
Content-disposition: $disposition;filename=$filename

$contents
DATA


sub exportFile {
    my @args = @_;
    my $file = shift @args;
    my $filename = 'tmp';
    if ( $file =~ /.*\/([^\/]+)/ ){
        $filename = $1;
    }
    my $contents = miRkwood::FileUtils::slurp_file ( $file );
    return ( $filename, $contents, 'attachment' );
}

