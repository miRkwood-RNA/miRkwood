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
my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($id_job);
my $path_to_reads_clouds_archive = miRkwood::Results->create_reads_archive( $absolute_job_dir );



my ( $filename, $contents, $disposition );

( $filename, $contents, $disposition ) = exportFile( $path_to_reads_clouds_archive );

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
    my $contents = miRkwood::FileUtils::slurp_bin_file ( $file );
    return ( $filename, $contents, 'attachment' );
}
