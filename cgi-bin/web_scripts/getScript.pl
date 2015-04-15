#!/usr/bin/perl
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::WebTemplate;
use miRkwood::FileUtils;

my %allowed_files = ( 'mirkwood-bam2bed.pl' => 1 );

my $cgi = CGI->new();

my $file = $cgi->param('file');

my $script_file = '';

if ( $allowed_files{ $file } ){
    $script_file = File::Spec->catfile( miRkwood::Paths::get_scripts_path(), $file );
}

my ( $filename, $contents, $disposition );

( $filename, $contents, $disposition ) = exportFile( $script_file );

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

