#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
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
use PipelineMiRNA::OpenDocument;

my $local_dir = dirname( abs_path($0) );
my $rootdir = File::Spec->catdir( $local_dir, '..' );
my $chaine = "";
my $id_job = $cgi->param('run_id');    # récupération id job

my $valid = PipelineMiRNA::WebFunctions->is_valid_jobID($id_job);

if ($valid) {
    my $odt = PipelineMiRNA::OpenDocument->generate_report($id_job);
    print $cgi->redirect( -uri => $odt );
}
else {
    die('Error');
}
