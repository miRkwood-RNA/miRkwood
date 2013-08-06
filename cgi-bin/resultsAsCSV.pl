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

my $local_dir = dirname( abs_path($0) );
my $rootdir = File::Spec->catdir( $local_dir, '..' );
my $chaine = "";
my $id_job = $cgi->param('run_id');    # récupération id job
my $data= $cgi->param('data');  # identifiant des candidats sélectionnés
my @tab =  split( ',',$data  );
my $valid = PipelineMiRNA::WebFunctions->is_valid_jobID($id_job);

if ($valid) {
    my %myResults =
      PipelineMiRNA::WebFunctions->get_structure_for_jobID($id_job);
    my $HTML_results =
      PipelineMiRNA::WebFunctions->resultstruct2pseudoXML( \%myResults );
    my $csv = PipelineMiRNA::WebFunctions->resultstruct2csv( \%myResults , \@tab);

    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/csv
Content-disposition: attachment;filename=Results-$id_job.csv

$csv
DATA
}
else {
    die('Error');
}
