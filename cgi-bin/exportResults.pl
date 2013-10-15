#!/usr/bin/perl -w
use strict;
use warnings;

use File::Basename qw(dirname);
use File::Spec;
use CGI::Carp qw(fatalsToBrowser);
use CGI;
use Switch;

my $cgi = CGI->new;
use Cwd qw( abs_path );
use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory
use PipelineMiRNA::WebFunctions;


my $id_job = $cgi->param('run_id');    # récupération id job
my $export_type = $cgi->param('type');
my $data= $cgi->param('data');  # identifiant des candidats sélectionnés
my @sequences_to_export =  split( ',',$data  );

my $valid = PipelineMiRNA::WebFunctions->is_valid_jobID($id_job);

if ($valid) {
    switch ($export_type) {
        case 'gff'        { exportAsGFF() }
        case 'odf'        { exportAsODF() }
        case 'csv'        { exportAsCSV() }
        case 'fasta'      { }
        case 'dotbracket' { }
        else              { die('Error: the export type is not supported'); }
    }
}
else {
    die('Error: the job ID provided is not valid');
}

sub exportAsGFF {
    use PipelineMiRNA::GFF;
    my $gff = PipelineMiRNA::GFF->generate_GFF_from_ID($id_job, \@sequences_to_export);
    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/gff
Content-disposition: attachment;filename=Results-$id_job.gff

$gff
DATA
}

sub exportAsODF {
    #TODO: Filter selected sequences.
    use PipelineMiRNA::OpenDocument;
    my $odt = PipelineMiRNA::OpenDocument->generate_report($id_job);
    print $cgi->redirect( -uri => $odt );
}

sub exportAsCSV {
    my %myResults =  PipelineMiRNA::WebFunctions->get_structure_for_jobID($id_job);
    my $csv = PipelineMiRNA::WebFunctions->resultstruct2csv( \%myResults , \@sequences_to_export);

    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/csv
Content-disposition: attachment;filename=Results-$id_job.csv

$csv
DATA
}
