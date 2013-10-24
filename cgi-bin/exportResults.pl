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
use PipelineMiRNA::Results;


my $id_job = $cgi->param('run_id');    # récupération id job
my $export_type = $cgi->param('type');
my $data= $cgi->param('data');  # identifiant des candidats sélectionnés
my @sequences_to_export =  split( ',',$data  );

my $valid = PipelineMiRNA::Results->is_valid_jobID($id_job);

if ($valid) {
    switch ($export_type) {
        case 'gff'        { exportAsGFF($id_job) }
        case 'odf'        { exportAsODF($id_job) }
        case 'csv'        { exportAsCSV($id_job) }
        case 'fas'        { exportAsFasta($id_job) }
        case 'dot'        { exportAsDotBracket($id_job) }
        else              { die("Error: the export type '$export_type' is not supported"); }
    }
}
else {
    die('Error: the job ID provided is not valid');
}

sub exportAsGFF {
    my $id_job = shift @_;
    use PipelineMiRNA::GFF;
    my $gff = PipelineMiRNA::GFF->generate_GFF_from_ID($id_job, \@sequences_to_export);
    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/gff
Content-disposition: attachment;filename=Results-$id_job.gff

$gff
DATA
}

sub exportAsODF {
    my $id_job = shift @_;
    #TODO: Filter selected sequences.
    use PipelineMiRNA::OpenDocument;
    my $odt = PipelineMiRNA::OpenDocument->generate_report($id_job);
    print $cgi->redirect( -uri => $odt );
}

sub exportAsCSV {
    my $id_job = shift @_;
    my %myResults =  PipelineMiRNA::Results->get_structure_for_jobID($id_job);
    my $csv = PipelineMiRNA::Results->resultstruct2csv( \%myResults , \@sequences_to_export);

    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/csv
Content-disposition: attachment;filename=Results-$id_job.csv

$csv
DATA
}

sub exportAsFasta {
    my $id_job = shift @_;
    my %myResults =  PipelineMiRNA::Results->get_structure_for_jobID($id_job);
    my $fasta = PipelineMiRNA::Results->export('fas', \%myResults , \@sequences_to_export);
    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/txt
Content-disposition: attachment;filename=Results-$id_job.txt

$fasta
DATA
}

sub exportAsDotBracket {
    my $id_job = shift @_;
    my %myResults =  PipelineMiRNA::Results->get_structure_for_jobID($id_job);
    my $dotbracket = PipelineMiRNA::Results->export('dot', \%myResults , \@sequences_to_export);
    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/txt
Content-disposition: attachment;filename=Results-$id_job.txt

$dotbracket
DATA
}
