#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use File::Spec;
use feature 'switch';
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::Results;
use miRkwood::WebTemplate;

my $cgi = CGI->new;

my $id_job = $cgi->param('run_id');    # récupération id job
my $export_type = $cgi->param('type');
my $data= $cgi->param('data');  # identifiant des candidats sélectionnés
my @sequences_to_export =  split( ',',$data  );

my $valid = miRkwood::Results->is_valid_jobID($id_job);

if ($valid) {

    given ($export_type) {
        when (/gff/) { exportAsGFF($id_job) }
        when (/odf/) { exportAsODF($id_job) }
        when (/csv/) { exportAsCSV($id_job) }
        when (/fas/) { exportAsFasta($id_job) }
        when (/dot/) { exportAsDotBracket($id_job) }
        default { miRkwood::WebTemplate::web_die("The export type '$export_type' is not supported"); }
    }
}
else {
    miRkwood::WebTemplate::web_die('The job ID provided is not valid');
}

sub exportAsGFF {
    my $id_job = shift @_;
    my %myResults =  miRkwood::Results->get_structure_for_jobID($id_job);
    my $gff = miRkwood::Results->export('gff', \%myResults , \@sequences_to_export);
    print <<"DATA" or miRkwood::WebTemplate::web_die("Error when printing content: $!");
Content-type: text/gff
Content-disposition: attachment;filename=Results-$id_job.gff

$gff
DATA
}

sub exportAsODF {
    my $id_job = shift @_;
    require miRkwood::OpenDocument
        or miRkwood::WebTemplate::web_die("Fail to import OpenDocument. Perl module ODF::lpOD is likely missing.");
    my $odf_factory = miRkwood::OpenDocument->new($id_job, \@sequences_to_export);
    my $odt = $odf_factory->get_report();
    print $cgi->redirect( -uri => $odt );
}

sub exportAsCSV {
    my $id_job = shift @_;
    my %myResults =  miRkwood::Results->get_structure_for_jobID($id_job);
    my $csv = miRkwood::Results->resultstruct2csv( \%myResults , \@sequences_to_export);

    print <<"DATA" or miRkwood::WebTemplate::web_die("Error when printing content: $!");
Content-type: text/csv
Content-disposition: attachment;filename=Results-$id_job.csv

$csv
DATA
}

sub exportAsFasta {
    my $id_job = shift @_;
    my %myResults =  miRkwood::Results->get_structure_for_jobID($id_job);
    my $fasta = miRkwood::Results->export('fas', \%myResults , \@sequences_to_export);
    print <<"DATA" or miRkwood::WebTemplate::web_die("Error when printing content: $!");
Content-type: text/txt
Content-disposition: attachment;filename=Results-$id_job.fa

$fasta
DATA
}

sub exportAsDotBracket {
    my $id_job = shift @_;
    my %myResults =  miRkwood::Results->get_structure_for_jobID($id_job);
    my $dotbracket = miRkwood::Results->export('dot', \%myResults , \@sequences_to_export);
    print <<"DATA" or miRkwood::WebTemplate::web_die("Error when printing content: $!");
Content-type: text/txt
Content-disposition: attachment;filename=Results-$id_job.txt

$dotbracket
DATA
}
