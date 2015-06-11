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
use miRkwood::ResultsExporterMaker;

my $cgi = CGI->new;

my $id_job = $cgi->param('run_id');    # récupération id job
my $export_type = $cgi->param('type');
my $data= $cgi->param('data');  # identifiant des candidats sélectionnés
my $header_type = $cgi->param('header_type');
my $mirna_type = $cgi->param('mirna_type');
my @sequences_to_export =  split( ',',$data  );

my $valid = miRkwood::Results->is_valid_jobID($id_job);

if (! $valid) {
    miRkwood::WebTemplate::web_die('The job ID provided is not valid');
}

my $exporter;
my %myResults =  miRkwood::Results->get_structure_for_jobID($id_job, $mirna_type);

given ($export_type) {
    when (/gff/) {
        $exporter = miRkwood::ResultsExporterMaker->make_gff_results_exporter();
    }
    when (/odf/) {
        $exporter = miRkwood::ResultsExporterMaker->make_opendocument_results_exporter();
        $exporter->{'cgi'} = $cgi;
    }
    when (/csv/) {
        $exporter = miRkwood::ResultsExporterMaker->make_csv_results_exporter($header_type . $mirna_type);
    }
    when (/fas/) {
        $exporter = miRkwood::ResultsExporterMaker->make_fasta_results_exporter();
    }
    when (/dot/) {
        $exporter = miRkwood::ResultsExporterMaker->make_dotbracket_results_exporter();
    }
    when (/reads/) {
        $exporter = miRkwood::ResultsExporterMaker->make_reads_clouds_results_exporter( $mirna_type );
    }    
    default { miRkwood::WebTemplate::web_die("The export type '$export_type' is not supported"); }
}
$exporter->initialize($id_job, \%myResults, \@sequences_to_export);
print $exporter->export_for_web();
