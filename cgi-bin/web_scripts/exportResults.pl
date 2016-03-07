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
my $pipeline_type = $cgi->param('pipeline');
my $mirna_type = $cgi->param('mirna_type');
my @sequences_to_export =  split( /,/,$data  );

my $valid = miRkwood::Results->is_valid_jobID($id_job);

if (! $valid) {
    miRkwood::WebTemplate::web_die('The job ID provided is not valid');
}

my $exporter;
my %myResults =  miRkwood::Results->get_structure_for_jobID($id_job, $mirna_type);

given ($export_type) {
    when (/gff/) {
        $exporter = miRkwood::ResultsExporterMaker->make_gff_results_exporter( $mirna_type );
    }
    when (/odf/) {
        $exporter = miRkwood::ResultsExporterMaker->make_opendocument_results_exporter( $mirna_type );
        $exporter->{'cgi'} = $cgi;
    }
    when (/org/) {
        $exporter = miRkwood::ResultsExporterMaker->make_org_results_exporter( $mirna_type );
    }
    when (/csv/) {
        $exporter = miRkwood::ResultsExporterMaker->make_csv_results_exporter($pipeline_type, $mirna_type );
    }
    when (/fas/) {
        $exporter = miRkwood::ResultsExporterMaker->make_fasta_results_exporter( $mirna_type );
    }
    when (/dot/) {
        $exporter = miRkwood::ResultsExporterMaker->make_dotbracket_results_exporter( $mirna_type );
    }
    when (/reads/) {
        $exporter = miRkwood::ResultsExporterMaker->make_reads_clouds_results_exporter( $mirna_type );
    }
    when (/pdf/) {
        # Create ORG file, mandatory for PDF creation
        my $dir_for_org_file = '/tmp/';
        my $exporter_org = miRkwood::ResultsExporterMaker->make_org_results_exporter( $mirna_type );
        $exporter_org->initialize($id_job, \%myResults, \@sequences_to_export);
        $exporter_org->export_on_disk( $dir_for_org_file );
        my $org_file = File::Spec->catfile($dir_for_org_file, $exporter_org->get_filename());

        # Call exporter for PDF
        $exporter = miRkwood::ResultsExporterMaker->make_pdf_results_exporter( $mirna_type, $org_file );
    }

    default { miRkwood::WebTemplate::web_die("The export type '$export_type' is not supported"); }
}
$exporter->initialize($id_job, \%myResults, \@sequences_to_export);
print $exporter->export_for_web();
