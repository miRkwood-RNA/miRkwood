#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::ResultsExporterMaker');
}
require_ok('miRkwood::ResultsExporterMaker');


## Set up ##

# Necessary as the headers are fetched in the config
use miRkwood;
use miRkwood::Results;
my $config_file = input_file('run_options.cfg');
miRkwood->CONFIG_FILE($config_file);

my $candidates_YML_file = input_file('basic_candidates.yml');
ok(
    my %results = miRkwood::Results->deserialize_results($candidates_YML_file),
    'Can call deserialize_results'
);


## CSV ##

ok( my $csv_exporter = miRkwood::ResultsExporterMaker->make_csv_results_exporter( 'abinitio' ),
    'can call make_csv_results_exporter');
isa_ok($csv_exporter, 'miRkwood::ResultsExporter::CSVExporter');
ok($csv_exporter->initialize('', \%results, ['1-1']),
    'can initialize CSVExporter');
ok(
    my $output_csv = $csv_exporter->perform_export(),
    'can export with CSVExporter'
);
my $expected_csv = slurp_file( input_file('ResultsExporter.csv.output') );
is( $output_csv, $expected_csv, 'CSVExporter returns the correct value' );


## FASTA ##

ok( my $fasta_exporter = miRkwood::ResultsExporterMaker->make_fasta_results_exporter(),
    'can call make_fasta_results_exporter');
isa_ok($fasta_exporter, 'miRkwood::ResultsExporter::FastaExporter');
ok( $fasta_exporter->initialize('', \%results, ['1-1']),
    'can initialize FastaExporter');
ok(
    my $output_fasta = $fasta_exporter->perform_export(),
    'can export with FastaExporter'
);

my $expected_fasta =
  slurp_file( input_file( 'workspace', '1', '1', 'seq.txt' ) );

is( $output_fasta, $expected_fasta,
    'FastaExporter returns the correct value' );


## Dot-bracket ##

ok( my $dotbracket_exporter = miRkwood::ResultsExporterMaker->make_dotbracket_results_exporter(),
    'can call make_dotbracket_results_exporter');
isa_ok($dotbracket_exporter, 'miRkwood::ResultsExporter::DotBracketExporter');
ok( $dotbracket_exporter->initialize('', \%results, ['1-1']),
    'can initialize DotBracketExporter');
ok(
    my $output_dotbracket = $dotbracket_exporter->perform_export(),
    'can export with DotBracketExporter'
);

my $expected_dotbracket =
  slurp_file( input_file('ResultsExporter.dotbracket.output') );

is( $output_dotbracket, $expected_dotbracket,
    'DotBracketExporter returns the correct value' );


## GFF ##

ok( my $gff_exporter = miRkwood::ResultsExporterMaker->make_gff_results_exporter(),
    'can call make_gff_results_exporter');
isa_ok($gff_exporter, 'miRkwood::ResultsExporter::GFFExporter');
ok( $gff_exporter->initialize('', \%results, ['1-1']),
    'can initialize GFFExporter');
ok(
    my $output_gff = $gff_exporter->perform_export(),
    'can export with GFFExporter'
);

my $expected_gff =
  slurp_file( input_file('ResultsExporter.gff.output') );

is( $output_gff, $expected_gff,
    'GFFExporter returns the correct value' );


## HTML ##

ok( my $html_exporter = miRkwood::ResultsExporterMaker->make_html_results_exporter( 'abinitio' ),
    'can call make_html_results_exporter');
isa_ok($html_exporter, 'miRkwood::ResultsExporter::HTMLExporter');
ok( $html_exporter->initialize('', \%results, ['1-1']),
    'can initialize HTMLExporter');
ok(
    my $output_html = $html_exporter->perform_export(),
    'can export with HTMLExporter'
);
my $expected_html =
  slurp_file( input_file('ResultsExporter.html.output') );

is( $output_html, $expected_html,
    'HTMLExporter returns the correct value' );
    
