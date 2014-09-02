#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::Results');
}
require_ok('miRkwood::Results');

my $candidates_dir = input_file('candidates');

ok(
    my %results = miRkwood::Results->deserialize_results($candidates_dir),
    'Can call deserialize_results'
);

my @keys       = keys %results;
my $identifier = $keys[0];
is( $identifier, '1-1', 'deserialize_results correctly deserialized data' );

ok( my $has_candidates = miRkwood::Results->has_candidates( \%results ),
    'can call has_candidates' );
ok( $has_candidates, 'has_candidates ok' );

# Necessary as the headers are fetched in the config
my $config_file = input_file('run_options.cfg');
miRkwood->CONFIG_FILE($config_file);

ok(
    my $output_csv =
      miRkwood::Results->resultstruct2csv( \%results, ['1-1'] ),
    'can call resultstruct2csv'
);
my $expected_csv = slurp_file( input_file('Results.resultstruct2csv.output') );
is( $output_csv, $expected_csv, 'resultstruct2csv returns the correct value' );

ok(
    my $output_xml =
      miRkwood::Results->resultstruct2pseudoXML( \%results ),
    'can call resultstruct2pseudoXML'
);
my $expected_xml =
  slurp_file( input_file('Results.resultstruct2pseudoXML.output') );
chomp $expected_xml;
is( $output_xml, $expected_xml,
    'resultstruct2pseudoXML returns the correct value' );

ok(
    my $output_export_fasta =
      miRkwood::Results->export( 'fas', \%results, ['1-1'] ),
    'can call export for FASTA'
);
my $expected_export_fasta = slurp_file( input_file( 'workspace', '1', '1', 'seq.txt' ) );

is( $output_export_fasta, $expected_export_fasta,
    'export for FASTA return the correct value' );

ok(
    my $output_export_vienna =
      miRkwood::Results->export( 'dot', \%results, ['1-1'] ),
    'can call export for Vienna'
);
my $expected_export_vienna =
  slurp_file( input_file('Results.export.vienna.output') );

is( $output_export_vienna, $expected_export_vienna,
    'export for Vienna return the correct value' );

ok(
    my $output_export_gff =
      miRkwood::Results->export( 'gff', \%results, ['1-1'] ),
    'can call export for GFF'
);
my $expected_export_gff = slurp_file( input_file('Results.export.gff.output') );

is( $output_export_gff, $expected_export_gff,
    'export for GFF return the correct value' );

ok(
    my $number_of_results_output =
      miRkwood::Results->number_of_results( \%results ),
    'can call number_of_results'
);

is( $number_of_results_output, 1,
    'number_of_results return the correct value' );
