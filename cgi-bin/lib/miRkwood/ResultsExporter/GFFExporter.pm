package miRkwood::ResultsExporter::GFFExporter;

# ABSTRACT: Abstract class for exporting results.

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

sub get_content_type {
    my ( $self, @args ) = @_;
    return 'text/gff';
}

sub get_file_extension {
    my ( $self, @args ) = @_;
    return 'gff';
}

sub get_header {
    my ( $self, @args ) = @_;
    my $gff_header = "##gff-version 3
# miRNA precursor sequences found by miRkwood have type 'miRNA_primary_transcript'.
# Note, these sequences do not represent the full primary transcript,
# rather a predicted stem-loop portion that includes the precursor.
";
    return $gff_header . "\n";
}

sub export_candidate {
    my ( $self, @args ) = @_;
    my $candidate = shift @args;
    return $candidate->candidate_as_gff();
}

1;
