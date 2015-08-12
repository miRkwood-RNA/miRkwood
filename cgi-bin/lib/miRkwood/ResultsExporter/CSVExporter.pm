package miRkwood::ResultsExporter::CSVExporter;

# ABSTRACT: Class for exporting results in CSV

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

sub get_content_type {
    my ( $self, @args ) = @_;
    return 'text/csv';
}

sub get_file_extension {
    my ( $self, @args ) = @_;
    return 'csv';
}

sub get_header {
    my ( $self, @args ) = @_;
    return join( ',', $self->get_csv_headers() ) . "\n";
}

sub export_candidate {
    my ( $self, @args ) = @_;
    my $candidate = shift @args;
    my $result;
    for my $header ($self->get_csv_headers()) {
        my $contents = ${$candidate}{$header};
        if ( !defined $contents ) {
            $contents = q{};
        }
        $result .= "$contents,";
    }
    $result .= "\n";
    return $result;
}

1;
