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

sub get_csv_headers {
    my ( $self, @args ) = @_;
    my @optional_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @csv_headers     = (
        'name', 'identifier', 'start_position', 'end_position', 'quality', '%GC',
        @optional_fields, 'structure_stemloop', 'sequence'
    );
    return @csv_headers;
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
        if ($header eq "reads"){
            $contents = 0;
            foreach my $key (keys( %{$candidate->{'reads'}} )){
                $contents += $candidate->{'reads'}{$key};
            }
        }        
        if ( !defined $contents ) {
            $contents = q{};
        }
        $result .= "$contents,";
    }
    $result .= "\n";
    return $result;
}

1;
