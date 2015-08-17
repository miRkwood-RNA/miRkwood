package miRkwood::ResultsExporter::HTMLExporterSmallRNAseqKnown;

# ABSTRACT: Class for exporting results in HTML table

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::HTMLExporter';

use miRkwood::Candidate;
use miRkwood::Utils;

sub get_headers {
    my ( $self, @args ) = @_;
    my @optional_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @headers = qw{identifier mirbase_id position strand mirna_sequence mirna_length quality nb_reads};
    return @headers;
}


sub get_header {
    my ( $self, @args ) = @_;
    my $output = '<tr>';
    for my $header ( ('chromosome'), $self->get_headers() ) {
        if ( $header eq 'mirbase_id' ){
            $output .= "<th>miRBase identifier</th>\n";
        }
        elsif ( $header eq 'mirna_sequence' ){
            $output .= "<th>miRNA</th>\n";
        }
        elsif ( $header eq 'nb_reads' ){
            $output .= "<th>reads</th>\n";
        }
        elsif ( $header eq 'mirna_length' ){
            $output .= "<th>miRNA length</th>\n";
        }
        else {
            $output .= "<th>$header</th>\n";
        }
    }
    $output .= "</tr>\n";
    return $output;
}


1;
