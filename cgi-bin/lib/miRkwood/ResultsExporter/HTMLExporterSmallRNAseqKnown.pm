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
    my @headers = qw{identifier precursor_name position strand length quality reads};
    return @headers;
}


sub get_header {
    my ( $self, @args ) = @_;
    my $output .= "<tr>";
    for my $header ( ('chromosome'), $self->get_headers() ) {
        if ( $header eq 'precursor_name' ){
            $output .= "<th>miRBase name</th>\n";
        }
        else {
            $output .= "<th>$header</th>\n";
        }
    }
    $output .= "</tr>\n";
    return $output;
}


1;
