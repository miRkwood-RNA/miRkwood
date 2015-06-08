package miRkwood::ResultsExporter::HTMLExporterGenomic;

# ABSTRACT: Class for exporting results in HTML table

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::HTMLExporter';

use miRkwood::Candidate;
use miRkwood::Utils;

sub get_headers {
    my ( $self, @args ) = @_;
    my @optional_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @headers =
      ( 'identifier', 'position', 'length', 'strand', 'quality', 'mfe', 'mfei', 'amfe', @optional_fields );
    return @headers;
}


1;
