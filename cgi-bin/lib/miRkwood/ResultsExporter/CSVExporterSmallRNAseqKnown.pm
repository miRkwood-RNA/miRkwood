package miRkwood::ResultsExporter::CSVExporterSmallRNAseqKnown;

# ABSTRACT: Class for exporting results in CSV

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::CSVExporter';


sub get_csv_headers {
    my ( $self, @args ) = @_;
    my @optional_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @csv_headers     = (
        'name', 'identifier', 'precursor_name', 'start_position', 'end_position',
        'quality', '%GC', 'reads', 'structure_stemloop', 'sequence'
    );
    return @csv_headers;
}

1;
