package miRkwood::ResultsExporter::CSVExporterSmallRNAseqKnown;

# ABSTRACT: Class for exporting results in CSV

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::CSVExporter';


sub get_csv_headers {
    my ( $self, @args ) = @_;
    my @csv_headers = (qw{name identifier mirbase_id start_position end_position strand},
        qw{mirna_sequence mirna_length quality %GC nb_reads structure_stemloop sequence});
    return @csv_headers;
}

1;
