package miRkwood::ResultsExporter::CSVExporterSmallRNAseqNew;

# ABSTRACT: Class for exporting results in CSV

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::CSVExporter';


sub get_csv_headers {
    my ( $self, @args ) = @_;
    my @optional_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @csv_headers     = (
        'name',
        'start_position',
        'end_position',
        'strand',
        'mirna_sequence',
        'mirna_length',
        'mirna_depth',
        'nb_alignments_for_miRNA',
        'weight',
        'quality',
        'nb_reads',
        'reads_distribution',
        '%GC',
        'mfe',
        'mfei',
        'amfe',
        @optional_fields,
        'structure_stemloop',
        'sequence'
    );
    return @csv_headers;
}

1;
