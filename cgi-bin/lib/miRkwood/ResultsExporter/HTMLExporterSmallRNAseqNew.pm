package miRkwood::ResultsExporter::HTMLExporterSmallRNAseqNew;

# ABSTRACT: Class for exporting results in HTML table

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::HTMLExporter';

use miRkwood::Candidate;
use miRkwood::Utils;

sub get_headers {
    my ( $self, @args ) = @_;
    my @optional_candidate_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @optional_mirna_fields = miRkwood::Candidate->get_optional_mirna_fields();
    my @headers = ( 
        'position',
        'strand',
        'nb_reads',
        'reads_distribution',
        'mfei',
        @optional_candidate_fields,
        'mirna_sequence',
        'mirna_length',
        'weight',
        @optional_mirna_fields );
    return @headers;
}

sub get_header {
    my ( $self, @args ) = @_;
    my $output = '<tr>';
    my $left_border_cell = "class='left_border'";
    for my $header ( ('chromosome'), $self->get_headers() ) {
        if ( $header eq 'chromosome' ){
            $output .= "<th>chr</th>\n";
        }
        elsif ( $header eq 'mirna_sequence' ){
            $output .= "<th $left_border_cell>miRNA</th>\n";
        }
        elsif ( $header eq 'mirna_length' ){
            $output .= "<th>length</th>\n";
        }
        elsif ( $header eq 'nb_reads' ){
            $output .= "<th>reads</th>\n";
        }
        elsif ( $header eq 'reads_distribution' ){
            $output .= "<th>reads distribution</th>\n";
        }
        elsif ( $header eq 'mfei' ){
            $output .= "<th>&nbsp; mfei &nbsp;</th>\n";
        }
        else {
            $output .= "<th>$header</th>\n";
        }
    }
    $output .= "</tr>\n";
    return $output;
}


1;
