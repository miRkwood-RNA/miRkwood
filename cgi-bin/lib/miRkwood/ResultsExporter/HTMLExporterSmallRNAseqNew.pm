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
    my $left_border_cell = "class='left_border'";
    my @optional_candidate_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @optional_mirna_fields = miRkwood::Candidate->get_optional_mirna_fields();
    my $nb_mirna_fields = scalar( @optional_mirna_fields ) + 3;

    my $output = '<tr>';

    for my $header ( ('chromosome'), $self->get_headers() ) {
        # candidate fields
        if ( $header eq 'chromosome' ){
            $output .= "<th rowspan='2'>chr</th>\n";
        }
        elsif ( $header eq 'position' ){
            $output .= "<th rowspan='2'>position</th>\n";
        }
        elsif ( $header eq 'strand' ){
            $output .= "<th rowspan='2'>strand</th>\n";
        }
        elsif ( $header eq 'nb_reads' ){
            $output .= "<th rowspan='2'>reads</th>\n";
        }
        elsif ( $header eq 'reads_distribution' ){
            $output .= "<th rowspan='2'>reads distribution</th>\n";
        }
        elsif ( $header eq 'mfei' ){
            $output .= "<th rowspan='2'>&nbsp; mfei &nbsp;</th>\n";
        }
        elsif ( $header eq 'shuffles' ){
            $output .= "<th rowspan='2'>shuffles</th>\n";
        }
        # miRNA fields
        elsif ( $header eq 'mirna_sequence' ){
            $output .= "<th colspan='$nb_mirna_fields' $left_border_cell>miRNA</th>\n";
            $output .= "</tr>\n";
            $output .= "<tr>\n";
            $output .= "<th class='left_border'>sequence</th>\n";
        }
        elsif ( $header eq 'mirna_length' ){
            $output .= "<th>length</th>\n";
        }
        else{
            $output .= "<th>$header</th>\n";
        }

    }

    $output .= "</tr>\n";
    return $output;
}


1;
