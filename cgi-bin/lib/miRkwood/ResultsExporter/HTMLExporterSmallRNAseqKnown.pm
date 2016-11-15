package miRkwood::ResultsExporter::HTMLExporterSmallRNAseqKnown;

# ABSTRACT: Class for exporting results in HTML table

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::HTMLExporter';

use miRkwood::Candidate;
use miRkwood::Utils;

sub get_headers {
    my ( $self, @args ) = @_;
    my @headers = qw{position strand quality precursor_name mirbase_id nb_reads mirna_sequence mirna_length};
    return @headers;
}


sub get_header {
    my ( $self, @args ) = @_;
    my $left_border_cell = "class='left_border'";
    my $output = '<tr>';
    for my $header ( ('chromosome'), $self->get_headers() ) {
        if ( $header eq 'chromosome' ){
            $output .= "<th rowspan='2'>chr</th>\n";
        }
        elsif ( $header eq 'mirbase_id' ){
            $output .= "<th rowspan='2'>miRBase ID</th>\n";
        }
        elsif ( $header eq 'precursor_name' ){
            $output .= "<th rowspan='2'>miRBase name</th>\n";
        }
        elsif ( $header eq 'mirna_sequence' ){
            $output .= "<th colspan='2' $left_border_cell>miRNA</th>\n";
            $output .= "</tr>\n";
            $output .= "<tr>\n";
            $output .= "<th class='left_border'>sequence</th>\n";
        }
        elsif ( $header eq 'mirna_length' ){
            $output .= "<th>length</th>\n";
        }
        elsif ( $header eq 'nb_reads' ){
            $output .= "<th rowspan='2'>reads</th>\n";
        }
        else {
            $output .= "<th rowspan='2'>$header</th>\n";
        }
    }
    $output .= "</tr>\n";
    return $output;
}


1;
