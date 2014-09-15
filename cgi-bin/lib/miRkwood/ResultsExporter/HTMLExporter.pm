package miRkwood::ResultsExporter::HTMLExporter;

# ABSTRACT: Class for exporting results in HTML table

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

use miRkwood::Candidate;

sub get_headers {
    my ( $self, @args ) = @_;
    my @optional_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @headers =
      ( 'position', 'length', 'strand', 'quality', @optional_fields );
    return @headers;
}

sub get_header {
    my ( $self, @args ) = @_;
    my $output .= "<tr>";
    for my $header ( ('name'), $self->get_headers() ) {

        $output .= "<th>$header</th>\n";
    }
    $output .= "</tr>\n";
    return $output;
}

sub export_candidate {
    my ( $self, @args ) = @_;
    my $candidate = shift @args;
    my $output .= '<tr>';
    my $anchor   = "${$candidate}{'name'}-${$candidate}{'position'}";
    my $contents = "<a href='#$anchor'>${$candidate}{'name'}</a>";
    $output .= "<td>$contents</td>\n";
    for my $header ($self->get_headers()) {
        my $td_content = "";
        my $contents   = ${$candidate}{$header};
        if ( !defined $contents ) {
            $contents = q{};
        }
        $output .= "<td>$contents</td>\n";
    }
    $output .= "\n</tr>\n";
}

sub perform_export{
    my ( $self, @args ) = @_;

    my %results = %{$self->{'results'}};
    my $output = '';

    $output .= "<table>\n<tbody>";

    $output .= $self->get_header();

    my @keys = $self->get_sorted_keys();
    foreach my $key (@keys) {
        if ( $self->is_sequence_to_export($key)){
            my $candidate = $results{$key};
            $output .= $self->export_candidate($candidate);
        }
    }
    $output .= "</tbody>\n</table>\n";
    return $output;
}

1;