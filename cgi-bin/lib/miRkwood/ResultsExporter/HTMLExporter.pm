package miRkwood::ResultsExporter::HTMLExporter;

# ABSTRACT: Class for exporting results in HTML table

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

use miRkwood::Candidate;
use miRkwood::Utils;

sub get_headers {
    my ( $self, @args ) = @_;
    my @optional_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @headers =
      ( 'identifier', 'position', 'length', 'strand', 'quality', @optional_fields );
    return @headers;
}

sub get_header {
    my ( $self, @args ) = @_;
    my $type = shift @args;
    my @headers = $self->get_headers();
    if ( $type eq 'Known' ){
        @headers = miRkwood::Utils::delete_element_in_array( 'alignment', \@headers );
        @headers = miRkwood::Utils::delete_element_in_array( 'shuffles', \@headers );
        @headers = miRkwood::Utils::delete_element_in_array( 'mfe', \@headers );
        @headers = miRkwood::Utils::delete_element_in_array( 'mfei', \@headers );
        @headers = miRkwood::Utils::delete_element_in_array( 'amfe', \@headers );
        push @headers, 'precursor_name';
    }    
    my $output .= '<tr>';
    for my $header ( ('name'), @headers ) {

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
    my @headers = $self->get_headers();
    if ( defined( $candidate->{'precursor_name'} ) ){
        @headers = miRkwood::Utils::delete_element_in_array( 'alignment', \@headers );
        @headers = miRkwood::Utils::delete_element_in_array( 'shuffles', \@headers );
        push @headers, 'precursor_name';
    }
    for my $header ( @headers ) {
        my $td_content = '';
        my $contents   = ${$candidate}{$header};
        if ($header eq 'reads'){
            $contents = 0;
            foreach my $key (keys( %{$candidate->{'reads'}} )){
                $contents += $candidate->{'reads'}{$key};
            }
        }
        elsif ($header eq 'quality'){
            $contents = '<center><font color="#FF8000">';
            for (my $i = 0; $i < ${$candidate}{"quality"}; $i++){
                $contents .= "*";
            }
            $contents .= '</font></center>';
        }
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
    my @keys = $self->get_sorted_keys();
    my $nb_results = scalar( @keys );

    my $output = '';
    $output .= "<h3>$nb_results candidates found.</h3>\n";

    $output .= "<table>\n<tbody>";

    my $type = 'New';
    if ( defined( $results{$keys[0]}->{'precursor_name'} ) ){
        $type = 'Known';
    }
    $output .= $self->get_header( $type );

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
