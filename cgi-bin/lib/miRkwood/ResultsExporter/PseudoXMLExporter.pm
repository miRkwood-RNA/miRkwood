package miRkwood::ResultsExporter::PseudoXMLExporter;

# ABSTRACT: Class for exporting results in pseudo-XML

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

sub export_candidate {
    my ( $self, @args ) = @_;
    my $candidate = shift @args;
    return $candidate->candidate_as_pseudoXML() . "\n";
}

sub perform_export{
    my ( $self, @args ) = @_;

    my %results = %{$self->{'results'}};
    my $output = '';

    $output .= "<results id='all'>\n";
    my @keys = $self->get_sorted_keys();
    foreach my $key (@keys) {
        if ( $self->is_sequence_to_export($key)){
            my $candidate = $results{$key};
            $output .= $self->export_candidate($candidate);
        }
    }
    $output .= "</results>\n";
    $output .= "<results id='all2'>\n";
    @keys = sort {
        ( $results{$b}->{'quality'} cmp $results{$a}->{'quality'} )
          || (
            $results{$a}->{'start_position'} <=> $results{$b}->{'start_position'} )
    } keys %results;
    foreach my $key (@keys) {
        my $candidate = $results{$key};
        $output .= $candidate->candidate_as_pseudoXML() . "\n";
    }
    $output .= "</results>";

    return $output;
}

1;
