package miRkwood::ResultsExporter::DotBracketExporter;

# ABSTRACT: Class for exporting results in Dot Bracket format

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

sub get_content_type {
    my ( $self, @args ) = @_;
    return 'text/txt';
}

sub get_file_extension {
    my ( $self, @args ) = @_;
    return 'txt';
}

sub export_candidate {
    my ( $self, @args ) = @_;
    my $candidate = shift @args;
    return $candidate->candidateAsVienna();
}

1;