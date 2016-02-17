package miRkwood::ResultsExporter::OrgExporter;

# ABSTRACT: Class for exporting results in Org-mode

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

sub get_content_type {
    my ( $self, @args ) = @_;
    return 'text/txt';
}

sub get_file_extension {
    my ( $self, @args ) = @_;
    return 'org';
}

sub export_candidate {
    my ( $self, @args ) = @_;
    my $candidate = shift @args;
    return $candidate->candidateAsOrg();
}

1;
