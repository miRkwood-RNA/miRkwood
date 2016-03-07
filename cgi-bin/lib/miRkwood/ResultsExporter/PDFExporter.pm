package miRkwood::ResultsExporter::PDFExporter;

# ABSTRACT: Class for exporting results in PDF

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

sub get_content_type {
    my ( $self, @args ) = @_;
    return 'text/txt';
}

sub get_file_extension {
    my ( $self, @args ) = @_;
    return 'pdf';
}

sub create_PDF_from_ORG {
    my ( $self, @args ) = @_;
    my $pdf_file = shift @args;
    my $org_file = shift @args;
    my $cmd = "pandoc -o $pdf_file $org_file";
    system( $cmd );
    return; 
}

1;
