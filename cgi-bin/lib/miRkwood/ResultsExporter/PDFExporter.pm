package miRkwood::ResultsExporter::PDFExporter;

# ABSTRACT: Class for exporting results in PDF

use strict;
use warnings;
use File::Spec;
use miRkwood::FileUtils;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my $mirna_type = shift @args;
    my $org_file = shift @args;
    my $self = {
        'identifier'          => undef,
        'results'             => undef,
        'sequences_to_export' => undef,
        'mirna_type'          => $mirna_type,
        'org_file'            => $org_file
    };
    bless $self, $class;
    return $self;
}

sub get_org_file {
    my ( $self, @args ) = @_;
    return $self->{'org_file'};
}

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
    my $org_file = $self->get_org_file();
    my $cmd = "pandoc -o $pdf_file $org_file";
    system( $cmd );
    return; 
}

sub export_for_web {
    my ( $self, @args ) = @_;
    my $filename = $self->get_filename();
    my $pdf_dir = miRkwood::Results->jobId_to_jobPath( $self->get_identifier() );
    my $pdf_file = File::Spec->catfile( $pdf_dir, $filename );
    my $content_type = $self->get_content_type();
    my $content_disposition = $self->get_content_disposition();

    $self->create_PDF_from_ORG( $pdf_file );
    my $contents = miRkwood::FileUtils::slurp_bin_file ( $pdf_file );
    unlink $self->get_org_file();
    unlink $pdf_file;

    my $answer = <<"DATA"
Content-type: $content_type
Content-disposition: $content_disposition;filename=$filename

$contents
DATA
;
    return $answer;

}

1;
