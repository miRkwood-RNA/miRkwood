package miRkwood::ResultsExporter::ReadsExporter;

# ABSTRACT: Class for exporting reads clouds

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

use miRkwood;
use miRkwood::FileUtils;


sub new {
    my ( $class, @args ) = @_;
    my $mirna_type = shift @args;
    my $self = {
        identifier => undef,
        results => undef,
        sequences_to_export => undef,
        mirna_type => $mirna_type,
    };
    bless $self, $class;
    return $self;
}

sub get_content_type {
    my ( $self, @args ) = @_;
    return 'text/plain';
}

sub get_file_extension {
    my ( $self, @args ) = @_;
    return 'tar.gz';
}

sub perform_export {
    my ( $self, @args ) = @_;

    my %results = %{$self->{'results'}};

    my @keys = @{$self->get_sequences_to_export()};
    # if all candidates are selected we sent an empty tab
    # so now we have to get all candidates
    if ( scalar(@keys) == 0 ){
        @keys = $self->get_sorted_keys();
    }
    my $archive = miRkwood::Results->create_reads_archive( $self->{'identifier'}, $self->{'mirna_type'}, \@keys );

    my $contents = miRkwood::FileUtils::slurp_bin_file ( $archive );

    return $contents;

}

1;
