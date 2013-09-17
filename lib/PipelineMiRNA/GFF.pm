package PipelineMiRNA::GFF;

# ABSTRACT: Exporting pipeline results as GFF

use strict;
use warnings;

# Export to Generic Feature Format (GFF)

use PipelineMiRNA::Paths;
use File::Spec;
use PipelineMiRNA::WebFunctions;
use Data::Dumper;

=method generate_GFF

Generate a GFF file from the results structure.

Usage:
my $gff = PipelineMiRNA::GFF->generate_GFF(\%results);

=cut

sub generate_GFF {
    my ( $self, @args ) = @_;
    my $results = shift @args;
    my %results = %{$results};

    my $output = '##gff-version 3';

    my @keys = sort keys %results;

    foreach my $key (@keys) {
        my $value = $results{$key};
        my ( $start, $end ) = split( m/[-]/xms, ${$value}{'position'} );
        $output .= "\n" .    # BEGIN
          ${$value}{'name'} . "\t" .    # seqid
          '.' . "\t" .                  # source
          'miRNA' . "\t" .              # type
          $start . "\t" .               # start
          $end . "\t" .                 # end
          '.' . "\t" .                  # score
          '.' . "\t" .                  # strand
          '.' . "\t" .                  # phase
          '.' . "\t" .                  # attributes
          q{};
    }
    return $output;
}

=method generate_GFF_from_ID

Generate the GFF for a given jobID

Usage:
my $gff = PipelineMiRNA::GFF->generate_GFF($id_job);

=cut

sub generate_GFF_from_ID {
    my ( $self, @args ) = @_;
    my $jobId   = shift @args;
    my %results = PipelineMiRNA::WebFunctions->get_structure_for_jobID($jobId);
    return $self->generate_GFF( \%results );
}

1;
