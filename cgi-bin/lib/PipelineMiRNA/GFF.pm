package PipelineMiRNA::GFF;

# ABSTRACT: Exporting pipeline results as GFF

use strict;
use warnings;

# Export to Generic Feature Format (GFF)

use PipelineMiRNA::Paths;
use File::Spec;
use PipelineMiRNA::Results;
use Data::Dumper;

=method generate_GFF

Generate a GFF file from the results structure.

Usage:
my $gff = PipelineMiRNA::GFF->generate_GFF(\%results);

=cut

sub generate_GFF {
    my ( $self, @args ) = @_;
    my $results = shift @args;
    my @sequences_to_export = shift @args;

    my %results = %{$results};

    my $output = '##gff-version 3';

    my @keys = sort keys %results;
    foreach my $key(@keys)
    {
        if (  $key ~~ \@sequences_to_export )
        {
            my $value = $results{$key};
            $output .= $self->make_gff_line_for_candidate($value);
        }
    }
    return $output;
}

=method make_gff_line_for_candidate

Generate a GFF line from a candidate structure.

Usage:
my $gff_line = PipelineMiRNA::GFF->make_gff_line_for_candidate($value);(\%candidate);

=cut

sub make_gff_line_for_candidate {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my ( $start, $end ) = split( m/[-]/xms, $candidate{'position'} );
    my $text .= "\n" .              # BEGIN
      $candidate{'name'} . "\t" .   # seqid
      '.' . "\t" .                  # source
      'miRNA' . "\t" .              # type
      $start . "\t" .               # start
      $end . "\t" .                 # end
      '.' . "\t" .                  # score
      '.' . "\t" .                  # strand
      '.' . "\t" .                  # phase
      '.' . "\t" .                  # attributes
      q{};
    return $text;
}

=method generate_GFF_from_ID

Generate the GFF for a given jobID

Usage:
my $gff = PipelineMiRNA::GFF->generate_GFF($id_job);

=cut

sub generate_GFF_from_ID {
    my ( $self, @args ) = @_;
    my $jobId   = shift @args;
    my @sequences_to_export = shift @args;
    my %results = PipelineMiRNA::Results->get_structure_for_jobID($jobId);
    return $self->generate_GFF( \%results, \@sequences_to_export);
}

1;
