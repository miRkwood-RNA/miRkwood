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
Filtering on the provided array of name
Usage:
my $gff = PipelineMiRNA::GFF->generate_GFF(\%results, \@seqs_to_export);

=cut

sub generate_GFF {
    my ( $self, @args ) = @_;
    my $results_ref = shift @args;
    my $sequences_to_export_ref = shift @args;

    my %results = %{$results_ref};

    my @sequences_to_export;
    if (! eval { @sequences_to_export = @{$sequences_to_export_ref} } ){
        @sequences_to_export = ();
    }
    my $no_seq_selected = ! (scalar @sequences_to_export);

    my $output = '##gff-version 3' . "\n";

    my @keys = sort keys %results;
    foreach my $key(@keys)
    {
        if ( ($key ~~ @sequences_to_export) || ( $no_seq_selected ) )
        {
            my $value = $results{$key};
            $output .= PipelineMiRNA::Candidate->candidate_as_gff($value);
        }
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
    my $sequences_to_export_ref = shift @args;
    my %results = PipelineMiRNA::Results->get_structure_for_jobID($jobId);
    return $self->generate_GFF( \%results, $sequences_to_export_ref);
}

1;
