package miRkwood::CandidateHandler;

# ABSTRACT: Code to manipulate around Candidate objects

use strict;
use warnings;

use miRkwood::Candidate;
use YAML::XS;

=method retrieve_candidate_information

Check correctness and get the result for a given candidate

Arguments:
- $job - the job directory
- $id - the candidate identifier

Returns:
 A miRkwood::Candidate instance
=cut

sub retrieve_candidate_information {
    my ( $self, @args ) = @_;
    my $job = shift @args;
    my $id = shift @args;
    my $candidate_filepath = $self->get_candidate_filepath($job, $id);
    if ( !-e $candidate_filepath ){
        die("Unvalid candidate information in $candidate_filepath ");
    }else{
        return miRkwood::Candidate->new_from_serialized($candidate_filepath);
    }
}

=method get_candidate_filepath

Get the candidate filepath given its identifier and the job directory

=cut

sub get_candidate_filepath {
    my ( $self, @args ) = @_;
    my $job = shift @args;
    my $id = shift @args;
    return File::Spec->catfile($job, 'candidates', $self->make_candidate_filename($id));
}

=method make_candidate_filename

Return the candidate filename based on the identifier.

=cut

sub make_candidate_filename {
    my ( $self, @args ) = @_;
    my $identifier = shift @args;
    return $identifier . '.yml';
}

=method serialize_candidate_information

Serialize the given candidate on disk

Arguments:
- $serialization_path - the filepath to serialize to
- %candidate - the candidate

=cut

sub serialize_candidate_information {
    my ( $self, @args ) = @_;
    my $serialization_path = shift @args;
    my $candidate_object = shift @args;
    my $candidate_base_filename = $self->make_candidate_filename($candidate_object->{'identifier'});
    my $serialization_file = File::Spec->catfile($serialization_path, $candidate_base_filename);
    my %converted_hash = %{$candidate_object};
    return YAML::XS::DumpFile($serialization_file, %converted_hash);
}



1;