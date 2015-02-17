package miRkwood::Results;

# ABSTRACT: Code directly tied to the results data structure

use strict;
use warnings;

use feature 'switch';
use Time::gmtime;
use File::Spec;

use miRkwood;
use miRkwood::Paths;
use miRkwood::Candidate;
use miRkwood::CandidateHandler;;
use miRkwood::Utils;

=method make_job_id

Return a jobId (based on the current time)

=cut

sub make_job_id {
	my ( $self, @args ) = @_;
	my $type = shift @args;     # type should be 'BAM' or 'Fasta'
	my $now = gmctime();
	$now =~ s/[: ]//g;
	$now = $type . substr( $now, 3 );
	return $now;
}

=method jobId_to_jobPath

Get the job path from a job identifier

=cut

sub jobId_to_jobPath {
	my ( $self, @args ) = @_;
	my $id_job      = shift @args;
	my $dirJob_name = 'job' . $id_job;
	my $results_dir = miRkwood::Paths->get_results_filesystem_path();
	my $jobPath     = File::Spec->catdir( $results_dir, $dirJob_name );
	return $jobPath;
}

=method is_job_finished

Return whether a job is finished or not

=cut

sub is_job_finished {
    my ( $self, @args ) = @_;
    my $id_job      = shift @args;
    my $job_dir     = $self->jobId_to_jobPath($id_job);
    my $is_finished_file = File::Spec->catfile( $job_dir, 'finished' );
    return (-e $is_finished_file);
}


=method get_candidates_dir


=cut

sub get_candidates_dir {
	my ( $self, @args ) = @_;
	my $id_job         = shift @args;
	my $results_dir    = $self->jobId_to_jobPath($id_job);
	my $candidates_dir = File::Spec->catdir( $results_dir, 'candidates' );
	return $candidates_dir;
}

=method is_valid_jobID

Test whether a jobID is valid - ie if there are results for it.

=cut

sub is_valid_jobID {
	my ( $self, @args ) = @_;
	my $id_job    = shift @args;
	my $full_path = $self->jobId_to_jobPath($id_job);
	return ( -e $full_path );
}

=method get_structure_for_jobID

Get the results structure of a given job identifier

Usage:
my %results = miRkwood::Results->get_structure_for_jobID($jobId);

=cut

sub get_structure_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
	my $job_dir = $self->jobId_to_jobPath($jobId);
	miRkwood->CONFIG_FILE(
		miRkwood::Paths->get_job_config_path($job_dir) );
	my $candidates_dir = $self->get_candidates_dir($jobId);
	return $self->deserialize_results($candidates_dir);
}

=method get_basic_structure_for_jobID

=cut

sub get_basic_structure_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
    my $type    = shift @args;  # $type should be 'New' or 'Known'
	my $job_dir = $self->jobId_to_jobPath($jobId);
    my $yml_file = '';
    if ( $type eq 'Known' ){
        $yml_file = 'basic_known_candidates.yml';
    }
    else{
        $yml_file = 'basic_candidates.yml';
    }
	miRkwood->CONFIG_FILE(
		miRkwood::Paths->get_job_config_path($job_dir) );
	my $candidates_file = File::Spec->catfile( $job_dir, $yml_file);
	return miRkwood::get_yaml_file( $candidates_file );
}

sub get_basic_pseudoXML_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
    my $type    = shift @args;  # $type should be 'New' or 'Known'

	my $results = $self->get_basic_structure_for_jobID($jobId, $type);

	my $output = '';

    $output .= "<results id='all'>\n";
    my @candidates = sort {
        ( $a->{'name'} cmp  $b->{'name'} )
          || (
            $a->{'start_position'} <=> $b->{'start_position'} )
    } @{$results};

    foreach my $candidate (@candidates) {
        $output .= $self->convert_basic_to_pseudoXML($candidate, $type) . "\n";
    }
    $output .= "</results>\n";

    $output .= "<results id='all2'>\n";
    @candidates = sort {
        ( $b->{'quality'} cmp $a->{'quality'} )
          || (
            $a->{'start_position'} <=> $b->{'start_position'} )
    } @candidates;
    foreach my $candidate (@candidates) {
        $output .= $self->convert_basic_to_pseudoXML($candidate, $type) . "\n";
    }
    $output .= '</results>';

    return $output;
}

sub convert_basic_to_pseudoXML {
	my ( $self, @args ) = @_;
	my $candidate = shift @args;
	my %candidate = %{$candidate};
    my $type      = shift @args;  # $type should be 'New' or 'Known'

	my @headers;
    my @fields_to_truncate = qw{mfe mfei amfe};
    my $result = '<Sequence';
	my @optional_fields = miRkwood::Candidate->get_optional_candidate_fields();

    if ( $type eq 'Known' ){
        push @headers, 'precursor_name';
        @optional_fields = miRkwood::Utils::delete_element_in_array( 'alignment', @optional_fields );
    }

    push @headers, ( 'name', 'position', 'length', 'strand', 'quality', @optional_fields, 'image', 'identifier' ); 
    for my $header (@headers) {
        my $contents = $candidate->{$header};
        if (grep { $header eq $_ } @fields_to_truncate){
            $contents = miRkwood::Utils::restrict_num_decimal_digits($contents, 3);
        }
        if ( !defined $contents ) {
            $contents = q{};
        }
        if ( $header eq 'shuffles' && $contents == 1){
            $contents = q{};
        }
        $result .= " $header='$contents'";
    }
    $result .= '></Sequence>';
}

=method has_candidates

Parse and serialize the results structure of $job_dir

Usage:
miRkwood::Results->has_candidates( \%myResults );

=cut

sub has_candidates {
	my ( $self, @args ) = @_;
	my $results = shift @args;
	my %results = %{$results};
	return ( keys %results > 0 );
}

=method deserialize_results

Retrieve the results in the given directory

Usage:
my %results = miRkwood::Results->deserialize_results( $candidates_dir );

=cut

sub deserialize_results {
	my ( $self, @args ) = @_;
	my $candidates_dir = shift @args;
	my %myResults      = ();
	opendir DIR, $candidates_dir;    #ouverture répertoire job
	my @files;
	@files = readdir DIR;
	closedir DIR;
	foreach my $file (@files)        # parcours du contenu
	{
		my $full_file = File::Spec->catfile( $candidates_dir, $file );
		if (   $file ne '.'
			&& $file ne '..' )
		{
			my $candidate;
			if (
				eval {
					$candidate =
					  miRkwood::Candidate->new_from_serialized($full_file);
				}
			  )
			{
				my $identifier = $candidate->get_identifier();
				$myResults{$identifier} = $candidate;
			}
		}
	}
	return %myResults;
}

=method number_of_results

return total number of candidates 

=cut

sub number_of_results {
	my ( $self, @args ) = @_;
	my $results = shift @args;
	my %results = %{$results};
	my $size    = scalar keys %results;
	return $size;
}

sub number_of_results_bis {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
    my $type    = shift @args;  # $type should be 'New' or 'Known'
	my $results = $self->get_basic_structure_for_jobID($jobId, $type);
	my $size    = scalar @{$results};
	return $size;
}

1;
