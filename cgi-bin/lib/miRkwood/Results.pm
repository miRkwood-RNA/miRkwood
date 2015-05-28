package miRkwood::Results;

# ABSTRACT: Code directly tied to the results data structure

use strict;
use warnings;

use feature 'switch';
use Time::gmtime;
use File::Spec;
use YAML::XS;

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
my %results = miRkwood::Results->get_structure_for_jobID($jobId, $type);
$type should be 'Known', 'New' (for BEDPipeline) or '' (for other pipelines)

=cut

sub get_structure_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
	my $mirna_type = shift @args;   # should be 'Known', 'New' (for BEDPipeline) or '' (for other pipelines)
	my $job_dir = $self->jobId_to_jobPath($jobId);
	miRkwood->CONFIG_FILE(
		miRkwood::Paths->get_job_config_path($job_dir) );
	my $candidates_dir = '';
    if ( $mirna_type eq 'Known' ){
        $candidates_dir = miRkwood::Paths::get_known_candidates_dir_from_job_dir($job_dir);
    }
    elsif ($mirna_type eq 'New' ){
        $candidates_dir = miRkwood::Paths::get_new_candidates_dir_from_job_dir($job_dir);
    }
    else{
        $candidates_dir = miRkwood::Paths::get_dir_candidates_path_from_job_dir($job_dir);
    }
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
        @optional_fields = miRkwood::Utils::delete_element_in_array( 'alignment', \@optional_fields );
        @optional_fields = miRkwood::Utils::delete_element_in_array( 'mfe', \@optional_fields );
        @optional_fields = miRkwood::Utils::delete_element_in_array( 'mfei', \@optional_fields );
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
        if ($header eq 'reads'){
            $contents = 0;
            foreach my $key (keys( %{$candidate->{'reads'}} )){
                $contents += $candidate->{'reads'}{$key};
            }
        }
        $result .= " $header='$contents'";
    }
    $result .= '></Sequence>';

    return $result;
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

sub count_reads_in_basic_yaml_file {
    my ( $self, @args ) = @_;
    my $yaml = shift @args;
    my $nb_reads = 0;
    my $nb_reads_unq = 0;

    my @attributes = YAML::XS::LoadFile($yaml);
    my $reads;

    for (my $i = 0; $i < scalar(@{$attributes[0]}); $i++){
        foreach my $positions ( keys%{$attributes[0][$i]{'reads'}} ){
            $reads->{ $positions } = $attributes[0][$i]{'reads'}{$positions};
        }
    }
    foreach my $positions (keys%{$reads}){
        $nb_reads += $reads->{$positions};
    }
    $nb_reads_unq = scalar (keys%{$reads});

    return ($nb_reads, $nb_reads_unq);
}

sub make_reads_barchart {
    my ( $self, $total_witdh, $percentage_CDS_reads, $percentage_other_reads, $percentage_multi_reads, $percentage_known_miRNAs_reads, $percentage_new_miRNAs_reads ) = @_;

    my $width_CDS_reads          = int($percentage_CDS_reads * $total_witdh / 100 + 0.5 );
    my $width_other_reads        = int($percentage_other_reads * $total_witdh / 100 + 0.5 );
    my $width_multi_reads        = int($percentage_multi_reads * $total_witdh / 100 + 0.5 );
    my $width_known_miRNAs_reads = int($percentage_known_miRNAs_reads * $total_witdh / 100 + 0.5 );
    my $width_new_miRNAs_reads   = int($percentage_new_miRNAs_reads * $total_witdh / 100 + 0.5 );
    my $width_orphans_reads      = 100 - $percentage_CDS_reads - $percentage_other_reads - $percentage_multi_reads - $percentage_known_miRNAs_reads - $percentage_new_miRNAs_reads;

    my $barchart = <<"END_TXT";
<div style='width:${total_witdh}px'>
    <table id="barchart_table">
        <tr>
            <td id="CDS" style="width:${width_CDS_reads}px"></td>
            <td id="other" style="width:${width_other_reads}px"></td>
            <td id="multimapped" style="width:${width_multi_reads}px;"></td>
            <td id="known_miRNAs" style="width:${width_known_miRNAs_reads}px;"></td>
            <td id="new_miRNAs" style="width:${width_new_miRNAs_reads}px;"></td>
            <td id="orphans" style="width:${width_orphans_reads}px;"></td>
        </tr>
    </table>

    <table id="barchart_legend">
        <tr>
            <td><span id="CDS">&nbsp;&nbsp;&nbsp;</span> CDS</td>
            <td><span id="other">&nbsp;&nbsp;&nbsp;</span> tRNA/rRNA/snoRNA</td>
            <td><span id="multimapped">&nbsp;&nbsp;&nbsp;</span> multiply mapped reads</td>
        </tr><tr>
            <td><span id="known_miRNAs">&nbsp;&nbsp;&nbsp;</span> knowns miRNAs</td>
            <td><span id="new_miRNAs">&nbsp;&nbsp;&nbsp;</span> new miRNAs</td>
            <td><span id="orphans">&nbsp;&nbsp;&nbsp;</span> orphans reads</td>
        </tr>  
    </table>
</div>
END_TXT

    return $barchart;

}

1;
