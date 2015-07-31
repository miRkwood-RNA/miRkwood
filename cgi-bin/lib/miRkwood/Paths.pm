package miRkwood::Paths;

# ABSTRACT: Managing paths and path construction

use strict;
use warnings;
use Cwd;
use File::Spec;
use File::Basename;

use miRkwood;
use miRkwood::Results;

=method get_config

Get the configuration file contents.

=cut

sub get_config {
    my ($self, @args) = @_;
    return miRkwood->PIPELINE_CONFIG();
}

=method get_job_config_path

Given a job directory, return the path to the job configuration file

=cut

sub get_job_config_path {
    my ($self, @args) = @_;
    my ($dir_job) = shift @args;
    my $job_config_path = File::Spec->catfile( $dir_job, 'run_options.cfg' );
    return $job_config_path;
}

=method get_data_path

Return the project data directory

=cut

sub get_data_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(miRkwood->MIRKWOOD_PATH(), 'data');
}

=method get_local_programs_path

Return the project local programs directory

=cut

sub get_local_programs_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(miRkwood->MIRKWOOD_PATH(), 'programs');
}

=method get_scripts_path

Return the project static directory

=cut

sub get_scripts_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(miRkwood->MIRKWOOD_PATH(), 'scripts');
}

=method get_lib_path

Return the project static directory

=cut

sub get_lib_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(miRkwood->MIRKWOOD_PATH(), 'lib');
}

=method get_results_filesystem_path

Return the path to the results, as seen from the get_results_filesystem_path

=cut

sub get_results_filesystem_path {
    my ($self, @args) = @_;
    my %config = $self->get_config();
    return $config{'filesystem_results'};
}

=method get_workspace_path

Return the path to the job workspace

=cut

sub get_workspace_path {
    my ($self, @args) = @_;
    my ($job_dir) = @args;
    return File::Spec->catdir($job_dir, 'workspace');
}

=method get_candidate_paths

Return the path for a given candidate

=cut

sub get_candidate_paths {
    my ($self, @args) = @_;
    my ($job_dir,  $dir, $subDir) = @args;
    my $workspace = $self->get_workspace_path($job_dir);
    my $candidate_dir = File::Spec->catdir($workspace, $dir, $subDir);
    return $candidate_dir;
}

=method get_dir_candidates_path

Return the path to the candidates directory
Parameter : job id

=cut

sub get_dir_candidates_path {
    my (@args) = @_;
    my $job_id = shift @args;
    my $job_dir = miRkwood::Results->jobId_to_jobPath( $job_id );
    return File::Spec->catdir($job_dir, 'candidates');
}

=method get_dir_candidates_path

Return the path to the candidates directory
Parameter : job directory

=cut
sub get_dir_candidates_path_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir($job_dir, 'candidates');
}

=method get_new_candidates_dir_from_job_dir

Return the path to the new candidates directory
Parameter : job id

=cut
sub get_new_candidates_dir {
    my (@args) = @_;
    my $job_id = shift @args;
    return File::Spec->catdir( get_dir_candidates_path($job_id), 'new');
}

=method get_new_candidates_dir_from_job_dir

Return the path to the new candidates directory
Parameter : job directory

=cut
sub get_new_candidates_dir_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir( get_dir_candidates_path_from_job_dir($job_dir), 'new');
}

=method get_known_candidates_dir

Return the path to the new candidates directory
Parameter : job id

=cut
sub get_known_candidates_dir {
    my (@args) = @_;
    my $job_id = shift @args;
    return File::Spec->catdir( get_dir_candidates_path($job_id), 'known');
}

=method get_known_candidates_dir_from_job_dir

Return the path to the new candidates directory
Parameter : job directory

=cut
sub get_known_candidates_dir_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir( get_dir_candidates_path_from_job_dir($job_dir), 'known');
}

=method get_dir_reads_path

Return the path to the reads directory
Parameter : job id

=cut

sub get_dir_reads_path {
    my (@args) = @_;
    my $job_id = shift @args;
    my $job_dir = miRkwood::Results->jobId_to_jobPath( $job_id );
    return File::Spec->catdir($job_dir, 'reads');
}

=method get_dir_reads_path_from_job_dir

Return the path to the reads directory
Parameter : job directory

=cut

sub get_dir_reads_path_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir($job_dir, 'reads');
}

=method get_new_reads_dir

Return the path to the new reads directory
Parameter : job id

=cut
sub get_new_reads_dir {
    my (@args) = @_;
    my $job_id = shift @args;
    return File::Spec->catdir( get_dir_reads_path($job_id), 'new');
}

=method get_new_reads_dir_from_job_dir

Return the path to the new reads directory
Parameter : job directory

=cut
sub get_new_reads_dir_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir( get_dir_reads_path_from_job_dir($job_dir), 'new');
}

=method get_known_reads_dir

Return the path to the new reads directory
Parameter : job id

=cut
sub get_known_reads_dir {
    my (@args) = @_;
    my $job_id = shift @args;
    return File::Spec->catdir( get_dir_reads_path($job_id), 'known');
}

=method get_known_reads_dir_from_job_dir

Return the path to the new reads directory
Parameter : job id

=cut
sub get_known_reads_dir_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir( get_dir_reads_path_from_job_dir($job_dir), 'known');
}

=method get_dir_reads_path_from_job_dir_and_mirna_type

Return the path to the 'new' or 'known' reads directory
Parameters : - job dir
             - mirna type ('New' or 'Known')

=cut
sub get_dir_reads_path_from_job_dir_and_mirna_type{
    my (@args) = @_;
    my $job_dir = shift @args;
    my $mirna_type = shift @args;
    if ( $mirna_type eq 'Known' ){
        return get_known_reads_dir_from_job_dir( $job_dir );
    }
    else{
        return get_new_reads_dir_from_job_dir( $job_dir );
    }
}

=method get_bed_file

Return the path to the BED file corresponding to given type
'', '_filtered', '_CDS, '_miRNAs', '_otherRNA', '_multimapped'
for the given job ID.
Parameter : job id and bed type

=cut
sub get_bed_file {
    my (@args) = @_;
    my $jobID    = shift @args;
    my $bed_type = shift @args; # should be '', '_filtered', '_CDS, '_miRNAs', '_otherRNA', '_multimapped'
    my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($jobID);
    my $run_options_file = miRkwood::Paths->get_job_config_path($absolute_job_dir);
    miRkwood->CONFIG_FILE($run_options_file);
    my $cfg      = miRkwood->CONFIG();
    my $bed_name = $cfg->param('job.bed');

    return File::Spec->catdir( $absolute_job_dir, $bed_name. $bed_type . '.bed');
}

=method get_workspace_chromosome_dir

  Return the path to a given chromosome directory
  in the workspace directory
  Parameters : - workspace directory
               - chromosome id

=cut
sub get_workspace_chromosome_dir {
    my (@args) = @_;
    my $workspace_dir = shift @args;
    my $chromosome = shift @args;
    return File::Spec->catdir( $workspace_dir, $chromosome);
}

=method get_workspace_candidate_dir

  Return the path to a Candidate directory
  in the workspace directory
  Parameters : - workspace directory
               - chromosome id
               - cluster positions
               - strand  

=cut
sub get_workspace_candidate_dir {
    my (@args) = @_;
    my $workspace_dir = shift @args;
    my $chromosome = shift @args;
    my $cluster = shift @args;
    my $strand = shift @args;

    my $chromosome_dir = get_workspace_chromosome_dir( $workspace_dir, $chromosome );
    return File::Spec->catdir( $chromosome_dir, $cluster . $strand );
}

=method create_folder

  Method to create a folder if it not already exists
  Parameter : the folder name

=cut
sub create_folder {
    my (@args) = @_;
    my $folder = shift @args;
    if ( !-e $folder ) {
        mkdir $folder or die "ERROR when creating $folder : $!";
    }
    return $folder;
}

1;
