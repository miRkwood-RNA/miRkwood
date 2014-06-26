package PipelineMiRNA::Paths;

# ABSTRACT: Managing paths and path construction

use strict;
use warnings;
use Cwd;
use File::Spec;
use File::Basename;

use PipelineMiRNA;

=method get_config

Get the configuration file contents.

=cut

sub get_config {
    my ($self, @args) = @_;
    return PipelineMiRNA->PIPELINE_CONFIG();
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
    return File::Spec->catdir(PipelineMiRNA->MIRKWOOD_PATH(), 'data');
}

=method get_local_programs_path

Return the project local programs directory

=cut

sub get_local_programs_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(PipelineMiRNA->MIRKWOOD_PATH(), 'programs');
}

=method get_scripts_path

Return the project static directory

=cut

sub get_scripts_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(PipelineMiRNA->MIRKWOOD_PATH(), 'scripts');
}

=method get_lib_path

Return the project static directory

=cut

sub get_lib_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(PipelineMiRNA->MIRKWOOD_PATH(), 'lib');
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

1;
