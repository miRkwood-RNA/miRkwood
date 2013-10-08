package PipelineMiRNA::Paths;

# ABSTRACT: Managing paths and path construction

use strict;
use warnings;
use Cwd;
use File::Spec;
use File::Basename;
use YAML;

=method get_config

Get the YAML configuration file contents.

=cut

sub get_config {
    my ($self, @args) = @_;
    my $config_file = File::Spec->catfile(File::Basename::dirname(__FILE__), 'pipeline.cfg');
    return %{YAML::LoadFile($config_file )};
}

=method get_server_root_dir

Return the project root directory, as considered by the server

=cut

sub get_server_root_dir {
    my ($self, @args) = @_;
    my %config = $self->get_config();
    return $config{'server_root'};
}

=method get_root_dir

Return the project root directory

=cut

sub get_root_dir {
    my ($self, @args) = @_;
    my %config = $self->get_config();
    return $config{'absolute_root'};
}

=method get_results_dir_name

Return the project results (relatively to the root)

=cut

sub get_results_dir_name {
    my ($self, @args) = @_;
    my %config = $self->get_config();
    return $config{'results'};
}

=method get_server_path

Append the given path to the server root

=cut

sub get_server_path {
    my ($self, @args) = @_;
    return File::Spec->catdir( $self->get_server_root_dir(), @args );
}

=method get_absolute_path

Append the given path to the absolute root

=cut

sub get_absolute_path {
    my ($self, @args) = @_;
    return File::Spec->catdir( $self->get_root_dir(), @args );
}

=method get_candidate_paths

Return both the server and absolute paths for a given candidate

=cut

sub get_candidate_paths {
    my ($self, @args) = @_;
    my ($job,  $dir, $subDir) = @args;
    my $candidate_dir = File::Spec->catdir($job,  $dir, $subDir);
    my $full_candidate_dir = $self->get_absolute_path($candidate_dir);
    return ($candidate_dir, $full_candidate_dir);
}

1;
