package PipelineMiRNA::Paths;

use strict;
use warnings;
use Cwd;
use File::Spec;
use File::Basename;
use YAML;

=method get_root_dir

Get the YAML configuration file contents.

=cut

sub get_config {
    my ($self, @args) = @_;
    my $config_file = File::Spec->catfile(File::Basename::dirname(__FILE__), 'pipeline.cfg');
    return %{YAML::LoadFile($config_file )};
}

=method get_root_dir

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

=method get_server_path

=cut

sub get_server_path {
    my ($self, @args) = @_;
    return File::Spec->catdir( $self->get_server_root_dir(), @args );
}

=method get_absolute_path

=cut

sub get_absolute_path {
    my ($self, @args) = @_;
    return File::Spec->catdir( $self->get_root_dir(), @args );
}

1;
