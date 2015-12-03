package miRkwood::WebPaths;

# ABSTRACT: Managing Web paths and path construction

use strict;
use warnings;
use File::Spec;

use miRkwood;
use miRkwood::Paths;

=method get_web_config

Get the configuration file contents.

=cut

sub get_web_config {
    my ($self, @args) = @_;
    return miRkwood->WEB_CONFIG();
}

=method get_static_path

Return the project static directory

=cut

sub get_static_path {
    my ($self, @args) = @_;
    my %config = $self->get_web_config();
    return $config{'static'};
}

=method get_web_root

Return the web root
(ie in http://myserver.org/web_root/<pipeline.pl>)

=cut

sub get_web_scripts {
    my ($self, @args) = @_;
    my %config = $self->get_web_config();
    return $config{'web_scripts'};
}

=method get_css_path

Return the project CSS directory

=cut

sub get_css_path {
    my ($self, @args) = @_;
    my %config = $self->get_web_config();
    return File::Spec->catdir($config{'html_root'}, 'style');
}

=method get_html_path

Return the html directory

=cut

sub get_html_path {
    my ($self, @args) = @_;
    my %config = $self->get_web_config();
    return $config{'html_root'};
}

=method get_server_css_path

Return the Server CSS directory

=cut

sub get_server_css_path {
    my ($self, @args) = @_;
    my %config = $self->get_web_config();
    return $config{'server_css'};
}

=method get_js_path

Return the project JavaScript directory

=cut

sub get_js_path {
    my ($self, @args) = @_;
    my %config = $self->get_web_config();
    return File::Spec->catdir($config{'html_root'}, 'js');
}

=method get_server_scripts_path

Return the Server scripts directory

=cut
sub get_server_scripts_path {
    my ($self, @args) = @_;
    my %config = $self->get_web_config();
    return $config{'server_scripts'};    
}

=method get_server_libs_path

Return the Server libs directory

=cut
sub get_server_libs_path {
    my ($self, @args) = @_;
    my %config = $self->get_web_config();
    return $config{'libs'};   
}

=method get_results_web_path

Return the path to the results, as seen from the web

=cut

sub get_results_web_path {
    my ($self, @args) = @_;
    my %config = $self->get_web_config();
    return $config{'web_results'};
}

=method filesystem_to_web_path

Convert a filesystem path to a web path

=cut

sub filesystem_to_relative_path {
    my ($self, @args) = @_;
    my $path = shift @args;
    my $filesystem_path = miRkwood::Paths->get_results_filesystem_path();
    my $web_path        = $self->get_results_web_path();
    $path =~ s/$filesystem_path/$web_path/g;
    return $path;
}

1;
