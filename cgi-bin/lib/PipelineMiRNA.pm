package PipelineMiRNA;

use strict;
use warnings;

use Config::Simple;
use File::Basename;
use CGI::Carp qw(carpout);
use YAML::XS;

# ABSTRACT: A Pipeline for microRNAs

our $LOG_FH;

=method LOGFH

Set and get the logging architecture

=cut

sub LOGFH {
    my ( $self, @args ) = @_;
    if ( @args == 1 ) {
        my $log_file = shift @args;
        open( my $LOG, '>>', $log_file )
          || die "Error when opening log file $log_file: $!";
        $LOG_FH = $LOG;
        carpout($LOG);
    }
    return $LOG_FH;
}

our $DEBUG;

=method DEBUG

Set and get the debug verbosity

=cut

sub DEBUG {
    my ( $self, @args ) = @_;
    if ( @args == 1 ) {
        $DEBUG = shift @args;
    }
    return $DEBUG;
}

our $CONFIG_FILE;
our $CONFIG;

=method CONFIG_FILE

=cut

sub CONFIG_FILE {
    my ( $self, @args ) = @_;
    if ( @args == 1 ) {
        $CONFIG_FILE = shift @args;
        $CONFIG = Config::Simple->new( syntax => 'ini' );
        $CONFIG->read($CONFIG_FILE);
    }
    return $CONFIG_FILE;
}

=method CONFIG

=cut

sub CONFIG {
    my ( $self, @args ) = @_;
    if ( @args == 1 ) {
        $CONFIG = shift @args;
        $CONFIG->write($CONFIG_FILE);
    }
    return $CONFIG;
}

our $PIPELINE_CONFIG_FILE;
our %PIPELINE_CONFIG;

=method PIPELINE_CONFIG_FILE

=cut

sub PIPELINE_CONFIG_FILE {
    my (@args) = @_;
    if ( @args == 1 ) {
        $PIPELINE_CONFIG_FILE = shift @args;
        my $yaml = get_yaml_file($PIPELINE_CONFIG_FILE);
        %PIPELINE_CONFIG = %{$yaml};
    }
    return $PIPELINE_CONFIG_FILE;
}

=method PIPELINE_CONFIG

=cut

sub PIPELINE_CONFIG {
    my ( $self, @args ) = @_;
    if ( @args == 1 ) {

        # We do not allow to overwrite the config yet
    }
    return %PIPELINE_CONFIG;
}

my $default_pipeline_cfg_file =
  File::Spec->catfile( File::Basename::dirname(__FILE__),
    'PipelineMiRNA', 'pipeline.cfg' );
PIPELINE_CONFIG_FILE($default_pipeline_cfg_file);

our $PROGRAMS_CONFIG_FILE;
our %PROGRAMS_CONFIG;

=method PIPELINE_CONFIG_FILE

=cut

sub PROGRAMS_CONFIG_FILE {
    my (@args) = @_;
    if ( @args == 1 ) {
        $PROGRAMS_CONFIG_FILE = shift @args;
        my $yaml = get_yaml_file($PROGRAMS_CONFIG_FILE);
        %PROGRAMS_CONFIG = %{$yaml};
    }
    return $PROGRAMS_CONFIG_FILE;
}

=method PROGRAMS_CONFIG

=cut

sub PROGRAMS_CONFIG {
    my ( $self, @args ) = @_;
    return %PROGRAMS_CONFIG;
}

my $default_programs_cfg_file =
  File::Spec->catfile( File::Basename::dirname(__FILE__),
    'PipelineMiRNA', 'programs.cfg' );
PROGRAMS_CONFIG_FILE($default_programs_cfg_file);

=method get_yaml_file

Check whether a given file is suitable for YAML deserialization,
and returns its content if so.

=cut

sub get_yaml_file {
    my @args      = @_;
    my $yaml_file = shift @args;
    ( -e $yaml_file )  or die("Error, YAML file $yaml_file does not exist");
    ( !-z $yaml_file ) or die("Error, YAML file $yaml_file is empty");
    ( -r $yaml_file )  or die("Error, YAML file $yaml_file is not readable");
    my $yaml;
    if ( !eval { $yaml = YAML::XS::LoadFile($yaml_file); } ) {
        die("Error, YAML file $yaml_file is malformed");
    }
    else {
        return $yaml;
    }
}

1;
