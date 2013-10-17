package PipelineMiRNA;

use strict;
use warnings;

use Config::Simple;

# ABSTRACT: A Pipeline for microRNAs

our $LOG_FH;

=method LOGFH

Set and get the logging architecture

=cut

sub LOGFH{
    my ($self, @args) = @_;
    if (@args == 1){
        my $log_file = shift @args;
        open( my $LOG, '>>', $log_file ) || die "Error when opening log file $log_file: $!";
        $LOG_FH = $LOG;
    }
    return $LOG_FH;
}


our $CONFIG_FILE;
our $CONFIG;

=method CONFIG_FILE

=cut

sub CONFIG_FILE{
    my ($self, @args) = @_;
    if (@args == 1){
        $CONFIG_FILE = shift @args;
        $CONFIG = new Config::Simple(syntax=>'ini');
        $CONFIG->read($CONFIG_FILE);
    }
    return $CONFIG_FILE;
}

=method CONFIG

=cut

sub CONFIG{
    my ($self, @args) = @_;
    if (@args == 1){
        $CONFIG = shift @args;
        $CONFIG->write($CONFIG_FILE);
    }
    return $CONFIG;
}

1;
