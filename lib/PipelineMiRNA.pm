package PipelineMiRNA;

use strict;
use warnings;

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

1;
