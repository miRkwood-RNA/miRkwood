package miRkwood::FastaPipeline;

# ABSTRACT: FASTA Pipeline object

use strict;
use warnings;

use parent 'miRkwood::Pipeline';

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($job_dir) = @args;
    my $self = {
        job_dir => $job_dir,
        sequences => undef
    };
    bless $self, $class;
    return $self;
}

1;