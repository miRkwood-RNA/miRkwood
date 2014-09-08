package miRkwood::SequenceJob;

# ABSTRACT: Job processing a given sequence

use strict;
use warnings;


=method new

Constructor

my $sequence_job = miRkwood::SequenceSubJob->new($sequence_dir, $name, $sequence, '+');

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($directory, $name, $sequence, $strand) = @args;
    my $self = {
    };
    bless $self, $class;
    return $self;
}

1;