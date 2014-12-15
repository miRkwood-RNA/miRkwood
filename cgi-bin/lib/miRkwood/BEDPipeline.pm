package miRkwood::BEDPipeline;

# ABSTRACT: BED Pipeline object

use strict;
use warnings;

use parent 'miRkwood::Pipeline';

use Log::Message::Simple qw[msg error debug];

sub new {
    my ( $class, @args ) = @_;
    my ($job_dir, $bed_file, $genome_file) = @args;
    my $self = {
        job_dir => $job_dir,
        bed_file => $bed_file,
        genome_file => $genome_file,
        sequences => undef
    };
    bless $self, $class;
    return $self;
}

=method init_sequences

=cut

sub init_sequences {
    my ($self, @args) = @_;
    debug( "~~~ Waiting for the new Clusters.pm ~~~", miRkwood->DEBUG() );
    #~ debug( "Sequences have been extracted from BED file.", miRkwood->DEBUG() );
    #~ my $clustering = miRkwood::Clusters->new($self->{'bam_file'}, $self->{'genome_file'});
    #~ my @sequences_array = $clustering->get_clustered_sequences_from_bam();
    #~ $self->{'sequences'} = \@sequences_array;
    return;
}

1;
