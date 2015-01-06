package miRkwood::BamPipelineSebastien;

# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use parent 'miRkwood::PipelineSebastien';

use Log::Message::Simple qw[msg error debug];

use miRkwood::ClustersSebastien;

=method new

Constructor

=cut

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
# SEB BEGIN
sub init_sequences {
    my ($self, @args) = @_;
    debug( "Extracting sequences from genome using BAM clusters", miRkwood->DEBUG() );
    my $clustering = miRkwood::ClustersSebastien->new($self->{'genome_file'});
    my $reads = $clustering->get_read_distribution_from_bed($self->{'bed_file'});
    my $sequences = $clustering->get_windows($reads, 2);
    $self->{'sequences'} = $sequences;
    $self->{'clustering'} = $clustering;
    return;
}
# SEB END
1;
