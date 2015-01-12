package miRkwood::BamPipelineSebastien;

# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use parent 'miRkwood::Pipeline';

use Log::Message::Simple qw[msg error debug];

use miRkwood::ClustersSebastien;

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($job_dir, $bed_file, $genome_file) = @args;
    my %genome_db = miRkwood::Utils::multifasta_to_hash( $genome_file );
    my $self = {
        job_dir => $job_dir,
        bed_file => $bed_file,
        genome_file => $genome_file,
        genome_db => \%genome_db,
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
    my ($reads, $parsed_bed) = $clustering->get_read_distribution_from_bed($self->{'bed_file'});
    my $sequences = $clustering->get_windows($reads, 2);
    $self->{'sequences'} = $sequences;
    $self->{'parsed_reads'} = $parsed_bed;
    $self->{'clustering'} = $clustering;
    return;
}
# SEB END
1;
