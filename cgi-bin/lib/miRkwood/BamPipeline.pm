package miRkwood::BamPipeline;

# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use parent 'miRkwood::Pipeline';

use Log::Message::Simple qw[msg error debug];

use miRkwood::Clusters;

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($job_dir, $bam_file, $genome_file) = @args;
    my %genome_db = miRkwood::Utils::multifasta_to_hash( $genome_file );
    my $self = {
        job_dir => $job_dir,
        bam_file => $bam_file,
        genome_file => $genome_file,
        genome_db   => \%genome_db,     # hash containing the genome, to avoid reading the file each time we need it
        sequences => undef
    };
    bless $self, $class;
    return $self;
}

=method init_sequences

=cut

sub init_sequences {
    my ($self, @args) = @_;
    debug( "Extracting sequences from genome using BAM clusters", miRkwood->DEBUG() );
    my $clustering = miRkwood::Clusters->new($self->{'bam_file'}, $self->{'genome_file'});
    my @sequences_array = $clustering->get_clustered_sequences_from_bam();
    $self->{'sequences'} = \@sequences_array;
    return;
}

sub serialize_candidates {
    my ($self, @args) = @_;
    my $candidates = shift @args;
    my @candidates_array = @{$candidates};

    foreach my $candidate (@candidates_array ) {

        $candidate->turn_relative_positions_into_absolute_positions();
        $candidate = $candidate->get_reads_from_bam_file($self->{'bam_file'});
        miRkwood::CandidateHandler::print_reads_clouds( $candidate, $self->get_new_reads_dir() );
        miRkwood::CandidateHandler->serialize_candidate_information( $self->get_candidates_dir(), $candidate );

        push $self->{'basic_candidates'}, $candidate->get_basic_informations();
    }
    return;
}

1;
