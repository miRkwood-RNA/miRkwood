package miRkwood::BamPipelineSebastien;

# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use parent 'miRkwood::Pipeline';

use Log::Message::Simple qw[msg error debug];

use miRkwood::ClusterBuilder;

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($job_dir, $bed_file, $genome_file) = @args;
    my %genome_db = miRkwood::Utils::multifasta_to_hash( $genome_file );
    my $self = {
        job_dir     => $job_dir,
        initial_bed => $bed_file,       # this is for the non filtered BED provided by the user
        bed_file    => '',              # this is for the final BED after all filtrations have been done
        mirna_bed   => '',              # this is for known miRNAs (ie miRNAs from miRBase). created during the filtration step. Maybe change that name ?
        genome_file => $genome_file,
        genome_db   => \%genome_db,     # hash containing the genome, to avoid reading the file each time we need it
        sequences   => undef
    };
    bless $self, $class;
    return $self;
}

=method init_sequences

=cut
# SEB BEGIN
sub init_sequences {
    my ($self, @args) = @_;
    debug( "Extracting sequences from genome using BED clusters", miRkwood->DEBUG() );
    my $clustering = miRkwood::ClusterBuilder->new($self->{'genome_db'}, $self->{'bed_file'});
    #~ my ($reads, $parsed_bed) = $clustering->get_read_distribution_from_bed($self->{'bed_file'});
    #~ my $sequences = $clustering->get_windows($reads, 2);
    #~ $self->{'sequences'} = $sequences;
    $self->{'sequences'} = $clustering->build_loci();
    $self->{'parsed_reads'} = $clustering->get_parsed_bed();
    return;
}
# SEB END
1;
