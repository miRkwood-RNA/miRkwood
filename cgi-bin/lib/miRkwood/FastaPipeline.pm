package miRkwood::FastaPipeline;

# ABSTRACT: FASTA Pipeline object

use strict;
use warnings;

use parent 'miRkwood::Pipeline';

use Log::Message::Simple qw[msg error debug];

use miRkwood::Maskers;
use miRkwood::Paths;
use miRkwood::Utils;

### Data ##
my $dirData = miRkwood::Paths->get_data_path();

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


=method init_sequences

Get the sequences to process from the job directory
This includes parsing, and masking if option selected.

 Usage : my @fasta_array = get_sequences($job_dir);
 Return: -

=cut

sub init_sequences {
    my ($self, @args) = @_;
    debug( 'Getting sequences', miRkwood->DEBUG() );
    my %filter = $self->get_masking_information();
    my $sequence_uploaded = $self->get_uploaded_sequences_file();

    open my $ENTREE_FH, '<', $sequence_uploaded
      or die "Error when opening sequences -$sequence_uploaded-: $!";
    debug( "Calling parse_multi_fasta() on $sequence_uploaded", miRkwood->DEBUG() );
    my @sequences_array = miRkwood::Utils::parse_multi_fasta($ENTREE_FH);
    close $ENTREE_FH;
    my @sequences;
    if (%filter){
        debug( 'Masking input sequences', miRkwood->DEBUG() );
        @sequences = miRkwood::Maskers::mask_sequences(\%filter, @sequences_array);
    } else {
        @sequences = @sequences_array;
    }
    $self->{'sequences'} = \@sequences;
    return;
}

=method get_masking_information

Get the masking information based on the job configuration

 Usage : my %filter = get_masking_information($job_dir);
 Input : The job directory
 Return: A hash (name => [ positions ])

=cut

sub get_masking_information {
    my ($self, @args) = @_;
    my %filter;
    my $cfg = miRkwood->CONFIG();
    if ( $cfg->param('options.filter') || $cfg->param('options.mask-trna') || $cfg->param('options.mask-rrna') ){
        my $masking_folder = miRkwood::Paths::create_folder( File::Spec->catdir($self->get_job_dir(), 'masks') );
        my $sequences = $self->get_uploaded_sequences_file();
    
        if ($cfg->param('options.filter')) {
            debug( 'Get masking information for coding regions', miRkwood->DEBUG() );
            my $plant = $cfg->param('job.plant');
            my $blast_database = File::Spec->catfile( $dirData, "$plant.fas" );
            my %blast_mask = miRkwood::Maskers::get_coding_region_masking_information( $sequences, $masking_folder, $blast_database );
            %filter = miRkwood::Utils::merge_hashes_of_arrays(\%filter, \%blast_mask);
        }
        if ($cfg->param('options.mask-trna')) {
            debug( 'Get masking information for tRNAs', miRkwood->DEBUG() );
            my %trna_mask = miRkwood::Maskers::get_trna_masking_information( $sequences, $masking_folder );
            %filter = miRkwood::Utils::merge_hashes_of_arrays(\%filter, \%trna_mask);
        }
        if ($cfg->param('options.mask-rrna')) {
            debug( 'Get masking information for ribosomal RNAs', miRkwood->DEBUG() );
            my %rrna_mask = miRkwood::Maskers::get_rnammer_masking_information( $sequences, $masking_folder );
            %filter = miRkwood::Utils::merge_hashes_of_arrays(\%filter, \%rrna_mask);
        }
    }
    return %filter;
}

1;
