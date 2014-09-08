package miRkwood::Pipeline;

# ABSTRACT: Pipeline object

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use miRkwood;
use miRkwood::CandidateHandler;
use miRkwood::FileUtils;
use miRkwood::Maskers;
use miRkwood::Paths;
use miRkwood::SequenceJob;
use miRkwood::Utils;

use Data::Dumper;

### Data ##
my $dirData = miRkwood::Paths->get_data_path();

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my $job_dir = shift @args;
    my $self = {
        job_dir => $job_dir,
        sequences => undef
    };
    bless $self, $class;
    return $self;
}

=method run_pipeline

Run the pipeline

 Usage : miRkwood::MainPipeline::fasta_pipeline( $idirJob );
 Input : The job directory
 Return: -

=cut

sub run_pipeline {
    my ($self, @args) = @_;
    $self->init_pipeline();
    $self->init_sequences();
    $self->run_pipeline_on_sequences();
    return;
}

=method setup_logging

 Usage : $self->setup_logging();
 Return: -

=cut

sub setup_logging {
    my ($self, @args) = @_;
    my $log_file = File::Spec->catfile( $self->{'job_dir'}, 'log.log' );
    $Log::Message::Simple::DEBUG_FH = miRkwood->LOGFH($log_file);
    miRkwood->DEBUG(1);
    return;
}

=method init_pipeline

Initialise the pipeline setup

 Usage : $self->init_pipeline();
 Return: -

=cut

sub init_pipeline {
    my ($self, @args) = @_;
    $self->setup_logging();
    my $run_options_file = $self->get_config_path();
    miRkwood->CONFIG_FILE($run_options_file);
    miRkwood::Programs::init_programs();
    mkdir $self->get_workspace_path();
    mkdir $self->get_candidates_dir();
    return;
}


=method init_sequences

Get the sequences to process from the job directory
This includes parsing, and masking if option selected.

 Usage : my @fasta_array = get_sequences($job_dir);
 Return: -

=cut

sub init_sequences {
    my ($self, @args) = @_;
    debug( "Getting sequences", miRkwood->DEBUG() );
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

=method get_sequences

Accessor to the sequences

 Usage : $self->get_sequences();

=cut

sub get_sequences {
    my ($self, @args) = @_;
    return @{$self->{'sequences'}};
}


=method get_job_dir

Accessor to the job directory

 Usage : $self->get_sequences();

=cut

sub get_job_dir {
    my ($self, @args) = @_;
    return $self->{'job_dir'};
}

=method get_candidates_dir

Return the path to the candidates directory

 Usage : $self->get_candidates_dir();

=cut

sub get_candidates_dir {
    my ($self, @args) = @_;
    my $candidates_dir = File::Spec->catdir( $self->get_job_dir(), 'candidates' );
}

=method get_uploaded_sequences_file

Return the path to the input sequences

 Usage : $self->get_uploaded_sequences_file();

=cut

sub get_uploaded_sequences_file {
    my ($self, @args) = @_;
    return File::Spec->catfile( $self->{'job_dir'}, 'input_sequences.fas' );
}

=method run_pipeline_on_sequences

Run the pipeline on the given sequences

 Usage : $self->run_pipeline_on_sequences();

=cut

sub run_pipeline_on_sequences {
    my ($self, @args) = @_;
    my @sequences_array = $self->get_sequences();
    my $sequences_count = scalar @sequences_array;
    debug( "$sequences_count sequences to process", miRkwood->DEBUG() );
    $self->compute_candidates();
    $self->process_tests();
    debug('miRkwood processing done', miRkwood->DEBUG() );
    $self->mark_job_as_finished();
    debug("Writing finish file", miRkwood->DEBUG() );
    return;
}

sub compute_candidates {
    my ($self, @args) = @_;
    my @sequences_array = $self->get_sequences();
    my $sequence_identifier = 0;
    foreach my $item (@sequences_array) {
        my ( $name, $sequence ) = @{$item};
        debug( "Considering sequence $sequence_identifier: $name",
               miRkwood->DEBUG() );
        $sequence_identifier++;
        my $sequence_dir =
          File::Spec->catdir( $self->get_workspace_path(), $sequence_identifier );
        mkdir $sequence_dir;
        my $sequence_job = miRkwood::SequenceJob->new($sequence_dir, $sequence_identifier, $name, $sequence);
        my $sequence_candidates = $sequence_job->run();
    }
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
    my $masking_folder = File::Spec->catdir($self->get_job_dir(), 'masks');
    mkdir $masking_folder;
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
    return %filter;
}




=method get_job_config_path

Given a job directory, return the path to the job configuration file

=cut

sub get_config_path {
    my ($self, @args) = @_;
    my $job_config_path = File::Spec->catfile( $self->get_job_dir(), 'run_options.cfg' );
    return $job_config_path;
}

=method get_workspace_path

Return the path to the job workspace

=cut

sub get_workspace_path {
    my ($self, @args) = @_;
    return File::Spec->catdir($self->get_job_dir(), 'workspace');
}

=method mark_job_as_finished

Mark the current job as finished

=cut

sub mark_job_as_finished {
    my ($self, @args) = @_;
    my $is_finished_file = File::Spec->catfile( $self->get_job_dir(), 'finished' );
    open( my $finish, '>', $is_finished_file )
        or die "Error when opening $is_finished_file: $!";
    close $finish;
    return (-e $is_finished_file);
}


=method process_tests

Perform the a posteriori tests for a given job

=cut

sub process_tests {
    my ($self, @args) = @_;
    my $workspace_dir = $self->get_workspace_path();
    my @sequence_dirs = miRkwood::FileUtils::get_dirs_from_directory($workspace_dir);
    foreach my $dir (@sequence_dirs)
    {
        my $sequence_dir = File::Spec->catdir( $workspace_dir, $dir );
        debug( "Entering sequence $sequence_dir", miRkwood->DEBUG() );

        my @candidate_dirs = miRkwood::FileUtils::get_dirs_from_directory($sequence_dir);
        foreach my $subDir (@candidate_dirs) {
            my $candidate_dir =
              File::Spec->catdir( $sequence_dir, $subDir );

            if (
                !eval {
                    miRkwood::CandidateHandler
                      ->serialize_candidate_from_run( $self->get_job_dir(), $dir,
                        $subDir, $self->get_candidates_dir() );
                }
              )
            {
                # Catching
                carp( "Serialization of $subDir failed" );
            }
            else {
                debug( "Done with serializing $subDir", miRkwood->DEBUG() );
                # All is well
            }
        } # foreach my $file (@files)
        debug( "Done with initial sequence $dir", miRkwood->DEBUG() );
    } # foreach my $dir (@dirs)
    debug( "Done with all the tests", miRkwood->DEBUG() );
    return 0;
}

1;