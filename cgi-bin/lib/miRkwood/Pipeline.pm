package miRkwood::Pipeline;

# ABSTRACT: Pipeline object

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use miRkwood;
use miRkwood::CandidateHandler;
use miRkwood::FileUtils;
use miRkwood::Paths;
use miRkwood::SequenceJob;
use miRkwood::Utils;


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
    if ( exists($self->{'bam_file'}) or exists($self->{'bed_file'})) {
        mkdir $self->get_reads_dir();
        mkdir $self->get_mirbase_reads_dir();
        mkdir $self->get_new_reads_dir();
    }
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

=method init_sequences

Abstract method.

=cut

sub init_sequences {
    my ($self, @args) = @_;
    die ('Unimplemented method init_sequences');
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

=method get_reads_dir

Return the path to the reads directory

 Usage : $self->get_reads_dir();

=cut

sub get_reads_dir {
    my ($self, @args) = @_;
    my $reads_dir = File::Spec->catdir( $self->get_job_dir(), 'reads' );
}

=method get_mirbase_reads_dir

Return the path to the directory
with reads corresponding to known miRNAs.

 Usage : $self->get_mirbase_reads_dir();

=cut

sub get_mirbase_reads_dir {
    my ($self, @args) = @_;
    my $mirbase_reads_dir = File::Spec->catdir( $self->get_reads_dir(), 'known' );
}

=method get_new_reads_dir

Return the path to the directory
with reads corresponding to new miRNAs.

 Usage : $self->get_new_reads_dir();

=cut

sub get_new_reads_dir {
    my ($self, @args) = @_;
    my $new_reads_dir = File::Spec->catdir( $self->get_reads_dir(), 'new' );
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
    $self->{'basic_candidates'} = [];
    my $sequences_count = scalar @sequences_array;
    debug( "$sequences_count sequences to process", miRkwood->DEBUG() );
    $self->compute_candidates();
    debug('miRkwood processing done', miRkwood->DEBUG() );
    $self->serialize_basic_candidates();
    $self->mark_job_as_finished();
    debug("Writing finish file", miRkwood->DEBUG() );
    return;
}

=method compute_candidates

=cut

sub compute_candidates {
    my ($self, @args) = @_;
    my @sequences_array = $self->get_sequences();
    my $sequence_identifier = 0;
    foreach my $item (@sequences_array) {
        my ( $name, $sequence ) = @{$item};
        debug( "Considering sequence $sequence_identifier: $name", miRkwood->DEBUG() );
        $sequence_identifier++;
        my $sequence_dir = $self->make_sequence_workspace_directory($sequence_identifier);
        my $sequence_job = miRkwood::SequenceJob->new($sequence_dir, $sequence_identifier, $name, $sequence);
        my $sequence_candidates = $sequence_job->run();
        $self->serialize_candidates($sequence_candidates);
    }
    return;
}

=method make_sequence_workspace_directory

Create the directory for the sequence (in the job workspace)
based on its identifier and returns the name

 Usage : my $dir = $self->make_sequence_workspace_directory();

=cut

sub make_sequence_workspace_directory {
    my ($self, @args) = @_;
    my $sequence_identifier = shift @args;
    my $sequence_dir = File::Spec->catdir( $self->get_workspace_path(), $sequence_identifier );
    mkdir $sequence_dir;
    return $sequence_dir;
}

sub serialize_candidates {
    my ($self, @args) = @_;
    my $candidates = shift @args;
    my @candidates_array = @{$candidates};
    foreach my $candidate (@candidates_array ) {
        if ( exists($self->{'bam_file'}) ){ # only for the standalone transcriptome version
            $candidate = $candidate->get_reads($self->{'bam_file'});
            miRkwood::CandidateHandler->print_reads_cloud( $self->get_reads_dir(), $self->{'genome_file'}, $candidate );
        }
        miRkwood::CandidateHandler->serialize_candidate_information( $self->get_candidates_dir(), $candidate );
        
        push $self->{'basic_candidates'}, $candidate->get_basic_informations();
    }
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

sub serialize_basic_candidates {
	
    my ($self, @args) = @_;
    my $serialization_file = File::Spec->catfile($self->get_job_dir(), "basic_candidates.yml");
    	
    return YAML::XS::DumpFile($serialization_file, $self->{'basic_candidates'});
}

1;
