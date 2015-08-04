package miRkwood::Pipeline;

# ABSTRACT: Pipeline object

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use miRkwood;
use miRkwood::CandidateHandler;
use miRkwood::FileUtils;
use miRkwood::Paths;
use miRkwood::Utils;
use miRkwood::SequenceJob;
use miRkwood::HairpinBuilder;
use miRkwood::PrecursorBuilder;
use miRkwood::BEDHandler;


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
    mkdir miRkwood::Paths::get_dir_candidates_path_from_job_dir( $self->get_job_dir() );
    $self->create_additional_directories();
    return;
}

=method create_additional_directories

=cut
sub create_additional_directories {
    my ($self, @args) = @_;
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
    select((select($Log::Message::Simple::DEBUG_FH), $| = 1)[0]);
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

=method count_clusters

Count the total number of clusters
(including those without results)

=cut
sub count_clusters {
    my ($self, @args) = @_;
    my $count = 0;
    foreach my $chromosome ( keys%{$self->{'sequences'}} ){
        $count += scalar ( @{$self->{'sequences'}{$chromosome}} );
    }
    return $count;
}

=method get_job_dir

Accessor to the job directory

 Usage : $self->get_job_dir();

=cut

sub get_job_dir {
    my ($self, @args) = @_;
    return $self->{'job_dir'};
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

    my $cfg = miRkwood->CONFIG();
    my $mode = $cfg->param('job.mode');
    $self->{'basic_candidates'} = [];
    my $sequences_count = 0;

    if ( $mode eq 'WebBAM' ){
        $sequences_count = $self->count_clusters();
    }
    else{
        my @sequences_array = $self->get_sequences();
        $sequences_count = scalar @sequences_array;
    }

    debug( "$sequences_count sequences to process", miRkwood->DEBUG() );

    $self->compute_candidates();
    debug('miRkwood processing done', miRkwood->DEBUG() );

    $self->serialize_basic_candidates( 'basic_candidates' );

    $self->mark_job_as_finished();

    debug('Writing finish file', miRkwood->DEBUG() );

    return;
}

=method compute_candidates

=cut
sub compute_candidates {
    my ($self, @args) = @_;

    my $loci = $self->{'sequences'};
    my $sequence_identifier = 0;

    my @sequences_array = $self->get_sequences();
    foreach my $item (@sequences_array) {
        my ( $name, $sequence ) = @{$item};
        debug( "- Considering sequence $sequence_identifier: $name", miRkwood->DEBUG() );
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
        miRkwood::CandidateHandler->serialize_candidate_information( miRkwood::Paths::get_dir_candidates_path_from_job_dir( $self->get_job_dir() ), $candidate );
        push $self->{'basic_candidates'}, $candidate->get_basic_informations();
    }
    return;
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
    return miRkwood::Paths->get_workspace_path( $self->get_job_dir() );
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
    my $type = shift @args;     # $type should be 'basic_candidates' or 'basic_known_candidates'

    my $serialization_file = File::Spec->catfile($self->get_job_dir(), "$type.yml");

    return YAML::XS::DumpFile($serialization_file, $self->{$type});

}

1;
