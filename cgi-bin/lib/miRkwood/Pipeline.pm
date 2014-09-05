package miRkwood::Pipeline;

# ABSTRACT: Pipeline object

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use miRkwood;
use miRkwood::CandidateHandler;
use miRkwood::Components;
use miRkwood::FileUtils;
use miRkwood::Maskers;
use miRkwood::MiRdup;
use miRkwood::Paths;
use miRkwood::PosterioriTests;
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
    my $sequence_dir_name = 0;
    foreach my $item (@sequences_array) {
        my ( $name, $sequence ) = @{$item};
        debug( "Considering sequence $sequence_dir_name: $name",
               miRkwood->DEBUG() );
        $sequence_dir_name++;
        my $sequence_dir =
          File::Spec->catdir( $self->get_workspace_path(), $sequence_dir_name );
        mkdir $sequence_dir;

        my $candidates_array =
          $self->get_raw_candidates( $sequence_dir, $name, $sequence );

        my %candidates_hash = $self->process_raw_candidates($candidates_array);

        create_directories( \%candidates_hash, $sequence_dir );
    }
    return;
}

=method get_raw_candidates

Get the candidates for the given 

 Usage : $self->get_raw_candidates($sequence_dir, $name, $sequence);

=cut

sub get_raw_candidates{
    my ($self, @args) = @_;
    my ($sequence_dir, $name, $sequence) = @args;
    my $sequence_job = miRkwood::SequenceJob->new($sequence_dir, $name, $sequence, '+');
    my $candidates = $sequence_job->get_raw_candidates_for_sequence();

    my @candidates_array1 = @{$candidates};
    my @candidates_array;

    my $cfg = miRkwood->CONFIG();
    if ( $cfg->param('options.strands') ) {
        debug( "Processing the other strand", miRkwood->DEBUG() );
        my $reversed_sequence =
          miRkwood::Utils::reverse_complement($self->{'sequence'});
        my $sequence_job2 = miRkwood::SequenceJob->new($sequence_dir, $name, $reversed_sequence, '-');
        my $candidates2 = $sequence_job->compute_candidates_for_sequence();
        my @candidates_array2 = @{$candidates2};
        @candidates_array = sort { $a->{start} <=> $b->{start} } ( @candidates_array1, @candidates_array2 );
    }
    else {
        @candidates_array = sort { $a->{start} <=> $b->{start} } ( @candidates_array1 );
    }
    return \@candidates_array;
}

sub process_raw_candidates{
    my ($self, @args) = @_;
    my $candidates_array = shift @args;

    my @candidates_array = @{$candidates_array};
    my $cfg = miRkwood->CONFIG();
    if ( $cfg->param('options.mfe') ) {
        debug('Select only sequences with MFEI < -0.6', miRkwood->DEBUG() );
        @candidates_array = grep { mfei_below_threshold($_, -0.6) } @candidates_array;
    }

    my %candidates_hash;
    if (@candidates_array) {
        debug("Merging candidates", miRkwood->DEBUG() );
        %candidates_hash = $self->merge_candidates( \@candidates_array );
    }
    else {
        %candidates_hash = ();
    }
    return %candidates_hash;
}

=method mfei_below_threshold

Return whether a given candidate has its mfei above a given threshold

=cut

sub mfei_below_threshold {
    my @args      = @_;
    my $candidate = shift @args;
    my $threshold = shift @args;
    my %current_candidate = %{ $candidate };
    my $mfei = $current_candidate{'mfei'};
    return $mfei < $threshold;
}

=method merge_candidates

Process the candidates and try merging them.
We assume the candidates array is already sorted by growing position

=cut

sub merge_candidates {
    my ($self, @args) = @_;
    my $candidates_array = shift @args;
    my @candidates_array = @{$candidates_array};
    my $nb_candidates = scalar @candidates_array;

    my @merged_candidates   = ();
    my %final_hash          = ();
    my %reference_candidate = %{ $candidates_array[0] };
    my %best_candidate      = %reference_candidate;
    my %current_candidate;
    for my $candidate_index ( 1 .. $#candidates_array ) {
        %current_candidate = %{ $candidates_array[$candidate_index] };
        my $start = $current_candidate{'start'};
        my $end   = $current_candidate{'end'};
        my ( $ref_start, $ref_end ) =
          ( $reference_candidate{'start'}, $reference_candidate{'end'} );
        if ( is_included( $start, $end, $ref_start, $ref_end ) ) {
            if ( $best_candidate{'mfei'} <= -0.8 ) {
                push @merged_candidates, {%current_candidate};
            }
            else {
                if ( $current_candidate{'mfei'} < $best_candidate{'mfei'} ) {
                    push @merged_candidates, {%best_candidate};
                    %best_candidate = %current_candidate;
                }
                else {
                    push @merged_candidates, {%current_candidate};
                }
            }

        }
        elsif ( is_overlapping( $start, $end, $ref_start, $ref_end ) ) {
            if ( $current_candidate{'mfei'} < $best_candidate{'mfei'} ) {
                push @merged_candidates, {%best_candidate};
                %best_candidate = %current_candidate;
            }
            else {
                push @merged_candidates, {%current_candidate};
            }
        }
        else {
            my $final_name = $best_candidate{'name'};
            $final_hash{$final_name}                 = {};
            $final_hash{$final_name}{'max'}          = {%best_candidate};
            $final_hash{$final_name}{'alternatives'} = [@merged_candidates];
            @merged_candidates                       = undef;
            @merged_candidates                       = ();
            %reference_candidate                     = %current_candidate;
            %best_candidate                          = %reference_candidate;
        }

    }    #foreach
    my $final_name = $best_candidate{'name'};
    $final_hash{$final_name}                 = {};
    $final_hash{$final_name}{'max'}          = {%best_candidate};
    $final_hash{$final_name}{'alternatives'} = [@merged_candidates];
    return %final_hash;
}


=method is_overlapping

Test whether one candidate is overlapping with the other
(based on their positions).

=cut

sub is_overlapping {
    my @args      = @_;
    my $start     = shift @args;
    my $end       = shift @args;
    my $ref_start = shift @args;
    my $ref_end   = shift @args;
    ( $start ne '' && $end ne '' && $ref_start ne '' && $ref_end ne '' ) or die('Not enough values provided');
    $ref_start <= $start or die("Positions should be ordered : $ref_start <= $start");
    return ( $start < ( $ref_start + $ref_end ) / 2 );
}

=method is_included

Test whether one candidate is included into the other
(based on their positions).

=cut

sub is_included {
    my @args      = @_;
    my $start     = shift @args;
    my $end       = shift @args;
    my $ref_start = shift @args;
    my $ref_end   = shift @args;
    ( $start ne '' && $end ne '' && $ref_start ne '' && $ref_end ne '' ) or die('Not enough values provided');
    $ref_start <= $start or die("Positions should be ordered : $ref_start <= $start");
    return ( $end <= $ref_end );
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


=method create_directories

Create the necessary directories.

=cut

sub create_directories {
    my @args                 = @_;
    my (%candidates_hash)    = %{ shift @args };
    my $current_sequence_dir = shift @args;
    my $candidate_counter    = 0;
    foreach my $key ( sort keys %candidates_hash ) {
        $candidate_counter++;
        my $candidate_dir =
          File::Spec->catdir( $current_sequence_dir, $candidate_counter );
        mkdir $candidate_dir;
        my $candidate_ref = $candidates_hash{$key}{'max'};
        populate_candidate_directory( $candidate_dir, $candidate_ref,
            $candidates_hash{$key}{'alternatives'} );
    }
    return;
}

=method populate_candidate_directory

Populate a candidate directory with the sequence, strand & so on.

=cut

sub populate_candidate_directory {
    my @args               = @_;
    my $candidate_dir      = shift @args;
    my %candidate          = %{ shift @args };
    my @alternatives_array = @{ shift @args };

    debug( "Writing candidate information in $candidate_dir", miRkwood->DEBUG() );
    #Writing seq.txt
    my $candidate_sequence = File::Spec->catfile( $candidate_dir, 'seq.txt' );
    open( my $SEQ_FH, '>', $candidate_sequence )
      or die "Error when opening $candidate_sequence: $!";
    print $SEQ_FH ">$candidate{'name'}\n$candidate{'dna'}\n";
    close $SEQ_FH;

    #Writing sequence information
    my $seq_info_file = File::Spec->catfile( $candidate_dir, 'sequence_information.txt' );
    open( my $SEQ_INFO_FH, '>', $seq_info_file )
      or die "Error when opening $seq_info_file: $!";
    print $SEQ_INFO_FH $candidate{'strand'} . "\t" .  $candidate{'start'} . "\t" . $candidate{'end'};
    close $SEQ_INFO_FH;

    process_outRNAFold( $candidate_dir, 'optimal', $candidate{'name'},
        $candidate{'dna'}, $candidate{'structure_optimal'}, $candidate{'energy_optimal'} );
    process_outRNAFold( $candidate_dir, 'stemloop', $candidate{'name'},
        $candidate{'dna'}, $candidate{'structure_stemloop'}, $candidate{'energy_stemloop'} );

    # Writing energy file
    my $energy_file = File::Spec->catfile( $candidate_dir, 'outMFEI.txt' );
    open( my $ENERGY_FH, '>', $energy_file )
        or die "Unable to open $energy_file: $!";
    my $content = $candidate{'name'} . "\t" . $candidate{'mfei'} . "\t" . $candidate{'energy_optimal'} . "\t" . $candidate{'amfe'};
    print $ENERGY_FH $content;
    close $ENERGY_FH or die "Unable to close: $!";

    #Writing alternativeCandidates.txt
    my $alternative_candidates =
      File::Spec->catfile( $candidate_dir, 'alternativeCandidates.txt' );
    if (@alternatives_array){
        open( my $OUT2, '>>', $alternative_candidates )
          or die "Error when opening $alternative_candidates: $!";
        foreach my $alternative (@alternatives_array) {
            print $OUT2
">$alternative->{'name'}\t$alternative->{'dna'}\t$alternative->{'structure_optimal'}\t$alternative->{'mfei'}\n";
        }
        close $OUT2;
    }

    # Writing VARNA image
    my $cfg = miRkwood->CONFIG();

    if ( $cfg->param('options.varna') ) {
        my $varna_image = File::Spec->catfile( $candidate_dir, 'image.png' );
        debug( "Generating image using VARNA in $varna_image", miRkwood->DEBUG() );
        miRkwood::Programs::run_varna_on_structure( $candidate{'dna'}, $candidate{'structure_stemloop'}, $varna_image )
          or carp('Problem during image generation using VARNA');
    }

    return;
}

=method process_outRNAFold

Writing (pseudo) rnafold output

=cut

sub process_outRNAFold {
    my ( $candidate_dir, $suffix, $nameSeq, $dna, $structure, $energy ) = @_;

    my $candidate_rnafold_output =
      File::Spec->catfile( $candidate_dir, "outRNAFold_$suffix.txt" );

    open( my $OUT2, '>', $candidate_rnafold_output )
      or die "Error when opening $candidate_rnafold_output: $!";
    print $OUT2 ">$nameSeq\n$dna\n$structure ($energy)\n";
    close $OUT2;

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
    debug( "A posteriori tests in $self->get_job_dir()", miRkwood->DEBUG() );
    
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

            debug( "Entering candidate $subDir", miRkwood->DEBUG() );
            process_tests_for_candidate( $candidate_dir, $subDir );
            debug( "Done with candidate $subDir", miRkwood->DEBUG() );
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

=method process_tests_for_candidate

Perform the a posteriori tests for a given candidate

=cut

sub process_tests_for_candidate {

    my @args = @_;
    my ( $candidate_dir, $file ) = @args;

    ####Traitement fichier de sortie outStemloop
    chmod 0777, $candidate_dir;

    my $seq_file = File::Spec->catfile( $candidate_dir, 'seq.txt' );
    my $candidate_rnafold_optimal_out =
      File::Spec->catfile( $candidate_dir, 'outRNAFold_optimal.txt' );
    my $candidate_rnafold_stemploop_out =
      File::Spec->catfile( $candidate_dir, 'outRNAFold_stemloop.txt' );

    my $cfg = miRkwood->CONFIG();

    ####calcul p-value randfold
    if ( $cfg->param('options.randfold') ) {
        debug( "Running test_randfold on $seq_file", miRkwood->DEBUG() );
        miRkwood::PosterioriTests::test_randfold( $candidate_dir,
            $seq_file );
    }

    if ( $cfg->param('options.align') ) {
        debug( "Running test_alignment on $candidate_rnafold_stemploop_out", miRkwood->DEBUG() );
        my $file_alignement =
          miRkwood::PosterioriTests::test_alignment( $candidate_dir,
            $candidate_rnafold_stemploop_out );
        post_process_alignments( $candidate_dir,
            $candidate_rnafold_stemploop_out,
            $file_alignement );
    }

    return;
}

=method post_process_alignments


=cut

sub post_process_alignments {
    my @args                            = @_;
    my $candidate_dir                   = shift @args;
    my $candidate_rnafold_stemploop_out = shift @args;
    my $file_alignement                 = shift @args;

    my @res =
      miRkwood::Components::get_data_from_rnafold_out(
        $candidate_rnafold_stemploop_out);
    my ( $name, $position, $DNASequence, $Vienna ) = @res;
    my %alignments;

    if (-z $file_alignement){
        return;
    }
    if (
        !eval {
            %alignments =
              miRkwood::Components::parse_custom_exonerate_output(
                $file_alignement);
        }
      )
    {
        # Catching exception
        carp("Exception when parsing exonerate output $file_alignement");
        return;
    }
    else {
        %alignments =
          miRkwood::Components::merge_alignments( \%alignments );
        my $tmp_file =
          File::Spec->catfile( $candidate_dir, "mirdup_validation.txt" );
        my %mirdup_results =
          miRkwood::MiRdup->validate_with_mirdup( $tmp_file, $name,
            $DNASequence, $Vienna, keys %alignments );
        my $mirdup_results_file =
          File::Spec->catfile( $candidate_dir, 'mirdup_results.yml' );
        YAML::XS::DumpFile( $mirdup_results_file, %mirdup_results );

        my $alignments_results_file =
          File::Spec->catfile( $candidate_dir, 'merged_alignments.yml' );
        YAML::XS::DumpFile( $alignments_results_file, %alignments );
    }
}

1;