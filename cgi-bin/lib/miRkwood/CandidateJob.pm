package miRkwood::CandidateJob;

# ABSTRACT: Job processing a candidate

use strict;
use warnings;

use miRkwood;
use miRkwood::CandidateHandler;
use miRkwood::Components;
use miRkwood::MiRdup;
use miRkwood::PosterioriTests;

use Log::Message::Simple qw[msg error debug];

=method new

Constructor

Usage:
  my $candidate_job = miRkwood::CandidateJob->new($candidate_dir);

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($directory, $identifier, $candidate_ref, $alternatives) = @args;
    my $self = {
        directory => $directory,
        identifier => $identifier,
        candidate => $candidate_ref,
        alternatives => $alternatives,
    };
    bless $self, $class;
    return $self;
}

sub get_directory {
    my ( $self, @args ) = @_;
    return $self->{'directory'};
}

sub run{
    my ( $self, @args ) = @_;
    $self->populate_candidate_directory();
    $self->process_tests_for_candidate();
    my $candidate_object = miRkwood::CandidateHandler->get_candidate_information_from_run($self->get_directory());
    $candidate_object->{'identifier'} = $self->{'identifier'};
    return $candidate_object;
}

=method populate_candidate_directory

Populate a candidate directory with the sequence, strand & so on.

=cut

sub populate_candidate_directory {
    my ( $self, @args ) = @_;
    my %candidate          = %{ $self->{'candidate'} };
    my @alternatives_array = @{ $self->{'alternatives'} };

    #Writing seq.txt
    my $candidate_sequence = File::Spec->catfile( $self->get_directory(), 'seq.txt' );
    open( my $SEQ_FH, '>', $candidate_sequence )
      or die "Error when opening $candidate_sequence: $!";
    print {$SEQ_FH} ">$candidate{'name'}\n$candidate{'dna'}\n";
    close $SEQ_FH;

    #Writing sequence information
    my $seq_info_file = File::Spec->catfile( $self->get_directory(), 'sequence_information.txt' );
    open( my $SEQ_INFO_FH, '>', $seq_info_file )
      or die "Error when opening $seq_info_file: $!";
    print {$SEQ_INFO_FH} $candidate{'strand'} . "\t" .  $candidate{'start'} . "\t" . $candidate{'end'};
    close $SEQ_INFO_FH;

    $self->process_outRNAFold( 'optimal', $candidate{'name'},
        $candidate{'dna'}, $candidate{'structure_optimal'}, $candidate{'energy_optimal'} );
    $self->process_outRNAFold( 'stemloop', $candidate{'name'},
        $candidate{'dna'}, $candidate{'structure_stemloop'}, $candidate{'energy_stemloop'} );

    # Writing energy file
    my $energy_file = File::Spec->catfile( $self->get_directory(), 'outMFEI.txt' );
    open( my $ENERGY_FH, '>', $energy_file )
        or die "Unable to open $energy_file: $!";
    my $content = $candidate{'name'} . "\t" . $candidate{'mfei'} . "\t" . $candidate{'energy_optimal'} . "\t" . $candidate{'amfe'};
    print $ENERGY_FH $content;
    close $ENERGY_FH or die "Unable to close: $!";

    #Writing alternativeCandidates.txt
    my $alternative_candidates =
      File::Spec->catfile( $self->get_directory(), 'alternativeCandidates.txt' );
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
        my $varna_image = File::Spec->catfile( $self->get_directory(), 'image.png' );
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
    my ( $self, @args ) = @_;
    my ( $suffix, $nameSeq, $dna, $structure, $energy ) = @args;

    my $candidate_rnafold_output =
      File::Spec->catfile( $self->get_directory(), "outRNAFold_$suffix.txt" );

    open( my $OUT2, '>', $candidate_rnafold_output )
      or die "Error when opening $candidate_rnafold_output: $!";
    print {$OUT2} ">$nameSeq\n$dna\n$structure ($energy)\n";
    close $OUT2;
    return;
}


=method process_tests_for_candidate

Perform the a posteriori tests for a given candidate

=cut

sub process_tests_for_candidate {
    my ( $self, @args ) = @_;
    my $dir = $self->get_directory();

    my $cfg = miRkwood->CONFIG();

    ####calcul p-value randfold
    if ( $cfg->param('options.randfold') ) {
        my $seq_file = File::Spec->catfile( $self->get_directory(), 'seq.txt' );
        debug( "Running test_randfold on $seq_file", miRkwood->DEBUG() );
        miRkwood::PosterioriTests::test_randfold( $self->get_directory(),
            $seq_file );
    }

    if ( $cfg->param('options.align') ) {
        my $candidate_rnafold_stemploop_out =
            File::Spec->catfile( $self->get_directory(), 'outRNAFold_stemloop.txt' );
        debug( "Running test_alignment on $candidate_rnafold_stemploop_out", miRkwood->DEBUG() );
        my $file_alignement =
          miRkwood::PosterioriTests::test_alignment( $self->get_directory(),
            $candidate_rnafold_stemploop_out );
        $self->post_process_alignments($candidate_rnafold_stemploop_out, $file_alignement );
    }

    return;
}


=method post_process_alignments


=cut

sub post_process_alignments {
    my ( $self, @args ) = @_;
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
          File::Spec->catfile( $self->get_directory(), "mirdup_validation.txt" );
        my %mirdup_results =
          miRkwood::MiRdup->validate_with_mirdup( $tmp_file, $name,
            $DNASequence, $Vienna, keys %alignments );
        my $mirdup_results_file =
          File::Spec->catfile( $self->get_directory(), 'mirdup_results.yml' );
        YAML::XS::DumpFile( $mirdup_results_file, %mirdup_results );

        my $alignments_results_file =
          File::Spec->catfile( $self->get_directory(), 'merged_alignments.yml' );
        YAML::XS::DumpFile( $alignments_results_file, %alignments );
    }
}

1;
