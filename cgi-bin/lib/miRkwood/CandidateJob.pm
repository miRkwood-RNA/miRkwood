package miRkwood::CandidateJob;

# ABSTRACT: Job processing a candidate

use strict;
use warnings;

use miRkwood;
use miRkwood::CandidateHandler;
use miRkwood::Components;
use miRkwood::MiRdup;
use miRkwood::PosterioriTests;
use miRkwood::Parsers;
use Carp;

use YAML::XS;
use Log::Message::Simple qw[msg error debug];

=method new

Constructor

Usage:
  my $candidate_job = miRkwood::CandidateJob->new($candidate_dir);

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($directory, $sequence_name, $identifier, $candidate_ref, $alternatives) = @args;
    my $self = {
        sequence_name => $sequence_name,
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
    my $candidate_test_info = $self->process_tests_for_candidate();
    my $candidate_information = $self->get_candidate_information();

    my %complete_candidate = (%{$candidate_test_info}, %{$candidate_information});
    my $candidate = miRkwood::Candidate->new(\%complete_candidate);
    $candidate->compute_alignment_quality();
    $candidate->compute_quality();
    return $candidate;
}

=method get_candidate_information


=cut

sub get_candidate_information {
    my ( $self, @args ) = @_;
    my $candidate_structure = $self->{'candidate'};

    $candidate_structure = $self->update_candidate_information_from_self($candidate_structure);
    $candidate_structure = $self->update_candidate_information_from_run($candidate_structure);

    return $candidate_structure;
}

sub update_candidate_information_from_self {
    my ( $self, @args ) = @_;
    my $candidate_structure = shift @args;
    $candidate_structure->{'name'} = $self->{'sequence_name'};
    $candidate_structure->{'identifier'} = $self->{'identifier'};
    $candidate_structure->{'alternatives'} = $self->convert_alternative_candidates();
    return $candidate_structure;
}

sub update_candidate_information_from_run{
    my ( $self, @args ) = @_;
    my $candidate_structure = shift @args;
    $candidate_structure->{'mfe'} = delete $candidate_structure->{'energy_stemloop'};
    $candidate_structure->{'image'} = $self->write_VARNA_if_needed();

    $candidate_structure->{'position'} = "$candidate_structure->{'start_position'}-$candidate_structure->{'end_position'}";
    $candidate_structure->{'length'} = $candidate_structure->{'end_position'} - $candidate_structure->{'start_position'} + 1;
    $candidate_structure->{'%GC'} = miRkwood::Utils::restrict_num_decimal_digits(
                            miRkwood::Utils::compute_gc_content($candidate_structure->{'sequence'}),
                            3);

    my $hairpin = miRkwood::Utils::make_ASCII_viz($candidate_structure->{'sequence'}, $candidate_structure->{'structure_stemloop'});
    $candidate_structure->{'hairpin'} = $hairpin;
    return $candidate_structure;
}

=method write_VARNA_if_needed


=cut

sub write_VARNA_if_needed {
    my ( $self, @args ) = @_;
    my $cfg = miRkwood->CONFIG();
    if ( $cfg->param('options.varna') ) {
        my $varna_image = File::Spec->catfile( $self->get_directory(), 'image.png' );
        debug( "Generating image using VARNA in $varna_image", miRkwood->DEBUG() );
        miRkwood::Programs::run_varna_on_structure( $self->{'candidate'}{'sequence'}, $self->{'candidate'}{'structure_stemloop'}, $varna_image )
          or carp('Problem during image generation using VARNA');
      return $varna_image
    }
    return '';
}

=method process_tests_for_candidate

Perform the a posteriori tests for a given candidate

=cut

sub process_tests_for_candidate {
    my ( $self, @args ) = @_;
    my $dir = $self->get_directory();

    my $cfg = miRkwood->CONFIG();
    my $result = {};

    if ( $cfg->param('options.randfold') ) {
        my $seq_file = $self->write_sequence_file();
        debug( "Running test_randfold on $seq_file", miRkwood->DEBUG() );
        $result->{'shuffles'} = miRkwood::PosterioriTests::test_randfold( $self->get_directory(),
            $seq_file );
    }

    if ( $cfg->param('options.align') ) {
        my $candidate_rnafold_stemploop_out = $self->write_RNAFold_stemloop_output();
        debug( "Running test_alignment on $candidate_rnafold_stemploop_out", miRkwood->DEBUG() );
        my $file_alignement =
          miRkwood::PosterioriTests::test_alignment( $self->get_directory(),
            $candidate_rnafold_stemploop_out );
        my ($mirdup_results, $alignments) = $self->post_process_alignments($file_alignement );
        if ($alignments) {
            $result->{'alignment_existence'} = 1;
        }
        if ($result->{'alignment_existence'}){
            $result->{'alignments'} = $alignments;
            $result->{'mirdup_validation'} = $mirdup_results;
        }
    } else {
        $result->{'alignment_existence'} = 0;
    }
    return $result;
}

=method write_RNAFold_stemloop_output

Writing (pseudo) rnafold output

=cut

sub write_RNAFold_stemloop_output {
    my ( $self, @args ) = @_;
    my $candidate_rnafold_output =
      File::Spec->catfile( $self->get_directory(), "outRNAFold_stemloop.txt" );
    open( my $OUT2, '>', $candidate_rnafold_output )
      or die "Error when opening $candidate_rnafold_output: $!";
    print {$OUT2} ">$self->{'candidate'}{'name'}\n$self->{'candidate'}{'sequence'}\n$self->{'candidate'}{'structure_stemloop'} ($self->{'candidate'}{'energy_stemloop'})\n";
    close $OUT2;
    return $candidate_rnafold_output;
}

=method write_sequence_file

Writing the candidate as a FASTA sequence

=cut

sub write_sequence_file {
    my ( $self, @args ) = @_;
    my $candidate_sequence = File::Spec->catfile( $self->get_directory(), 'seq.txt' );
    open( my $SEQ_FH, '>', $candidate_sequence )
        or die "Error when opening $candidate_sequence: $!";
     print {$SEQ_FH} ">$self->{'candidate'}{'name'}\n$self->{'candidate'}{'sequence'}\n";
    close $SEQ_FH;
    return $candidate_sequence;
}

=method post_process_alignments


=cut

sub post_process_alignments {
    my ( $self, @args ) = @_;
    my $file_alignement = shift @args;

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
        %alignments = miRkwood::Components::merge_alignments( \%alignments );
        my $tmp_file =
          File::Spec->catfile( $self->get_directory(), "mirdup_validation.txt" );
        my %mirdup_results =
          miRkwood::MiRdup->validate_with_mirdup( $tmp_file, $self->{'sequence_name'},
            $self->{'candidate'}{'sequence'}, $self->{'candidate'}{'structure_stemloop'}, keys %alignments );
        return (\%mirdup_results, \%alignments);
    }
}

=method convert_alternative_candidates

=cut

sub convert_alternative_candidates {
    my ( $self, @args ) = @_;
    my @alternatives = @{$self->{'alternatives'}};
    my %results;
    foreach my $alternative (@alternatives) {
        my $name = delete $alternative->{'name'};
        $results{$name} = $alternative;
    }
    return \%results;
}

1;
