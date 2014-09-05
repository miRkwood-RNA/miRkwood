package miRkwood::CandidateJob;

# ABSTRACT: Job processing a candidate

use strict;
use warnings;

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
    my ($directory) = @args;
    my $self = {
        directory => $directory,
    };
    bless $self, $class;
    return $self;
}

sub get_directory {
    my ( $self, @args ) = @_;
    return $self->{'directory'};
}


=method process_tests_for_candidate

Perform the a posteriori tests for a given candidate

=cut

sub process_tests_for_candidate {
    my ( $self, @args ) = @_;
    my $dir = $self->get_directory();
    my $seq_file = File::Spec->catfile( $self->get_directory(), 'seq.txt' );
    my $candidate_rnafold_optimal_out =
      File::Spec->catfile( $self->get_directory(), 'outRNAFold_optimal.txt' );
    my $candidate_rnafold_stemploop_out =
      File::Spec->catfile( $self->get_directory(), 'outRNAFold_stemloop.txt' );

    my $cfg = miRkwood->CONFIG();

    ####calcul p-value randfold
    if ( $cfg->param('options.randfold') ) {
        debug( "Running test_randfold on $seq_file", miRkwood->DEBUG() );
        miRkwood::PosterioriTests::test_randfold( $self->get_directory(),
            $seq_file );
    }

    if ( $cfg->param('options.align') ) {
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
