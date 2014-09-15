package miRkwood::PosterioriTests;

# ABSTRACT: Everything needed to run the a posteriori tests

use strict;
use warnings;

use miRkwood::Parsers;
use miRkwood::Programs;
use miRkwood::Components;

use File::Spec;
use Carp;

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

=method test_randfold

Run the Randfold a posteriori test

=cut

sub test_randfold {
    my ( $self, @args ) = @_;
    my $seq_file = shift @args;
    my $randfold_out = File::Spec->catfile( $self->get_directory(), 'randfold.out' );
    miRkwood::Programs::run_randfold( $seq_file, $randfold_out, 200)
      or die('Problem when running Randfold');
    my $res = miRkwood::Parsers::parse_pvalue($randfold_out);
    return $res;
}


=method test_alignment

Run the Alignment (exonerate) a posteriori test

=cut

sub test_alignment {
    my ( $self, @args ) = @_;
    my $candidate_rnafold_stemploop_out = shift @args;

    my $candidate_ct_stemloop_file =
      File::Spec->catfile( $self->get_directory(), 'outB2ct_stemloop.ct' );
    miRkwood::Programs::convert_to_ct( $candidate_rnafold_stemploop_out,
        $candidate_ct_stemloop_file )
      or die('Problem when converting to CT format');

    my $seqN = File::Spec->catfile( $self->get_directory(), 'seqWithN.txt' );
    miRkwood::Components::mask_CT_file($candidate_ct_stemloop_file, $seqN);
    my $exonerate_out = File::Spec->catfile( $self->get_directory(), 'alignement.txt' );
    miRkwood::Programs::run_exonerate( $seqN, $exonerate_out )
      or die('Problem when running Exonerate');
    my ($mirdup_results, $alignments) = $self->post_process_alignments($exonerate_out );
    return ($mirdup_results, $alignments);
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
              miRkwood::Parsers::parse_custom_exonerate_output(
                $file_alignement);
        }
      )
    {
        # Catching exception
        carp("Exception when parsing exonerate output $file_alignement");
        return;
    }
    else {
        %alignments = $self->merge_alignments( \%alignments );
        my $tmp_file =
          File::Spec->catfile( $self->get_directory(), "mirdup_validation.txt" );
        my %mirdup_results =
          miRkwood::MiRdup->validate_with_mirdup( $tmp_file, $self->{'sequence_name'},
            $self->{'candidate'}{'sequence'}, $self->{'candidate'}{'structure_stemloop'}, keys %alignments );
        return (\%mirdup_results, \%alignments);
    }
}


=method merge_alignments

Merge overlapping alignments.
Given ordered positions, merge in [a..b] all [c..d] if a<=c and d<=b+2

=cut

sub merge_alignments {
    my ( $self, @args ) = @_;
    my $alignments = shift @args;
    my %alignments = %{$alignments};

    my %merged_alignments;
    my ( $stocked_left, $stocked_right ) = ( -10, -10 );

    my @keys = sort {
        ( miRkwood::Utils::get_element_of_split( $a, '-', 0 )
              <=> miRkwood::Utils::get_element_of_split( $b, '-', 0 ) )
          || ( miRkwood::Utils::get_element_of_split( $a, '-', 1 )
            <=> miRkwood::Utils::get_element_of_split( $b, '-', 1 ) )
    } keys %alignments;
    my @stocked_hits;
    my $final_key;
    my $final_hit_count = -1;

    foreach my $current_key (@keys) {
        my ( $current_left, $current_right ) = split( /-/, $current_key );
        my $current_hit_count = scalar @{ $alignments{$current_key} };

        if ( $current_right > $stocked_right + 3 ) {

 # No merge ; drop the list of current hits in the hash (only if there are some)
            if (@stocked_hits) {
                push @{ $merged_alignments{$final_key} }, @stocked_hits;
            }

            # Reinitialise
            $final_hit_count = -1;
            @stocked_hits    = ();
            ( $stocked_left, $stocked_right ) =
              ( $current_left, $current_right );
        }
        if ( $current_hit_count > $final_hit_count ) {

    # This position holds more hits than the previous, so it will be our new key
            $final_hit_count = $current_hit_count;
            $final_key       = $current_key;
        }

        # Stock the current hits
        push @stocked_hits, @{ $alignments{$current_key} };
    }

    # Drop the remaining hits in the hash
    push @{ $merged_alignments{$final_key} }, @stocked_hits;
    return %merged_alignments;
}



1;
