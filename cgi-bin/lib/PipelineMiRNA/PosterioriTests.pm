package PipelineMiRNA::PosterioriTests;

# ABSTRACT: Everything needed to run the a posteriori tests

use strict;
use warnings;

use File::Spec;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Components;

=method test_randfold

Run the Randfold a posteriori test

=cut

sub test_randfold {
    my ( $candidate_dir, $seq_file ) = @_;
    my $randfold_out = File::Spec->catfile( $candidate_dir, 'randfold.out' );
    PipelineMiRNA::Programs::run_randfold( $seq_file, $randfold_out, 200)
      or die('Problem when running Randfold');
    chmod 777, $randfold_out;
    return -e $randfold_out;
}


=method test_alignment

Run the Alignment (exonerate) a posteriori test

=cut

sub test_alignment {
    my ( $candidate_dir, $candidate_rnafold_stemploop_out ) = @_;

    my $candidate_ct_stemloop_file =
      File::Spec->catfile( $candidate_dir, 'outB2ct_stemloop.ct' );
    PipelineMiRNA::Programs::convert_to_ct( $candidate_rnafold_stemploop_out,
        $candidate_ct_stemloop_file )
      or die('Problem when converting to CT format');

    my $seqN = File::Spec->catfile( $candidate_dir, 'seqWithN.txt' );
#    open( my $SEQN_FH, '>>', $seqN )
#      or die "Impossible d'ouvrir le fichier d'entree  : $!";
    PipelineMiRNA::Components::mask_CT_file($candidate_ct_stemloop_file, $seqN);
    my $exonerate_out = File::Spec->catfile( $candidate_dir, 'alignement.txt' );
    PipelineMiRNA::Programs::run_exonerate( $seqN, $exonerate_out )
      or die('Problem when running Exonerate');
    return $exonerate_out;
}

1;
