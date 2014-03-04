package PipelineMiRNA::PosterioriTests;

# ABSTRACT: Everything needed to run the a posteriori tests

use strict;
use warnings;

use File::Spec;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Components;

=method test_mfei

Run the MFEI a posteriori test

=cut

sub test_mfei {
    my ( $candidate_dir, $candidate_ct_file, $seq ) = @_;
    my $MFEI_output = File::Spec->catfile( $candidate_dir, 'outMFEI.txt' );
    PipelineMiRNA::Components::compute_energy( $candidate_ct_file, $MFEI_output,
        $seq );
    return -e $MFEI_output;
}

=method test_randfold

Run the Randfold a posteriori test

=cut

sub test_randfold {
    my ( $candidate_dir, $seq_file ) = @_;
    my $randfold_out = File::Spec->catfile( $candidate_dir, 'randfold.out' );
    PipelineMiRNA::Programs::run_randfold( $seq_file, $randfold_out, 200)
      or die("Problem when running Randfold");
    chmod 777, $randfold_out;
    return -e $randfold_out;
}


=method test_alignment

Run the Alignment (exonerate) a posteriori test

=cut

sub test_alignment {

    # creation sequence boucle terminale masquee avec des N pour
    # chaque sequence (repertoire ) et resultat alignement mirBASE
    my ( $candidate_dir, $CT_file ) = @_;
    my $seqN = File::Spec->catfile( $candidate_dir, 'seqWithN.txt' );
#    open( my $SEQN_FH, '>>', $seqN )
#      or die "Impossible d'ouvrir le fichier d'entree  : $!";
    PipelineMiRNA::Components::mask_CT_file($CT_file, $seqN);
    my $exonerate_out = File::Spec->catfile( $candidate_dir, 'alignement.txt' );
    PipelineMiRNA::Programs::run_exonerate( $seqN, $exonerate_out )
      or die("Problem when running Exonerate");
    return $exonerate_out;
}

1;
