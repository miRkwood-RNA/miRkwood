package PipelineMiRNA::PosterioriTests;

use strict;
use warnings;

use File::Spec;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Components;

sub test_mfei {
    my ( $candidate_dir, $candidate_ct_file, $seq ) = @_;
    my $MFEI_output = File::Spec->catfile( $candidate_dir, 'outMFEI.txt' );
    PipelineMiRNA::Components::compute_energy( $candidate_ct_file, $MFEI_output,
        $seq );
}

sub test_randfold {
    my ( $candidate_dir, $seq_file ) = @_;
    my $randfold_out = File::Spec->catfile( $candidate_dir, 'pvalue.txt' );
    PipelineMiRNA::Programs::run_randfold( $seq_file, $randfold_out )
      or die("Problem when running Randfold");
    chmod 777, $randfold_out;
}

sub test_selfcontain {
    my ( $candidate_dir, $seq_file ) = @_;
    my $selfcontain_out =
      File::Spec->catfile( $candidate_dir, 'selfContain.txt' );
    PipelineMiRNA::Programs::run_selfcontain( $seq_file, $selfcontain_out )
      or die("Problem when running Selfcontain");
    chmod 777, $selfcontain_out;
}

sub test_alignment {

    # creation sequence boucle terminale masquee avec des N pour
    # chaque sequence (repertoire ) et resultat alignement mirBASE
    my ( $candidate_dir, $seq_file ) = @_;
    my $seqN = File::Spec->catfile( $candidate_dir, 'seqWithN.txt' );
    open( my $SEQN_FH, '>>', $seqN )
      or die "Impossible d'ouvrir le fichier d'entree  : $!";
    open( my $SEQUENCE_FH, '<', $seq_file )
      or die "Impossible d'ouvrir le fichier d'entree  : $!";
    while ( my $line = <$SEQUENCE_FH> ) {
        if ( $line =~ /^>/ ) {
            print $SEQN_FH $line;
        }
        else {
            $line =~ s/U/T/gi;
            print $SEQN_FH $line;
        }
    }
    close $SEQN_FH;
    close $SEQUENCE_FH;
    my $exonerate_out = File::Spec->catfile( $candidate_dir, 'alignement.txt' );
    PipelineMiRNA::Programs::run_exonerate( $seqN, $exonerate_out )
      or die("Problem when running Exonerate");
}

1;
