######################################################################
#  @author : Mohcen BENMOUNAH                                        #
#                                                                    #
#                                                                    #
######################################################################
package PipelineMiRNA::Components;

use strict;
use warnings;

use File::Spec;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Utils;

sub filter_CDS {
    my ( $dirData, $dirJob, $plant ) = @_;

    my $uploaded_sequences =
      File::Spec->catfile( $dirJob, 'sequenceUpload.fas' );
    my $input_sequences = File::Spec->catfile( $dirJob, 'Sequences.fas' );

    my $blast_database = File::Spec->catfile( $dirData, "$plant.fas" );
    my $blast_output   = File::Spec->catfile( $dirJob,  'outBlast.txt' );
    my $blastx_options = '-outfmt 6 -max_target_seqs 1 -evalue 1E-5';
    PipelineMiRNA::Programs::run_blast(
                                        $uploaded_sequences, $blast_database,
                                        $blastx_options,     $blast_output
    ) or die('Problem when running Blastx');
    chmod 777, $blast_output;

    my %blast_seqs;

    # Looping in blast_output, indexing sequences found
    open( my $FOut, '<', $blast_output )
      || die "Problème à l\'ouverture : $!";
    while ( my $line = <$FOut> ) {
        my @name = split( '\t', $line );
        $blast_seqs{ $name[0] } = 1;
    }
    close $FOut || die "Problème à la fermeture : $!";

    # Looping in uploaded_sequences, copying only the files not found in Blast
    open( my $FSeq, '<', $uploaded_sequences )
      or die "Problème à l\'ouverture : $!";
    open( my $RES_FH, '>>', $input_sequences )
      or die "Problème à l\'ouverture : $!";
    PipelineMiRNA::Utils::filter_fasta( $FSeq, $RES_FH, \%blast_seqs );
    close $FSeq   or die "Problème à la fermeture : $!";
    close $RES_FH or die "Problème à la fermeture : $!";
    chmod 777, $input_sequences;

    return $input_sequences;
}

sub compute_energy {
    my ( $CT, $MFEI_output, $sequence_name ) = @_;

    # ouverture du fichier de sortie de Unafold
    open( my $FOut, '<', $CT ) || die "Unable to open CT file: $!";

    my $cg = 0;
    my ( $nameSeq, $mfe, $longueur );
    while ( my $line = <$FOut> ) {
        my @variable = split( '\s+', $line );

        if (
            $variable[2] =~ m{
                          ^              # Begin of line
                          ENERGY         # ENERGY
                        }smx
          )
        {

            # Récupération du nom de la séquence
            $nameSeq = substr $variable[5], 0;

            # Récupération du MFE
            $mfe      = substr $variable[4], 0;
            $longueur = $variable[1];
            $cg       = 0;
        }
        else {    # Not the ENERGY line
            if (
                $variable[2] =~ m{
                              [CG]     # G or C
                            }smxi
              )
            {
                $cg++;
            }
            if (
                ( $longueur == $variable[1] )

                && ( $nameSeq eq $sequence_name )
              ) # tester si la longueur actuelle est égale à la longueur de la séquence
            {

                my $num = ( $mfe / $longueur ) * 100;
                my $other = $num / ( ( $cg / $longueur ) * 100 );
                open( my $RES, '>>', $MFEI_output ) || die "Unable to open: $!";
                my $content =
                  $nameSeq . "\t" . $other . "\t" . $mfe . "\t" . $num;
                print $RES $content;
                close $RES or die "Unable to close: $!";
                chmod 777, $MFEI_output;
            }

        }    # /Else− not ENERGY line
    }    # /While
    close $FOut or die "Unable to close: $!";
    return;
}

1;
