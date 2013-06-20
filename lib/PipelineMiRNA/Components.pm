#!/usr/local/bin/perl
######################################################################
#  Lancement du script : perl filtrage.pl sequenceTest DBArapidopsis #
#  @author : Mohcen BENMOUNAH                                        #
#                                                                    #
#                                                                    #
######################################################################
use strict;
use warnings;

use File::Spec;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Utils;

my ( $idirData, $idirJob, $iplant ) = @ARGV;
main( $idirData, $idirJob, $iplant );

sub filter_CDS {
    my ( $dirData, $dirJob, $plant ) = @_;

    my $uploaded_sequences = File::Spec->catfile( $dirJob, 'sequenceUpload.fas' );
    my $input_sequences    = File::Spec->catfile( $dirJob, 'Sequences.fas' );

    my $blast_database = File::Spec->catfile( $dirData, "$plant.fas" );
    my $blast_output = File::Spec->catfile( $dirJob, 'outBlast.txt' );
    my $blastx_options = '-outfmt 6 -max_target_seqs 1 -evalue 1E-5';
    PipelineMiRNA::Programs::run_blast($uploaded_sequences,
                                       $blast_database,
                                       $blastx_options,
                                       $blast_output
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
