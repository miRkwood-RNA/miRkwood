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

my ( $idirData, $idirJob, $iplant ) = @ARGV;
main( $idirData, $idirJob, $iplant );

sub main {
    my ( $dirData, $dirJob, $plant ) = @_;

    my $uploaded_sequences = File::Spec->catfile( $dirJob, 'sequenceUpload.fas' );
    my $input_sequences    = File::Spec->catfile( $dirJob, 'Sequences.fas' );

    my $blast_database = File::Spec->catfile( $dirData, "$plant.fas" );
    my $blast_output = File::Spec->catfile( $dirJob, 'outBlast.txt' );
    my $blastx_options = '-outfmt 6 -max_target_seqs 1 -evalue 1E-5';
    PipelineMiRNA::Programs::run_blast($uploaded_sequences,
                                       $blast_database,
                                       $blastx_options,
                                       $blast_output)
      or die('Problem when running Blastx');
    chmod 777, $blast_output;

    ##Filtrage du ficher : Elimination des régions codantes###
    my @SeqNamesOUT =
      ();    # tableau contenant la liste des séquences issues du Blast
    my @SeqNames = ();    # tableau contenant la liste des séquences à traiter
    my @SeqDiff  = ();    # tableau contenant la liste des séquences non codantes
    my $i        = 0;

    my %blast_seqs;

    # Looping in blast_output, indexing sequences found
    open( my $FOut, '<', $blast_output )
      || die "Problème à l\'ouverture : $!";
    while ( my $line = <$FOut> ) {
        my @name = split( '\t', $line );
        $blast_seqs{$name[0]} = 1;
    }
    close $FOut || die "Problème à la fermeture : $!";

    # Looping in uploaded_sequences, copying only the files not found in Blast
    open( my $RES,  '>>', $input_sequences )
      or die "Problème à l\'ouverture : $!";
    open( my $FSeq, '<', $uploaded_sequences )
      || die "Problème à l\'ouverture : $!";
    my $lineSeq;
    while ( my $line = <$FSeq> ) {
        if ( grep { /^>/msx } $line ) {
            $lineSeq = substr $line, 1, -1;
        }
        if(!exists($blast_seqs{$lineSeq})) {
            printf $RES $line;
        }
    }
    close $FSeq || die "Problème à la fermeture : $!";
    close $RES  || die "Problème à la fermeture : $!";
    chmod 777, $input_sequences;
    return $input_sequences;
}
