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

my ( $idirBlast, $idirData, $idirJob, $iplant ) = @ARGV;
main( $idirBlast, $idirData, $idirJob, $iplant );

sub main {
    my ( $dirBlast, $dirData, $dirJob, $plant ) = @_;

    my $uploaded_sequences = File::Spec->catfile( $dirJob, 'sequenceUpload.fas' );
    my $input_sequences    = File::Spec->catfile( $dirJob, 'Sequences.fas' );

    my $blast_database = File::Spec->catfile( $dirData, "$plant.fas" );
    my $blast_output = File::Spec->catfile( $dirJob, 'outBlast.txt' );
    my $blastx_options = "-outfmt 6 -max_target_seqs 1 -evalue 1E-5";
    PipelineMiRNA::Programs::run_blast($uploaded_sequences,
                                       $blast_database,
                                       $blastx_options,
                                       $blast_output)
      or die("Problem when running Blastx");
    chmod 777, $blast_output;

    ##Filtrage du ficher : Elimination des régions codantes###
    my @SeqNamesOUT =
      ();    # tableau contenant la liste des séquences issues du Blast
    my @SeqNames = ();    # tableau contenant la liste des séquences à traiter
    my @SeqDiff  = ();    # tableau contenant la liste des séquences non codantes
    my $i        = 0;

# Ouverture du fichier outBlast.TXT et construction du tableau contenant la liste des séquences issues du blast
    open( my $FOut, '<', $blast_output )
      || die "Problème à l\'ouverture : $!";
    while ( my $line = <$FOut> ) {
        my @name = split( '\t', $line );
        push( @SeqNamesOUT, $name[0] );
    }
    close $FOut || die "Problème à la fermeture : $!";

# Ouverture du fichier sequenceUpload.fas et construction du tableau contenant la liste des séquences issues du blast
    open( my $FSeq, '<', $uploaded_sequences )
      || die "Problème à l\'ouverture : $!";
    while ( my $line = <$FSeq> ) {
        if ( grep { /^>/msx } $line ) {
            my $lineSeq = substr $line, 1, -1;
            push( @SeqNames, $lineSeq );
        }

    }
    close $FSeq || die "Problème à la fermeture : $!";

    #Remplissage du tableau @SeqDiff avec les séquences non codantes
    foreach my $SeqName (@SeqNames) {
        my $exists = 1;
        foreach my $SeqNameOUT (@SeqNamesOUT) {
            if ( $SeqName eq $SeqNameOUT ) {
                $exists = 0;
            }
        }
        if ( $exists ) {
            push( @SeqDiff, $SeqName );

            #print "\nUne seq ".$SeqName."\n";
        }
        $i = $i++;
    }

    # Création du fichier FASTA filtrée des régions codantes
    my $record = 0;
    $i      = 0;
    open( my $RES,  '>>', $input_sequences )
      or die "Problème à l\'ouverture : $!";
    open( my $FSeq2, '<', $uploaded_sequences )
      or die "Problème à l\'ouverture : $!";
    while ( my $line = <$FSeq2> ) {
        if ( grep { /^>/msx } $line ) {
            my $lineSeq = substr $line, 1, -1;
            if ( $SeqDiff[$i] eq $lineSeq ) {
                $record = 1;
                $i++;
            }
            else {
                $record = 0;
            }
        }
        if ( $record ) {
            printf $RES $line;
        }
    }
    close $FSeq2 || die "Problème à la fermeture : $!";
    close $RES   || die "Problème à la fermeture : $!";
    chmod 777, $input_sequences;
    return;
}
