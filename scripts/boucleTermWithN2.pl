#!/usr/local/bin/perl
######################################################################
# Script permettant de détecter les boucles terminales               #
#                                                                    #
#	Lancement du script : perl boucleTerm.pl Fichier.txt.ct  	     #
#  @author : Mohcen BENMOUNAH                                        #                                                                	 #										         #
#                                                                    #
######################################################################
use Class::Struct;
my ( $CT, $dirData ) = @ARGV;

# tableau associatif contenant le nom la séquence (clé) et un struct (valeur)
%tab = ();

# tableau temporaire contenant les élements de la colonne 1 du fichier ct
$tab1 = [];

# tableau temporaire contenant les élements de la colonne 5 du fichier ct
$tab2 = [];

# tableau temporaire contenant les bases de la séquence
$tab3 = [];

struct Sequence => {    # déclaration de la structure de données (Sequence)
    colonne1        => '@',
    colonne5        => '@',
    boucleterminale => '$',
    nbnucleotides   => '$',
    debBoucle       => '$',
    finBoucle       => '$',
    length          => '$',
    base            => '@',
};

#ouverture du fichier de sortie de Unafold
open( FOut, $CT )
  or die "Error when opening $CT: $!";
while ( my $line = <FOut> ) {

    my @variable = split( '\s+', $line );
    $nbSeq  = $variable[1];
    $numSeq = $variable[5];
    $base   = $variable[2];
    if ( $variable[2] =~ /^ENERGY/ ) {
        $longueur = $variable[1];
        $nomSeq   = substr $variable[5],
          0;                          # Récupération du nom de la séquence
        $tab1 = [];
        $tab2 = [];
        $tab3 = [];
    }
    else {
        push( @$tab1, $nbSeq );       # remplissage du tableau temporaire tab1
        push( @$tab2, $numSeq );      # remplissage du tableau temporaire	tab2
        push( @$tab3, $base );
        if ( $longueur == $variable[1]
          ) # tester si la longueur actuelle est égale à la longueur de la séquence
        {
            $mySequence = Sequence->new();    # initialisation du struct
            $mySequence->colonne1($tab1);
            $mySequence->colonne5($tab2);
            $mySequence->boucleterminale(0);
            $mySequence->nbnucleotides(0);
            $mySequence->length($longueur);
            $mySequence->base($tab3);
            $tab{$nomSeq} = $mySequence
              ;    # remplissage du tableau associatif avec le nom et le struct
        }
    }
}
close FOut or die "Problème à la fermeture : $!";

###parcours et traitement et remplissage du struct (nb nucléotides, position deb, position fin)####
foreach $sequence ( values %tab ) {
    $element5Copy  = -1;
    $nbnucleotides = 0;
    for ( $i = 0 ;
        $i < $sequence->length ; $i++ )    # parcours des deux colonnes 1 et 5
    {
        $element1 = $sequence->colonne1($i);
        $element5 = $sequence->colonne5($i);
        if ( $element5 == 0 )              # lorsque l'on rencontre un 0
        {
            $nbnucleotides =
              $nbnucleotides + 1;          #nombre des nucléotides de la boucle
        }
        else                               # tout ce qui n'est pas 0
        {
            if (   ( $i != 0 )
                && ( $sequence->colonne5( $i - 1 ) == 0 )
                && ( $element1 == $element5Copy )
              )    # ici le premier nombre non nul arpès une suite de 0 et
            {
                $sequence->boucleterminale(1);
                $sequence->nbnucleotides($nbnucleotides);
                $sequence->debBoucle( $element5Copy - $nbnucleotides );
                $sequence->finBoucle( $element5Copy - 1 );
                last;
            }
            $nbnucleotides = 0;
            $element5Copy  = $element5;
        }
    }
}
## lecture du tableau associatif et generation fichier boucleTermWithN.txt ##
open( RES, '>>' . $dirData . '/boucleTermWithN.txt' );
foreach $sequenceNom ( keys %tab ) {
    print RES ">" . $sequenceNom . "\n";
    for ( $i = 0 ;
        $i < $tab{$sequenceNom}->length ;
        $i++ )    # parcours des deux colonnes 1 et 5
    {
        if (   ( $i >= $tab{$sequenceNom}->debBoucle - 1 )
            && ( $i <= $tab{$sequenceNom}->finBoucle - 1 ) )
        {
            print RES 'N';
        }
        else {
            print RES $tab{$sequenceNom}->base($i);

        }
    }
    print RES "\n";
}
close RES;
