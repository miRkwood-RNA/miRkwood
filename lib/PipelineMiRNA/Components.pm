package PipelineMiRNA::Components;

# ABSTRACT: Components directly needed for the execution of the pipeline

use strict;
use warnings;

use YAML;
use File::Spec;
use Class::Struct;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Utils;
use PipelineMiRNA::Parsers;

struct Sequence => {   # déclaration de la structure de données (Sequence)
    colonne1        => '@',
    colonne5        => '@',
    boucleterminale => '$',
    nbnucleotides   => '$',
    debBoucle       => '$',
    finBoucle       => '$',
    length          => '$',
    base            => '@',
};

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
    chmod 0777, $input_sequences;
    chmod 0777, $blast_output;
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

sub mask_CT_file {
    my ( $CT, $boucleTermWithN_out ) = @_;

 # tableau associatif contenant le nom la séquence (clé) et un struct (valeur)
    my %tab = ();

    # tableau temporaire contenant les élements de la colonne 1 du fichier ct
    my $tab1 = [];

    # tableau temporaire contenant les élements de la colonne 5 du fichier ct
    my $tab2 = [];

    # tableau temporaire contenant les bases de la séquence
    my $tab3 = [];

    #ouverture du fichier de sortie de Unafold
    open( my $FOut, '<', $CT )
      or die "Error when opening $CT: $!";
    my ( $nbSeq, $numSeq, $base, $longueur, $nomSeq );
    while ( my $line = <$FOut> ) {

        my @variable = split( '\s+', $line );
        $nbSeq  = $variable[1];
        $numSeq = $variable[5];
        $base   = $variable[2];
        if ( $variable[2] =~ /^ENERGY/xms ) {
            $longueur = $variable[1];
            $nomSeq   = substr $variable[5],
              0;    # Récupération du nom de la séquence
            $tab1 = [];
            $tab2 = [];
            $tab3 = [];
        }
        else {
            push( @$tab1, $nbSeq );   # remplissage du tableau temporaire tab1
            push( @$tab2, $numSeq );  # remplissage du tableau temporaire   tab2
            push( @$tab3, $base );
            if ( $longueur == $variable[1]
              ) # tester si la longueur actuelle est égale à la longueur de la séquence
            {
                my $mySequence = Sequence->new();    # initialisation du struct
                $mySequence->colonne1($tab1);
                $mySequence->colonne5($tab2);
                $mySequence->boucleterminale(0);
                $mySequence->nbnucleotides(0);
                $mySequence->length($longueur);
                $mySequence->base($tab3);
                $tab{$nomSeq} = $mySequence
                  ; # remplissage du tableau associatif avec le nom et le struct
            }
        }
    }
    close $FOut or die "Problème à la fermeture : $!";

###parcours et traitement et remplissage du struct (nb nucléotides, position deb, position fin)####
    my ( $element5Copy, $nbnucleotides, $element5, $element1 );
    foreach my $sequence ( values %tab ) {
        $element5Copy  = -1;
        $nbnucleotides = 0;
        for ( my $i = 0 ;
            $i < $sequence->length ; $i++ )  # parcours des deux colonnes 1 et 5
        {
            $element1 = $sequence->colonne1($i);
            $element5 = $sequence->colonne5($i);
            if ( $element5 == 0 )            # lorsque l'on rencontre un 0
            {
                $nbnucleotides =
                  $nbnucleotides + 1;    #nombre des nucléotides de la boucle
            }
            else                         # tout ce qui n'est pas 0
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
    open( my $RES, '>>', $boucleTermWithN_out )
      or die "Error when opening $boucleTermWithN_out: $!";
    foreach my $sequenceNom ( keys %tab ) {
        print $RES '>' . $sequenceNom . "\n";
        for ( my $i = 0 ;
            $i < $tab{$sequenceNom}->length ;
            $i++ )    # parcours des deux colonnes 1 et 5
        {
            if (   ( $i >= $tab{$sequenceNom}->debBoucle - 1 )
                && ( $i <= $tab{$sequenceNom}->finBoucle - 1 ) )
            {
                print $RES 'N';
            }
            else {
                my $echo = $tab{$sequenceNom}->base($i);
                $echo =~ s/U/T/gi;
                print $RES $echo;
            }
        }
        print $RES "\n";
    }
    close $RES or die "Problème à la fermeture : $!";
    return $boucleTermWithN_out;
}

=method parse_custom_exonerate_output

Parse our custom Exonerate output
(as defined in Programs::run_exonerate)

Usage: parse_exonerate_alignment($alignment);

=cut

sub parse_custom_exonerate_output{
    my @args = @_;
    my $yaml_file = shift @args;
    (-e $yaml_file) or die("Error, YAML file $yaml_file does not exist");
    (! -z $yaml_file) or die("Error, YAML file $yaml_file is empty");

    my $yaml = YAML::LoadFile($yaml_file) or die("Error when parsing YAML file $yaml_file");
    my @contents = @{$yaml} or die("Error when parsing YAML file $yaml_file");

    my %results;
    foreach my $element (@contents){
         my $key = "$element->{'begin_target'}-$element->{'end_target'}";
         my $value = {
             'name' => $element->{'name'},
             'seq' => $element->{'seq'},
             'score' => $element->{'score'},
             'begin_target' => $element->{'begin_target'},
             'end_target' => $element->{'end_target'},
             'begin_query' => $element->{'begin_query'},
             'end_query' => $element->{'end_query'},
             'seq' => $element->{'seq'},
             'score' => $element->{'score'},
             'alignment' => parse_exonerate_alignment($element->{'alignment'}),
         };
         if (! exists $results{$key}) {
             $results{$key} = ();
         }
         push @{$results{$key}}, $value;
    }
    return %results;
}

=method parse_exonerate_alignment

Parse the alignment given in the output of Exonerate

Usage: parse_exonerate_alignment($alignment);

=cut

sub parse_exonerate_alignment {
    my @args = @_;
    my $alignment = shift @args;

    my $SPACE = q{ };
    my @top;
    my @middle;
    my @bottom;

    for (split /\n/mxs, $alignment ) {
        $_ =~ s/T/U/g;
        $_ =~ s/t/u/g;
        my ($first, $second, $label) = split($SPACE, $_);
        if ($label ne 'none') {
            push( @top,    $second );
            push( @bottom, $first );
            if ( uc $first eq uc $second ) {
                push( @middle, '|' );
            }else{
                push( @middle, $SPACE );
            }
        }
    }
    my $top = join('', @top);
    my $middle = join('', @middle);
    my $bottom= join('', @bottom);

    my $result = <<"END";
$top
$middle
$bottom
END
    return $result;
}

1;
