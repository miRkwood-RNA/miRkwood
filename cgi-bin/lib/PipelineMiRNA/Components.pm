package PipelineMiRNA::Components;

# ABSTRACT: Components directly needed for the execution of the pipeline

use strict;
use warnings;

use YAML::XS;
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

=method filter_CDS

Filter the sequences based on Blast results

=cut

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

    my %blast_seqs = %{PipelineMiRNA::Parsers::index_blast_output($blast_output)};

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

=method compute_energy

Parse the CT file and fill the MFEI output file.

=cut

sub compute_energy {
    my ( $CT_file, $MFEI_output_file, $sequence_name ) = @_;

    # ouverture du fichier de sortie de Unafold
    open( my $FOut, '<', $CT_file ) || die "Unable to open CT file: $!";

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
            if ( $longueur == $variable[1] )
            # tester si la longueur actuelle est égale à la longueur de la séquence
            {

                my $num = ( $mfe / $longueur ) * 100;
                my $other = PipelineMiRNA::Utils::restrict_num_decimal_digits($num / ( ( $cg / $longueur ) * 100 ),3);
                $num = PipelineMiRNA::Utils::restrict_num_decimal_digits($num, 3);
                open( my $RES, '>', $MFEI_output_file ) || die "Unable to open $MFEI_output_file: $!";
                my $content =
                  $nameSeq . "\t" . $other . "\t" . $mfe . "\t" . $num;
                print $RES $content;
                close $RES or die "Unable to close: $!";
            }

        }    # /Else− not ENERGY line
    }    # /While
    close $FOut or die "Unable to close: $!";
    return $MFEI_output_file;
}

=method mask_CT_file

Mask the CT file and outputting to boucleTermWithN_out file

=cut

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
    open( my $RES, '>', $boucleTermWithN_out )
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
    check_yaml_file($yaml_file);

    my $yaml = YAML::XS::LoadFile($yaml_file) or die("Error when parsing YAML file $yaml_file");
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
             'strand_target' => $element->{'strand_target'},
             'begin_query' => $element->{'begin_query'},
             'end_query' => $element->{'end_query'},
             'def_query' => $element->{'def_query'},
             'seq' => $element->{'seq_query'},
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
            push( @top,    lc $second );
            push( @bottom, lc $first );
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

=method get_data_from_rnafold_out

Retrieve sequence name, position, sequence & structure from
a RNAfold output.

=cut

sub get_data_from_rnafold_out {
    my @args                            = @_;
    my $candidate_rnafold_stemploop_out = shift @args;
    my @vienna_res = PipelineMiRNA::Parsers::parse_RNAfold_output(
        $candidate_rnafold_stemploop_out);
    my @position    = split( /__/, $vienna_res[0] );
    my $name        = $position[0];
    my $position    = $position[1];
    my $DNASequence = $vienna_res[1];
    my $Vienna      = $vienna_res[2];
    return ( $name, $position, $DNASequence, $Vienna );
}

=method check_yaml_file

Check whether a given file is suitable for YAML deserialization.

=cut

sub check_yaml_file {
    my @args      = @_;
    my $yaml_file = shift @args;
    ( -e $yaml_file )  or die("Error, YAML file $yaml_file does not exist");
    ( !-z $yaml_file ) or die("Error, YAML file $yaml_file is empty");
    ( -r $yaml_file )  or die("Error, YAML file $yaml_file is not readable");
    return;
}

=method merge_alignments

Merge overlapping alignments.
Given ordered positions, merge in [a..b] all [c..d] if a<=c and d<=b+2

=cut

sub merge_alignments {
    my (@args)     = @_;
    my $alignments = shift @args;
    my %alignments = %{$alignments};

    my %merged_alignments;
    my ( $stocked_left, $stocked_right ) = ( -10, -10 );

    my @keys = sort {
        ( PipelineMiRNA::Utils::get_element_of_split( $a, '-', 0 )
              <=> PipelineMiRNA::Utils::get_element_of_split( $b, '-', 0 ) )
          || ( PipelineMiRNA::Utils::get_element_of_split( $a, '-', 1 )
            <=> PipelineMiRNA::Utils::get_element_of_split( $b, '-', 1 ) )
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
