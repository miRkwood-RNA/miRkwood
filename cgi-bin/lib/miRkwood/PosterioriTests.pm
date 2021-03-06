package miRkwood::PosterioriTests;

# ABSTRACT: Everything needed to run the a posteriori tests

use strict;
use warnings;

use miRkwood::MiRdup;
use miRkwood::Parsers;
use miRkwood::Programs;
use miRkwood::Utils;

use Log::Message::Simple qw[msg error debug];

use File::Which;
use File::Spec;
use Carp;
use Class::Struct;

sub new {
    my ( $class, @args ) = @_;
    my ($directory) = @args;
    my $self = {
        directory => $directory,
    };
    bless $self, $class;
    return $self;
}

sub get_directory {
    my ( $self, @args ) = @_;
    return $self->{'directory'};
}

=method test_randfold

Run the Randfold a posteriori test

=cut

sub test_randfold {
    my ( $self, @args ) = @_;
    my $seq_file = shift @args;
    my $res = 0;
    if ( ! which( 'RNAshuffles' ) ){
        debug('[WARNING] RNAshuffles is not installed. Cannot compute thermodynamic stability with RNAshuffles.', miRkwood->DEBUG());
    }
    else {
        my $randfold_out = File::Spec->catfile( $self->get_directory(), 'randfold.out' );
        miRkwood::Programs::run_randfold( $seq_file, $randfold_out, 200)
          or die('Problem when running Randfold');
        $res = miRkwood::Parsers::parse_pvalue($randfold_out);
    }
    return $res;
}


=method test_alignment

Run the Alignment (piccolo) a posteriori test

=cut

sub test_alignment {
    my ( $self, @args ) = @_;
    my $candidate_rnafold_stemploop_out = shift @args;
    my $cfg = miRkwood->CONFIG();
    my $mode = $cfg->param('job.pipeline');

    my $candidate_ct_stemloop_file =
      File::Spec->catfile( $self->get_directory(), 'outB2ct_stemloop.ct' );
    miRkwood::Programs::convert_to_ct( $candidate_rnafold_stemploop_out,
        $candidate_ct_stemloop_file )
      or warn("Problem when converting to CT format. File outB2ct_stemloop.ct non created.\n");

    my $seqN = File::Spec->catfile( $self->get_directory(), 'seqWithN.txt' );
    if ( $mode eq 'smallRNAseq' || ! -e $candidate_ct_stemloop_file ){
        my $seqWithT = $self->{'candidate'}{'sequence'};
        $seqWithT =~ s/U/T/g;
        open(my $FILE, '>', $seqN) or die "ERROR while creating $seqN: $!";
        print $FILE ">$self->{'candidate'}{'name'}\n";
        print $FILE "$seqWithT\n";
        close $FILE;
    }
    else{
        $self->mask_CT_file($candidate_ct_stemloop_file, $seqN);
    }

    my $exonerate_out = File::Spec->catfile( $self->get_directory(), 'alignement.txt' );
    my $alignments;
    my $created_aln_file = miRkwood::Programs::run_piccolo( $seqN, $exonerate_out );
    if ( $created_aln_file ){
        $alignments = $self->post_process_alignments($exonerate_out );
    }
    else {
        debug('[WARNING] No alignment file was created by piccolo.', miRkwood->DEBUG());
    }

    return $alignments;
}

=method post_process_alignments


=cut

sub post_process_alignments {
    my ( $self, @args ) = @_;
    my $file_alignement = shift @args;
    my %alignments;

    if (-z $file_alignement){
        return;
    }
    if (
        !eval {
            %alignments =
              miRkwood::Parsers::parse_custom_exonerate_output(
                $file_alignement);
        }
      )
    {
        # Catching exception
        carp("Exception when parsing piccolo output $file_alignement");
        return;
    }
    else {
        %alignments = $self->merge_alignments( \%alignments );
        return \%alignments;
    }
}

sub validate_alignments_with_mirdup {
    my ( $self, @args ) = @_;
    my $alignments = shift @args;
    my $tmp_file = File::Spec->catfile( $self->get_directory(), 'mirdup_validation.txt' );
    my %mirdup_results =
      miRkwood::MiRdup->validate_with_mirdup( $tmp_file, $self->{'sequence_name'},
        $self->{'candidate'}{'sequence'}, $self->{'candidate'}{'structure_stemloop'}, keys %{$alignments} );
    return \%mirdup_results;
}


=method merge_alignments

Merge overlapping alignments.
Given ordered positions, merge in [a..b] all [c..d] if a<=c and d<=b+2

=cut

sub merge_alignments {
    my ( $self, @args ) = @_;
    my $alignments = shift @args;
    my %alignments = %{$alignments};
    my %merged_alignments;
    my @keys = sort {
        ( miRkwood::Utils::get_element_of_split( $a, '-', 0 )
              <=> miRkwood::Utils::get_element_of_split( $b, '-', 0 ) )
          || ( miRkwood::Utils::get_element_of_split( $a, '-', 1 )
            <=> miRkwood::Utils::get_element_of_split( $b, '-', 1 ) )
    } keys %alignments;

    my @stocked_hits;
    my $final_left = -1;
    my $final_right = -1;

    foreach my $current_key (@keys) {
        my ( $current_left, $current_right ) = split( /-/, $current_key );

        # First alignment position
        if ( $final_left == -1 ){
            $final_left = $current_left;
        }
        if ( $final_right == -1 ){
            $final_right = $current_right;
        }

        if ( $current_right > $final_right + 3 ) {    # Position is too far, we don't merge
            # Store the alignments
            if (@stocked_hits) {
                push @{ $merged_alignments{ "$final_left-$final_right" } }, @stocked_hits;
            }
            # Reinitialise
            $final_left = $current_left;
            $final_right = $current_right;
            @stocked_hits = ();
        }
        else{   # Merge
            if ( $current_right > $final_right ){
                $final_right = $current_right;
            }
        }
        push @stocked_hits, @{ $alignments{$current_key} };
    }

    # Drop the remaining hits in the hash
    push @{ $merged_alignments{ "$final_left-$final_right" } }, @stocked_hits;
    return %merged_alignments;
}

# Sequence data structure

struct Sequence => {
    colonne1        => '@',
    colonne5        => '@',
    boucleterminale => '$',
    nbnucleotides   => '$',
    debBoucle       => '$',
    finBoucle       => '$',
    length          => '$',
    base            => '@',
};

=method mask_CT_file

Mask the CT file and outputting to boucleTermWithN_out file

=cut

sub mask_CT_file {
    my ( $self, @args ) = @_;
    my ( $CT, $boucleTermWithN_out ) = @args;

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
            push( @{$tab1}, $nbSeq );   # remplissage du tableau temporaire tab1
            push( @{$tab2}, $numSeq );  # remplissage du tableau temporaire   tab2
            push( @{$tab3}, $base );
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
                  )    # ici le premier nombre non nul après une suite de 0 et
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

sub validate_mirna_with_mirdup {
    my ( $self, @args ) = @_;

    my $input_file_for_mirdup = File::Spec->catfile( $self->get_directory(), 'mirdup_on_mirna.txt' );
    open (my $FILE, '>', $input_file_for_mirdup) or die "ERROR while creating $input_file_for_mirdup : $!";
    print $FILE 'mirna' . "\t" . $self->{'candidate'}{'mirna_sequence'} . "\t" . $self->{'candidate'}{'sequence'};
    close $FILE or die "ERROR while closing $FILE : $!";

    my $result_file = miRkwood::Programs::run_mirdup_validation_on_file($input_file_for_mirdup);
    if ( ! -e $result_file || ! -r $result_file ){
        return miRkwood::MiRdup->parse_validation_output($input_file_for_mirdup);
    }
    return miRkwood::MiRdup->parse_validation_output($result_file);
}

1;
