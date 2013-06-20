#!/usr/bin/perl -w

use warnings;
use strict;

use File::Path 'rmtree';
use File::Basename;
use Cwd qw( abs_path );
use File::Copy;
use PipelineMiRNA::Utils;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Components;
use PipelineMiRNA::PosterioriTests;

use Log::Message::Simple qw[msg error debug
  carp croak cluck confess];

my $local_dir = dirname( abs_path($0) );
my $rootdir = File::Spec->catdir( $local_dir, ".." );

## Data ##
my $dirData = File::Spec->catdir( $rootdir, 'data' );    # chemin séquence

## Code ##
my ( $icheck, $imfei, $irandfold, $iSC, $ialign, $idirJob, $iplant ) = @ARGV;
main_entry( $icheck, $imfei, $irandfold, $iSC, $ialign, $idirJob, $iplant );

sub main_entry {
    my ( $check, $mfei, $randfold, $SC, $align, $dirJob, $plant ) = @_;

    my $sequences_input = File::Spec->catfile( $dirJob, 'Sequences.fas' );
    if ( $check eq "checked" ) {

        #Filtering CDS
        PipelineMiRNA::Components::filter_CDS( $dirData, $dirJob, $plant );
    }
    else {
        my $sequence_uploaded =
          File::Spec->catfile( $dirJob, 'sequenceUpload.fas' );
        File::Copy::move( $sequence_uploaded, $sequences_input );
    }

    ##Passage du multifasta -> fasta et appel au script Stemloop
    open my $ENTREE_FH, '<', $sequences_input
      or die "Impossible d'ouvrir le fichier d'entree : $!";
    my %tab = PipelineMiRNA::Utils::parse_multi_fasta($ENTREE_FH);
    close $ENTREE_FH;

    foreach my $name ( keys %tab ) {

        my $temp_file = File::Spec->catfile( $dirJob, 'tempFile.txt' );
        open( my $TEMPFILE_FH, '>', $temp_file )
          or die "Impossible d'ouvrir le fichier d'entree  : $!";
        chmod 0777, $temp_file;
        print $TEMPFILE_FH "$name\n$tab{$name}";
        my $name = substr $name, 1;
        my $sequence_dir = File::Spec->catdir( $dirJob, $name );
        mkdir $sequence_dir;

        my $rnalfold_output =
          File::Spec->catfile( $sequence_dir, 'RNALfold.out' );
        PipelineMiRNA::Programs::run_rnalfold( $temp_file, $rnalfold_output )
          or die("Problem when running RNALfold");

        ## Appel de RNAstemloop
        my $rnastemloop_out =
          File::Spec->catfile( $sequence_dir, 'rnastemloop.out' );
        PipelineMiRNA::Programs::run_rnastemloop( $rnalfold_output,
            $rnastemloop_out )
          or die("Problem when running RNAstemloop");
        unlink $temp_file;
        process_RNAstemloop($rnastemloop_out);
    }
    process_tests( $dirJob, $mfei, $randfold, $SC, $align );
}

sub process_RNAstemloop {
    my ($rnastemloop_out) = @_;
    my $current_sequence_dir = dirname($rnastemloop_out);

    my $rnaeval_out =
      File::Spec->catfile( $current_sequence_dir, 'rnaeval.out' );

    PipelineMiRNA::Programs::run_rnaeval( $rnastemloop_out, $rnaeval_out )
      or die("Problem when running RNAeval");

    my ( $nameSeq, $dna, $Vienna );
    open( my $stem, '<', $rnastemloop_out ) or die $!;
    open( my $eval, '<', $rnaeval_out )     or die $!;
    open( my $out,  '>', 'output.txt' )     or die $!;

    #my $line = <$stem>;
    my $line2;

    while ( my $line = <$stem> ) {

        if ( ( $line =~ /^>(.*)/ ) ) {    # nom sequence
            $nameSeq = $1;
        }
        elsif ( ( $line =~ /^[a-zA-Z]/ ) ) { # récupération de la sequence adn
            $dna = substr $line, 0, -1;
            $line2 = substr( <$eval>, 0, -1 );    # the sequence as well
            if ( $dna ne $line2 ) {

                #This should not happen.
            }
        }
        elsif ( ( $line =~ /(.*)/ ) ) {
            $Vienna = $1;
            $line2  = <$eval>;    # the structure as well, and the energy
            if (
                my ( $structure, $energy ) =
                $line2 =~ m{
                        ^                #Begin of line
                        ([\.()]+?)       #A sequence of ( ) .
                        \s+?             #Some whitespace
                        \(               #Opening parenthesis
                            ([-.\d]*?)   #
                        \)               #Closing parenthesis
                        \s*$             #Whitespace until the end
                    }smx
              )
            {
                if ( $Vienna ne $structure ) {

                    #This should not happen.
                }
                {
                    my $candidate_dir =
                      File::Spec->catdir( $current_sequence_dir, $nameSeq );
                    mkdir $candidate_dir;

                    #Writing seq.txt
                    my $candidate_sequence =
                      File::Spec->catfile( $candidate_dir, 'seq.txt' );
                    open( my $OUT, '>', $candidate_sequence ) || die "$!";
                    print $OUT ">$nameSeq\n$dna\n";
                    close $OUT;

                    #Writing (pseudo) rnafold output
                    my $candidate_rnafold_output =
                      File::Spec->catfile( $candidate_dir, 'outRNAFold.txt' );
                    open( my $OUT2, '>', $candidate_rnafold_output )
                      || die "$!";
                    print $OUT2 ">$nameSeq\n$dna\n$structure ($energy)\n";
                    close $OUT2;
                }
            }
        }    #elsif
    }    #while $line=<IN>

    close($stem);
    close($eval);
    return 0;
}

sub process_tests {
    my ( $dirJob, $mfei, $randfold, $SC, $align ) = @_;
    ##Traitement fichier de sortie outStemloop
    opendir DIR, $dirJob;    #ouverture répertoire job
    my @dirs;
    @dirs = readdir DIR;
    closedir DIR;

    foreach my $dir (@dirs)    # parcours du contenu
    {
        if (   $dir ne "."
            && $dir ne ".."
            && -d $dirJob . $dir )    #si fichier est un répertoire
        {

            my $sequence_dir = File::Spec->catdir( $dirJob, $dir );

            opendir DIR, $sequence_dir;    # ouverture du sous répertoire
            my @files;
            @files = readdir DIR;
            closedir DIR;
            foreach my $file (@files) {
                if (   $file ne "."
                    && $file ne ".."
                    && -d File::Spec->catdir( $sequence_dir, $file )
                  )    # si le fichier est de type repertoire
                {

                    ####Traitement fichier de sortie outStemloop
                    my $candidate_dir =
                      File::Spec->catdir( $sequence_dir, $file );
                    chmod 0777, $candidate_dir;

                    my $seq_file =
                      File::Spec->catfile( $candidate_dir, 'seq.txt' );
                    my $candidate_rnafold_out =
                      File::Spec->catfile( $candidate_dir, 'outRNAFold.txt' );

                    ####conversion en format CT
                    my $candidate_ct_file =
                      File::Spec->catfile( $candidate_dir, 'outB2ct.ct' );
                    PipelineMiRNA::Programs::convert_to_ct(
                        $candidate_rnafold_out, $candidate_ct_file )
                      or die("Problem when converting to CT format");

                    my $varna_image =
                      File::Spec->catfile( $candidate_dir, 'image.png' );
                    PipelineMiRNA::Programs::run_varna( $candidate_ct_file,
                        $varna_image )
                      or die("Problem suring image generation using VARNA");

                    ## traitement du fichier OutVienna pour la récupération des données(Format Vienna, séquence ADN)
                    my $out_Vienna = File::Spec->catfile( $candidate_dir,
                        'outViennaTraited.txt' );
                    open( my $TRAITED_FH, '>', $out_Vienna ) || die "$!";
                    open( my $INPUT_FH, '<', $candidate_rnafold_out )
                      || die "$!";
                    my ( $nameSeq, $dna, $Vienna );
                    while ( my $line = <$INPUT_FH> ) {
                        if ( ( $line =~ /^>(.*)/ ) ) {    # nom sequence
                            $nameSeq = $1;
                        }
                        elsif ( ( $line =~ /^[a-zA-Z]/ ) )
                        {    # récupération de la sequence adn
                            $dna = substr $line, 0, -1;
                        }
                        elsif ( ( $line =~ /(.*) / ) ) {
                            $Vienna = $1;
                            print $TRAITED_FH $nameSeq . "\t" 
                              . $dna . "\t"
                              . $Vienna
                              . "\n";    #récupération du format Vienna
                        }
                    }
                    close $INPUT_FH;
                    close $TRAITED_FH;
                    system("chmod 777 $out_Vienna");
                    ####calcul MFEI (appel script energie.pl)
                    if ( $mfei eq "mfeiChecked" ) {
                        PipelineMiRNA::PosterioriTests::test_mfei(
                            $candidate_dir, $candidate_ct_file, $file );
                    }
                    ####calcul p-value randfold
                    if ( $randfold eq "randfoldChecked" ) {
                        PipelineMiRNA::PosterioriTests::test_randfold(
                            $candidate_dir, $seq_file );
                    }
                    ####calcul self-contain
                    if ( $SC eq "SCChecked" ) {
                        PipelineMiRNA::PosterioriTests::test_selfcontain(
                            $candidate_dir, $seq_file );
                    }
                    if ( $align eq "alignChecked" ) {
                        PipelineMiRNA::PosterioriTests::test_alignment(
                            $candidate_dir, $seq_file );
                    }
                }
            }
        }
    }
    return 0;
}
