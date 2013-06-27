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
    my $debug = 0;
    debug('BEGIN execute_scripts', $debug);
    my $sequences_input = File::Spec->catfile( $dirJob, 'Sequences.fas' );
    if ( $check eq "checked" ) {
        debug("FilteringCDS", $debug);
        #Filtering CDS
        PipelineMiRNA::Components::filter_CDS( $dirData, $dirJob, $plant );
    }
    else {
        my $sequence_uploaded =
          File::Spec->catfile( $dirJob, 'sequenceUpload.fas' );
        debug("Moving file $sequence_uploaded to $sequences_input", $debug);
        File::Copy::move( $sequence_uploaded, $sequences_input );
    }

    ##Passage du multifasta -> fasta et appel au script Stemloop
    debug("Opening multifasta $sequences_input", $debug);
    open my $ENTREE_FH, '<', $sequences_input
      or die "Error when opening sequences -$sequences_input-: $!";
    my %tab = PipelineMiRNA::Utils::parse_multi_fasta($ENTREE_FH);
    close $ENTREE_FH;
    
    debug("Iterating over names", $debug);
    debug("Hop: %tab", $debug);
    foreach my $name ( keys %tab ) {
        debug("Considering $name", $debug);
        my $temp_file = File::Spec->catfile( $dirJob, 'tempFile.txt' );
        open( my $TEMPFILE_FH, '>', $temp_file )
          or die "Error when opening tempfile -$temp_file-: $!";
        chmod 0777, $temp_file;
        print $TEMPFILE_FH "$name\n$tab{$name}";
        my $name = substr $name, 1;
        my $sequence_dir = File::Spec->catdir( $dirJob, $name );
        mkdir $sequence_dir;

        my $rnalfold_output =
          File::Spec->catfile( $sequence_dir, 'RNALfold.out' );
        debug("Running RNAfold", $debug);
        PipelineMiRNA::Programs::run_rnalfold( $temp_file, $rnalfold_output )
          or die("Problem when running RNALfold: $!");

        ## Appel de RNAstemloop
        my $rnastemloop_out_optimal =
          File::Spec->catfile( $sequence_dir, 'rnastemloop_optimal.out' );
        my $rnastemloop_out_stemloop =
          File::Spec->catfile( $sequence_dir, 'rnastemloop_stemloop.out' );
        PipelineMiRNA::Programs::run_rnastemloop( $rnalfold_output,
            $rnastemloop_out_optimal, $rnastemloop_out_stemloop)
          or die("Problem when running RNAstemloop");
        unlink $temp_file;
        process_RNAstemloop($rnastemloop_out_optimal, 'optimal');
        process_RNAstemloop($rnastemloop_out_stemloop, 'stemloop');
    }
    process_tests( $dirJob, $mfei, $randfold, $SC, $align );
}

sub process_RNAstemloop {
    my ($rnastemloop_out, $suffix) = @_;
    my $current_sequence_dir = dirname($rnastemloop_out);

    my $rnaeval_out =
      File::Spec->catfile( $current_sequence_dir, "rnaeval_$suffix.out" );

    PipelineMiRNA::Programs::run_rnaeval( $rnastemloop_out, $rnaeval_out )
      or die("Problem when running RNAeval");

    my ( $nameSeq, $dna, $Vienna );
    open( my $stem, '<', $rnastemloop_out ) or die $!;
    open( my $eval, '<', $rnaeval_out )     or die $!;

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
                      File::Spec->catfile( $candidate_dir, "seq.txt" );
                    open( my $OUT, '>', $candidate_sequence ) || die "$!";
                    print $OUT ">$nameSeq\n$dna\n";
                    close $OUT;

                    #Writing (pseudo) rnafold output
                    my $candidate_rnafold_output =
                      File::Spec->catfile( $candidate_dir, "outRNAFold_$suffix.txt" );
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
                    my $candidate_rnafold_optimal_out =
                      File::Spec->catfile( $candidate_dir, 'outRNAFold_optimal.txt' );
                    my $candidate_rnafold_stemploop_out =
                      File::Spec->catfile( $candidate_dir, 'outRNAFold_stemloop.txt' );

                    ####conversion en format CT
                    my $candidate_ct_optimal_file =
                      File::Spec->catfile( $candidate_dir, 'outB2ct_optimal.ct' );
                    PipelineMiRNA::Programs::convert_to_ct(
                        $candidate_rnafold_optimal_out, $candidate_ct_optimal_file )
                      or die("Problem when converting to CT format");

                    my $candidate_ct_stemloop_file =
                      File::Spec->catfile( $candidate_dir, 'outB2ct_stemloop.ct' );
                    PipelineMiRNA::Programs::convert_to_ct(
                        $candidate_rnafold_stemploop_out, $candidate_ct_stemloop_file )
                      or die("Problem when converting to CT format");

                    my $varna_image =
                      File::Spec->catfile( $candidate_dir, 'image.png' );
                    PipelineMiRNA::Programs::run_varna( $candidate_ct_stemloop_file,
                        $varna_image )
                      or die("Problem during image generation using VARNA");

                    ## traitement du fichier OutVienna pour la récupération des données(Format Vienna, séquence ADN)
                    my $out_Vienna = File::Spec->catfile( $candidate_dir,
                        'outViennaTraited.txt' );
                    open( my $TRAITED_FH, '>', $out_Vienna ) || die "$!";
                    open( my $INPUT_FH, '<', $candidate_rnafold_optimal_out ) #TODO: Check if correct
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
                            $candidate_dir, $candidate_ct_optimal_file, $file );
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
                            $candidate_dir, $candidate_ct_stemloop_file );
                    }
                }
            }
        }
    }
    return 0;
}
