package PipelineMiRNA::MainPipeline;

use strict;
use warnings;

use CGI::Carp qw(fatalsToBrowser);
use File::Path 'rmtree';
use File::Basename;
use Cwd qw( abs_path );
use File::Copy;
use PipelineMiRNA::Paths;
use PipelineMiRNA::Utils;
use PipelineMiRNA::Parsers;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Components;
use PipelineMiRNA::PosterioriTests;
use Data::Dumper;

use Log::Message::Simple qw[msg error debug];


### Data ##
my $dirData = PipelineMiRNA::Paths->get_absolute_path( 'data' );
#

sub main_entry {
    my ( $check, $mfei, $randfold, $SC, $align, $dirJob, $plant ) = @_;
    my $debug = 1;
    my $log_file = File::Spec->catfile( $dirJob, 'log.log' );
    open( my $LOG, '>>', $log_file ) || die "Error when opening log file: $!";
    local $Log::Message::Simple::DEBUG_FH   = $LOG;

    debug('BEGIN execute_scripts', $debug);
    my $sequences_input = File::Spec->catfile( $dirJob, 'Sequences.fas' );
    if ( $check eq 'checked' ) {
        debug('FilteringCDS', $debug);
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
    debug("Calling parse_multi_fasta() on $sequences_input", $debug);
    my %tab = PipelineMiRNA::Utils::parse_multi_fasta($ENTREE_FH);
    close $ENTREE_FH;

    debug('Iterating over names', $debug);

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
        debug('Running RNAfold', $debug);
        PipelineMiRNA::Programs::run_rnalfold( $temp_file, $rnalfold_output )
          or die("Problem when running RNALfold: $!");

        ## Appel de RNAstemloop
        my $rnastemloop_out_optimal =
          File::Spec->catfile( $sequence_dir, 'rnastemloop_optimal.out' );
        my $rnastemloop_out_stemloop =
          File::Spec->catfile( $sequence_dir, 'rnastemloop_stemloop.out' );
        debug("Running RNAstemloop on $rnalfold_output", $debug);
        PipelineMiRNA::Programs::run_rnastemloop( $rnalfold_output,
            $rnastemloop_out_optimal, $rnastemloop_out_stemloop)
          or die("Problem when running RNAstemloop");
        unlink $temp_file;
        process_RNAstemloop_wrapper($rnastemloop_out_optimal, 'optimal');
        process_RNAstemloop_wrapper($rnastemloop_out_stemloop, 'stemloop');
    }
    process_tests( $dirJob, $mfei, $randfold, $SC, $align );
    return;
}

sub process_RNAstemloop_wrapper {
    my ($rnastemloop_out, $suffix) = @_;
    my $current_sequence_dir = dirname($rnastemloop_out);
    debug("Processing RNAstemloop output for $suffix $rnastemloop_out", 1);
    my $rnaeval_out =
      File::Spec->catfile( $current_sequence_dir, "rnaeval_$suffix.out" );

    debug("Running RNAeval in $rnaeval_out", 1);
    PipelineMiRNA::Programs::run_rnaeval( $rnastemloop_out, $rnaeval_out )
      or die("Problem when running RNAeval");

    open( my $stem, '<', $rnastemloop_out ) or die $!;
    open( my $eval, '<', $rnaeval_out )     or die $!;
    process_RNAstemloop($current_sequence_dir, $suffix, $stem, $eval);
    close($stem);
    close($eval);
    return 0;
}

sub process_RNAstemloop {
    my @args   = @_;
    my ($current_sequence_dir) = shift @args;
    my ($suffix) = shift @args;
    my ($stem)   = shift @args;
    my ($eval)   = shift @args;

    my $line2;
    my ( $nameSeq, $dna, $Vienna );
    while ( my $line = <$stem> ) {

        if ( ( $line =~ /^>(.*)/ ) ) {    # nom sequence
            $nameSeq = $1;
        }
        elsif ( ( $line =~ /^[a-zA-Z]/ ) ) { # récupération de la sequence adn
            $dna = substr $line, 0, -1;
            $line2 = substr( <$eval>, 0, -1 );    # the sequence as well

            if ( $dna ne $line2 ) {
                # Should not happen
            }

        }
        elsif ( ( $line =~ /(.*)/ ) ) {
            $Vienna = $1;
            $line2  = <$eval>;    # the structure as well, and the energy
            if (
                my ( $structure, $energy ) = PipelineMiRNA::Parsers::parse_Vienna_line($line2)
               )
            { # We have a structure
                if ( $Vienna ne $structure ) {
                }
                my $candidate_dir =
                  File::Spec->catdir( $current_sequence_dir, $nameSeq );
                mkdir $candidate_dir;

                #Writing seq.txt
                my $candidate_sequence =
                  File::Spec->catfile( $candidate_dir, 'seq.txt' );
                open( my $OUT, '>', $candidate_sequence )
                  or die "Error when opening $candidate_sequence: $!";
                print $OUT ">$nameSeq\n$dna\n";
                close $OUT;

                #Writing (pseudo) rnafold output
                my $candidate_rnafold_output =
                  File::Spec->catfile( $candidate_dir, "outRNAFold_$suffix.txt" );

                open( my $OUT2, '>', $candidate_rnafold_output )
                  or die "Error when opening $candidate_rnafold_output: $!";
                print $OUT2 ">$nameSeq\n$dna\n$structure ($energy)\n";
                close $OUT2;

            } else {
                debug("No structure found in $line2", 1);
            } # if $line2
        }else{
            # Should not happen
        } #if $line1
    } #while $line=<IN>
}

=method process_tests

Perform the a posteriori tests for a given job

=cut

sub process_tests {
    my ( $dirJob, $mfei, $randfold, $SC, $align ) = @_;
    debug("A posteriori tests in $dirJob", 1);
    ##Traitement fichier de sortie outStemloop
    opendir DIR, $dirJob;    #ouverture répertoire job
    my @dirs;
    @dirs = readdir DIR;
    closedir DIR;

    foreach my $dir (@dirs)    # parcours du contenu
    {
        debug("Considering $dir", 1);
        my $sequence_dir = File::Spec->catdir( $dirJob, $dir );
        if (   $dir ne '.'
            && $dir ne '..'
            && -d $sequence_dir )    #si fichier est un répertoire
        {
            debug("Entering sequence $sequence_dir", 1);
            opendir DIR, $sequence_dir;    # ouverture du sous répertoire
            my @files;
            @files = readdir DIR;
            closedir DIR;
            foreach my $file (@files) {
                debug("Considering $file", 1);
                my $candidate_dir =
                      File::Spec->catdir( $sequence_dir, $file );
                if (   $file ne '.'
                    && $file ne '..'
                    && -d $candidate_dir
                  )    # si le fichier est de type repertoire
                {
                    debug("Entering candidate $file", 1);
                    process_tests_for_candidate($candidate_dir, $file, $mfei, $randfold, $SC, $align);
                    debug("Done with candidate $file", 1);
                } # foreach my $file (@files)
            } # if directory
            debug("Done with initial sequence $dir", 1);
        } # foreach my $dir (@dirs)
    } #process_tests
    return 0;
}

=method process_tests_for_candidate

Perform the a posteriori tests for a given candidate

=cut

sub process_tests_for_candidate {

    my @args = @_;
    my ($candidate_dir, $file, $mfei, $randfold, $SC, $align ) = @args;

    ####Traitement fichier de sortie outStemloop
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
    debug("Converting optimal to CT in $candidate_ct_optimal_file", 1);
    PipelineMiRNA::Programs::convert_to_ct(
        $candidate_rnafold_optimal_out, $candidate_ct_optimal_file )
      or die('Problem when converting to CT format');

    my $candidate_ct_stemloop_file =
      File::Spec->catfile( $candidate_dir, 'outB2ct_stemloop.ct' );
    debug("Converting stemloop to CT in $candidate_ct_stemloop_file", 1);
    PipelineMiRNA::Programs::convert_to_ct(
        $candidate_rnafold_stemploop_out, $candidate_ct_stemloop_file )
      or die('Problem when converting to CT format');

    my $varna_image =
      File::Spec->catfile( $candidate_dir, 'image.png' );
    debug("Generating image using VARNA in $varna_image", 1);
    PipelineMiRNA::Programs::run_varna( $candidate_ct_stemloop_file,
        $varna_image )
      or die('Problem during image generation using VARNA');

    ## traitement du fichier OutVienna pour la récupération des données(Format Vienna, séquence ADN)
    my $out_Vienna = File::Spec->catfile( $candidate_dir,
        'outViennaTraited.txt' );
    open( my $TRAITED_FH, '>', $out_Vienna )
      or die "Error when opening $out_Vienna: $!";
    open( my $INPUT_FH, '<', $candidate_rnafold_optimal_out ) #TODO: Check if correct
      or die "Error when opening $candidate_rnafold_optimal_out: $!";
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
    chmod 777, $out_Vienna;
    ####calcul MFEI (appel script energie.pl)
    if ( $mfei eq 'mfeiChecked' ) {
        debug("Running test_mfei on $file", 1);
        PipelineMiRNA::PosterioriTests::test_mfei(
            $candidate_dir, $candidate_ct_optimal_file, $file );
    }
    ####calcul p-value randfold
    if ( $randfold eq 'randfoldChecked' ) {
        debug("Running test_randfold on $seq_file", 1);
        PipelineMiRNA::PosterioriTests::test_randfold(
            $candidate_dir, $seq_file );
    }
    ####calcul self-contain
    if ( $SC eq 'SCChecked' ) {
        debug("Running test_selfcontain on $seq_file", 1);
        PipelineMiRNA::PosterioriTests::test_selfcontain(
            $candidate_dir, $seq_file );
    }
    if ( $align eq 'alignChecked' ) {
        debug("Running test_alignment on $candidate_ct_stemloop_file", 1);
        PipelineMiRNA::PosterioriTests::test_alignment(
            $candidate_dir, $candidate_ct_stemloop_file );
    } # if file
    return;
}

1;
