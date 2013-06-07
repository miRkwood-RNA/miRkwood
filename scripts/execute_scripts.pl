#!/usr/bin/perl -w

use warnings;
use strict;

use File::Path 'rmtree';
use File::Basename;
use Cwd qw( abs_path );
use PipelineMiRNA::Utils;

my $local_dir = dirname( abs_path($0) );
my $rootdir = File::Spec->catdir( $local_dir, ".." );

my $dirScript = File::Spec->catdir( $rootdir, 'scripts' );   # chemin script
my $dirProgs  = File::Spec->catdir( $rootdir, 'programs' );  # chemin programmes
my $dirImages = File::Spec->catdir( $rootdir, 'images' );    # chemin images
my $dirData   = File::Spec->catdir( $rootdir, 'data' );      # chemin séquence

## Programs ##

my $vienna_dir   = File::Spec->catfile( $dirProgs,   'ViennaRNA-2.1.2' );
my $rnafold_bin  = File::Spec->catfile( $vienna_dir, 'Progs', 'RNAfold' );
my $rnalfold_bin = File::Spec->catfile( $vienna_dir, 'Progs', 'RNALfold' );
my $rnaeval_bin  = File::Spec->catfile( $vienna_dir, 'Progs', 'RNAeval' );
my $b2ct_bin     = File::Spec->catfile( $vienna_dir, 'Utils', 'b2ct' );

my $randfold_bin = File::Spec->catfile( $dirProgs, 'randfold-2.0', 'randfold' );
my $selfcontain_bin =
  File::Spec->catfile( $dirProgs, 'selfcontain_unix', 'selfcontain.py' );
my $exonerate_bin =
  File::Spec->catfile( $dirProgs, 'exonerate-2.2.0-i386', 'bin', 'exonerate' );
my $varna_bin        = File::Spec->catfile( $dirProgs, 'VARNAv3-9.jar' );
my $rnastemploop_bin = File::Spec->catfile( $dirProgs, 'RNAstemloop' );

my $dirBlast = File::Spec->catdir( $dirProgs, 'ncbi-blast-2.2.28+-src', 'c++',
    'GCC460-Debug', 'bin' );    # chemin Blast

## Data ##

my $mirbase_file = File::Spec->catfile( $dirData, 'MirbaseFile.txt' );
my $matrix_file  = File::Spec->catfile( $dirData, 'matrix' );

## Code ##
my $DEBUG = 0;
my ( $check, $mfei, $randfold, $SC, $align, $dirJob, $plant ) = @ARGV;

if ( $check eq "checked" ) {

    #appel script filterCDS.pl qui permet de filter les CDS
    my $filter_script = File::Spec->catfile( $dirScript, 'filterCDS.pl' );
    system("perl $filter_script $dirBlast $dirData $dirJob $plant");
}
else {
    system(
        'mv ' . $dirJob . 'sequenceUpload.fas ' . $dirJob . 'Sequences.fas' );
}

##Passage du multifasta -> fasta et appel au script Stemloop
open my $ENTREE_FH, '<', File::Spec->catfile( $dirJob, 'Sequences.fas' )
  or die "Impossible d'ouvrir le fichier d'entree  : $!";
my %tab = PipelineMiRNA::Utils::parse_multi_fasta($ENTREE_FH);
close $ENTREE_FH;

foreach my $name ( keys %tab ) {
    my $temp_file = File::Spec->catfile( $dirJob, 'tempFile.txt' );
    open( my $TEMPFILE_FH, '>', $temp_file )
      or die "Impossible d'ouvrir le fichier d'entree  : $!";
    system("chmod 777 $temp_file");
    print $TEMPFILE_FH $name . $tab{$name};
    my $name = substr $name, 1, -1;

    my $sequence_dir = File::Spec->catdir( $dirJob, $name );
    mkdir $sequence_dir;

    my $rnalfold_output = File::Spec->catfile( $sequence_dir, 'RNALfold.out' );
    my $rnalfold_cmd = "$rnalfold_bin < $temp_file > $rnalfold_output";
    if ($DEBUG) { print "$rnalfold_cmd\n" }
    system($rnalfold_cmd);

    ####conversion en format CT
  #    my $ct_file = File::Spec->catfile( $dirJob, $name, 'fichierOutB2ct.ct' );
  #    system("$b2ct_bin < $rnalfold_output > $ct_file");
  #    system("chmod 777 $ct_file");

    ## Appel de RNAstemloop
    my $rnastemloop_out =
      File::Spec->catfile( $sequence_dir, 'rnastemloop.out' );

    my $rnastemloop_cmd =
      "$rnastemploop_bin -i $rnalfold_output -o $rnastemloop_out";
    if ($DEBUG) { print "$rnastemloop_cmd\n" }
    system($rnastemloop_cmd);
    unlink $temp_file;
    process_RNAstemloop($rnastemloop_out);
    process_tests($dirJob);
}

sub process_RNAstemloop {
    my ($rnastemloop_out) = @_;
    my $sequence_dir = dirname($rnastemloop_out);

    my $rnaeval_out = File::Spec->catfile( $sequence_dir, 'rnaeval.out' );

    my $rnaeval_cmd = "$rnaeval_bin < $rnastemloop_out > $rnaeval_out";
    if ($DEBUG) { print "$rnaeval_cmd\n" }
    system($rnaeval_cmd);

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
                      File::Spec->catdir( $sequence_dir, $nameSeq );
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
    my $dirJob = shift;
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
                    my $b2ct_cmd =
                      "$b2ct_bin < $candidate_rnafold_out > $candidate_ct_file";
                    system($b2ct_cmd);
                    system("chmod 777 $candidate_ct_file");

#system('./test.sh '.$dirScript.' '.$dirJob.$dir.'/'.$file.'/ '.' > trash 2> errors');
                    my $varna_image =
                      File::Spec->catfile( $candidate_dir, 'image.png' );
                    system(
"/usr/bin/java -cp $varna_bin fr.orsay.lri.varna.applications.VARNAcmd -i $candidate_ct_file -o $varna_image > trash"
                    );

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
                        test_mfei( $candidate_ct_file, $sequence_dir, $file );
                    }
                    ####calcul p-value randfold
                    if ( $randfold eq "randfoldChecked" ) {
                        test_randfold( $candidate_dir, $seq_file );
                    }
                    ####calcul self-contain
                    if ( $SC eq "SCChecked" ) {
                        test_selfcontain( $candidate_dir, $seq_file );
                    }
                    ####### creation sequence boucle terminale masquee avec des N pour
                    ####### chaque sequence (repertoire ) et resultat alignement mirBASE
                    if ( $align eq "alignChecked" ) {
                        test_alignment( $candidate_dir, $seq_file );
                    }
                }
            }
        }
    }
    return 0;
}

sub test_mfei {

    #TODO: We do not need these three parameters
    #      but must restructure energie.pl first
    my ( $candidate_ct_file, $sequence_dir, $file ) = @_;
    my $energie_script = File::Spec->catfile( $dirScript, 'energie.pl' );
    my $energie_cmd =
      "perl $energie_script $candidate_ct_file $sequence_dir $file";

    #print "$energie_cmd\n";
    system($energie_cmd);
}

sub test_randfold {
    my ( $candidate_dir, $seq_file ) = @_;
    my $randfold_out = File::Spec->catfile( $candidate_dir, 'pvalue.txt' );
    system("$randfold_bin -d $seq_file 7 > $randfold_out");
    system("chmod 777 $randfold_out");
}

sub test_selfcontain {
    my ( $candidate_dir, $seq_file ) = @_;
    my $selfcontain_out =
      File::Spec->catfile( $candidate_dir, 'selfContain.txt' );
    system("python $selfcontain_bin -i $seq_file -n 100  > $selfcontain_out");
    system("chmod 777 $selfcontain_out");
}

sub test_alignment {
    my ( $candidate_dir, $seq_file ) = @_;
    my $seqN = File::Spec->catfile( $candidate_dir, 'seqWithN.txt' );
    open( my $SEQN_FH, '>>', $seqN )
      or die "Impossible d'ouvrir le fichier d'entree  : $!";
    open( my $SEQUENCE_FH, '<', $seq_file )
      or die "Impossible d'ouvrir le fichier d'entree  : $!";
    while ( my $line = <$SEQUENCE_FH> ) {
        if ( $line =~ /^>/ ) {
            print $SEQN_FH $line;
        }
        else {
            $line =~ s/U/T/gi;
            print $SEQN_FH $line;
        }
    }
    close $SEQN_FH;
    close $SEQUENCE_FH;
    my $exonerate_out = File::Spec->catfile( $candidate_dir, 'alignement.txt' );
    my $exonerate_cmd =
"$exonerate_bin -E --model affine:bestfit $mirbase_file $seqN -d $matrix_file --bestn 1 --score -3 -e -1 -o -1 > $exonerate_out";

    #print $exonerate_cmd;
    system($exonerate_cmd);
}
