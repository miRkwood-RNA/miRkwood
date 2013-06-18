package PipelineMiRNA::Programs;

use strict;
use warnings;

use PipelineMiRNA::Energie;

my $rootdir = '/home/jeanfred/Tuile/pipelineMiRNA-web/';

my $dirProgs = File::Spec->catdir( $rootdir, 'programs' );

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
my $dirData = File::Spec->catdir( $rootdir, 'data' );    # chemin séquence
my $mirbase_file = File::Spec->catfile( $dirData, 'MirbaseFile.txt' );
my $matrix_file  = File::Spec->catfile( $dirData, 'matrix' );

sub run_varna {
    # Run VARNA on the given CT file
    my ($candidate_ct_file, $varna_image) = @_;
    my $varna_cmd =
"/usr/bin/java -cp $varna_bin fr.orsay.lri.varna.applications.VARNAcmd -i $candidate_ct_file -o $varna_image > /dev/null 2>&1";
    system($varna_cmd);
    return ( -e $varna_image );
}

sub convert_to_ct {

    # Convert to CT format using b2ct
    my ( $rnafold_out, $ct_file ) = @_;
    my $b2ct_cmd = "$b2ct_bin < $rnafold_out > $ct_file";
    system($b2ct_cmd);
    chmod 0777, $ct_file;
    return (-e $ct_file);
}

sub run_rnalfold {
    my ( $input, $output ) = @_;
    my $rnalfold_cmd = "$rnalfold_bin < $input > $output";
    system($rnalfold_cmd);
    return ( -e $output );
}

sub run_rnaeval{
    my ( $input, $output ) = @_;
    my $rnaeval_cmd = "$rnaeval_bin < $input > $output";
    system($rnaeval_cmd);
    return ( -e $output );
}

sub run_randfold {
    my ( $input, $output ) = @_;
    my $randfold_cmd = "$randfold_bin -d $input 7 > $output";
    system($randfold_cmd);
    return ( -e $output );
}

sub run_selfcontain {
    my ( $input, $output ) = @_;
    my $num_contexts = 100;
    my $selfcontain_cmd = "python $selfcontain_bin -i $input -n $num_contexts  > $output";
    system($selfcontain_cmd);
    return ( -e $output );
}

sub run_exonerate {
    my ( $input, $output ) = @_;
    my $exonerate_cmd =
        "$exonerate_bin " . "-E "
      . "--model affine:bestfit $mirbase_file $input "
      . "-d $matrix_file "
      . "--bestn 1 "
      . "--score -3 "
      . "-e -1 -o -1 "
      . "> $output ";
    system($exonerate_cmd);
    return ( -e $output );
}

sub run_rnastemloop {
    my ( $input, $output ) = @_;
    my $rnastemloop_cmd = "$rnastemploop_bin -i $input -o $output";
    system($rnastemloop_cmd);
    return ( -e $output );
}

1;
