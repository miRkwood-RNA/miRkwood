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

