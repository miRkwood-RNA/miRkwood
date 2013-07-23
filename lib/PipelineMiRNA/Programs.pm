package PipelineMiRNA::Programs;

use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Copy;

my $module_path = abs_path(dirname(__FILE__));
my $rootdir = File::Spec->catdir($module_path, '..', '..');

my $dirProgs = File::Spec->catdir( $rootdir, 'programs' );

my $vienna_dir   = File::Spec->catfile( $dirProgs,   'ViennaRNA-2.1.2' );
my $vienna_progs_dir = File::Spec->catdir( $vienna_dir, 'Progs');
my $rnafold_bin  = File::Spec->catfile( $vienna_progs_dir, 'RNAfold' );
my $rnalfold_bin = File::Spec->catfile( $vienna_progs_dir, 'RNALfold' );
my $rnaeval_bin  = File::Spec->catfile( $vienna_progs_dir, 'RNAeval' );
my $b2ct_bin     = File::Spec->catfile( $vienna_dir, 'Utils', 'b2ct' );

my $randfold_bin = File::Spec->catfile( $dirProgs, 'randfold-2.0', 'randfold' );
my $selfcontain_bin =
  File::Spec->catfile( $dirProgs, 'selfcontain_unix', 'selfcontain.py' );
my $exonerate_bin =
  File::Spec->catfile( $dirProgs, 'exonerate-2.2.0-i386', 'bin', 'exonerate' );
my $varna_bin        = File::Spec->catfile( $dirProgs, 'VARNAv3-9.jar' );
my $rnastemploop_bin = File::Spec->catfile( $dirProgs, 'RNAstemloop' );
my $blastx_bin = 'blastx';
#File::Spec->catfile( $dirProgs, 'blastx' );

my $miRdup_jar = File::Spec->catfile( $dirProgs, 'miRdup_1.1'   , 'miRdup.jar' );

## Data ##
my $dirData = File::Spec->catdir( $rootdir, 'data' );    # chemin séquence
my $mirbase_file = File::Spec->catfile( $dirData, 'MirbaseFile.txt' );
my $matrix_file  = File::Spec->catfile( $dirData, 'matrix' );
my $miRdup_model_path =File::Spec->catdir( $dirData, 'mirdup');
my $miRdup_model_name = 'MirbaseFile.model';

sub run_varna {

    # Run VARNA on the given CT file
    my ( $candidate_ct_file, $varna_image ) = @_;
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
    return ( -e $ct_file );
}

sub run_rnalfold {
    my ( $input, $output ) = @_;
    my $options = '-L 400';
    my $rnalfold_cmd = "$rnalfold_bin $options < $input > $output";
    system($rnalfold_cmd);
    return ( -e $output );
}

sub run_rnaeval {
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
    my $selfcontain_cmd =
      "python $selfcontain_bin -i $input -n $num_contexts  > $output";
    system($selfcontain_cmd);
    return ( -e $output );
}

sub run_exonerate {
    my ( $input, $output ) = @_;
    my $output_fmt = 'name : %ti\n  begin: %tab\n  end  : %tae\n  score: %s\n  seq  : %tas\n';
    my $exonerate_cmd =
        "$exonerate_bin " . "-E "
      . "--model affine:bestfit $mirbase_file $input "
      . "-d $matrix_file "
      . '--bestn 1 '
      . '--score -3 '
      . '-e -1 -o -1 '
      . '--showvulgar no --showalignment no --verbose 0 '
      . "--ryo '$output_fmt'"
      . "> $output  2> /dev/null";
    system($exonerate_cmd);
    return ( -e $output );
}

sub run_rnastemloop {
    my ( $input, $output_stemloop, $output_optimal ) = @_;
    my $rnastemloop_cmd = "$rnastemploop_bin -i $input --output-stemloop $output_stemloop --output-optimal $output_optimal";
    system($rnastemloop_cmd);
    return ( -e $output_stemloop && -e $output_optimal);
}

sub run_blast {
    my ( $query, $database, $blastx_options, $output ) = @_;
    my $blastx_cmd =
        "$blastx_bin "
      . "-query $query "
      . "-db $database "
      . "$blastx_options "
      . "-out $output";
    system($blastx_cmd);
    return ( -e $output );
}

=method train_mirdup

Train a MiRdup model using the given data

 Usage : PipelineMiRNA::Programs::train_mirdup($matures_miRNA, $hairpins_precursors);
 Input : Matures miRNAs file, in FASTA format
         Hairpins precursors file, in FASTA format
 Return: -

=cut

sub train_mirdup {
    my @args                = @_;
    my $matures_miRNA       = shift @args;
    my $hairpins_precursors = shift @args;
    my $run_mirdup_cmd =
"java -Xms500m -Xmx1500m -jar $miRdup_jar -r $vienna_progs_dir -m $matures_miRNA -h $hairpins_precursors";
    system($run_mirdup_cmd);
    return;
}

=method run_mirdup_prediction_on_sequence

Predict a miRNA given a pre-miRNA using MiRdup model.

 Usage : PipelineMiRNA::Programs::run_mirdup_prediction_on_sequence($sequence, $result_dir, $name);
 Input : 
 Return: -

=cut

sub run_mirdup_prediction_on_sequence {
    my @args        = @_;
    my $sequence    = shift @args;
    my $output_dir  = shift @args;
    my $output_name = shift @args;
    my $output =
        $output_name
      . '.generatedmirnas.'
      . $miRdup_model_name
      . '.miRdupOutput.txt';
    my $run_mirdup_cmd =
"java -jar $miRdup_jar -r $vienna_progs_dir/ -d $miRdup_model_name -predict -u $sequence -f $output_name ";
    chdir($miRdup_model_path) or die "$!";
    system($run_mirdup_cmd);

    my @old_files = glob "$miRdup_model_path/$output_name*";

    foreach my $old_file (@old_files) {
        move( $old_file, $output_dir )
          or die "Could not move $old_file to $output_dir: $!\n";
    }
    my $output_file = File::Spec->catfile( $output_dir, $output );
    return ( -e $output_file );
}

1;
