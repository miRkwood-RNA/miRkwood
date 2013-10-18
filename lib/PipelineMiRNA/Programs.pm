package PipelineMiRNA::Programs;

# ABSTRACT: Executing external programs

use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Copy;
use Log::Message::Simple qw[msg error debug];

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
my $exonerate_dir =
  File::Spec->catdir( $dirProgs, 'exonerate-2.2.0-i386', 'bin' );
my $exonerate_bin = 'exonerate';
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
my $miRdup_model_name = 'plant.model';

=method list_programs

List all the binaries needed in the pipeline.

=cut

sub list_programs {
    my @args = @_;
    return (
        $rnafold_bin,  $rnalfold_bin, $rnaeval_bin,
        $randfold_bin, $varna_bin,    $rnastemploop_bin
    );
}

=method list_unavailable_programs

List the missing binaries among those needed in the pipeline.

=cut

sub list_unavailable_programs {
    my @args       = @_;
    my @programs   = list_programs();
    my @unexisting = grep { !-f $_ } @programs;
    return @unexisting;
}

=method run_varna

Run VARNA on the given CT file
Return whether the output file exists.

=cut

sub run_varna {
    my ( $candidate_ct_file, $varna_image ) = @_;
    my $varna_cmd =
"/usr/bin/java -cp $varna_bin fr.orsay.lri.varna.applications.VARNAcmd -titleSize 0 -i $candidate_ct_file -o $varna_image > /dev/null 2>&1";
    system($varna_cmd);
    return ( -e $varna_image );
}

=method convert_to_ct

Convert to CT format using b2ct
Return whether the output file exists.

=cut

sub convert_to_ct {
    my ( $rnafold_out, $ct_file ) = @_;
    my $b2ct_cmd = "$b2ct_bin < $rnafold_out > $ct_file";
    system($b2ct_cmd);
    chmod 0777, $ct_file;
    return ( -e $ct_file );
}

=method run_rnalfold

Run RNAfold on the given FASTA file
Return whether the output file exists.

=cut

sub run_rnalfold {
    my ( $input, $output ) = @_;
    my $options = '-L 400';
    my $rnalfold_cmd = "$rnalfold_bin $options < $input > $output";
    system($rnalfold_cmd);
    return ( -e $output );
}

=method run_rnaeval

Run RNAeval on the file
Return whether the output file exists.

=cut

sub run_rnaeval {
    my ( $input, $output ) = @_;
    my $rnaeval_cmd = "$rnaeval_bin < $input > $output";
    system($rnaeval_cmd);
    return ( -e $output );
}

=method run_randfold

Run Randfold on the given file
Return whether the output file exists.

=cut

sub run_randfold {
    my ( $input, $output ) = @_;
    my $randfold_cmd = "$randfold_bin -d $input 7 > $output";
    system($randfold_cmd);
    return ( -e $output );
}

=method run_selfcontain

Run SelfContain on the given sequence file
Return whether the output file exists.

=cut

sub run_selfcontain {
    my ( $input, $output ) = @_;
    my $num_contexts = 100;
    my $selfcontain_cmd =
      "python $selfcontain_bin -i $input -n $num_contexts  > $output";
    debug($selfcontain_cmd, 1);
    system($selfcontain_cmd);
    return ( -e $output );
}

=method run_exonerate

Run Exonerate on the given file
Return whether the output file exists.

=cut

sub run_exonerate {
    my ( $input, $output ) = @_;
    my $output_fmt = '- name : %qi\n'
                    .'  begin_target  : %tab\n'
                    .'  end_target    : %tae\n'
                    .'  strand_target : %tS\n'
                    .'  begin_query   : %qab\n'
                    .'  end_query     : %qae\n'
                    .'  def_query     : %qd\n'
                    .'  seq_query     : %qas\n'
                    .'  score: %s\n'
                    .'  alignment: |{\n'
                    .'    %Pqs %Pts %Pl}\n';
    my $exonerate_cmd =
        "cd $exonerate_dir && "
      . "$exonerate_bin " . "-E "
      . "--model affine:bestfit $mirbase_file $input "
      . "-d $matrix_file "
      . '--bestn 1 '
      . '--score -3 '
      . '-e -1 -o -1 '
      . '--showvulgar no --showalignment no --verbose 0 '
      . "--ryo '$output_fmt'"
      . "> $output  2> /dev/null";

    debug($exonerate_cmd, 1);
    system($exonerate_cmd);
    return ( -e $output );
}

=method run_rnastemloop

Run RNAstemloop on the given file
Return whether the output files exist.

=cut

sub run_rnastemloop {
    my ( $input, $output_stemloop, $output_optimal ) = @_;
    my $rnastemloop_cmd = "$rnastemploop_bin -i $input --output-stemloop $output_stemloop --output-optimal $output_optimal";
    system($rnastemloop_cmd);
    return ( -e $output_stemloop && -e $output_optimal);
}

=method run_blast

Run BLAST.

Usage:
PipelineMiRNA::Programs::run_blast($query_file, $blast_database_file,
                                   $blastx_options, $blast_output_file)

Return whether the output file exists.

=cut

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

sub run_mirdup_prediction_on_sequence_file {
    my @args                   = @_;
    my $prediction_source_file = shift @args;

    my $run_mirdup_cmd =
"java -jar $miRdup_jar -r $vienna_progs_dir/ -d $miRdup_model_name -predict -i $prediction_source_file -f out  > /dev/null 2>&1";
    chdir($miRdup_model_path) or die "$!";
    system($run_mirdup_cmd);

    my $output_file =
        $prediction_source_file
      . '.miRdup.predictions.txt';
    return $output_file;
}

1;