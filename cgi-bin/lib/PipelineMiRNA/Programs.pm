package PipelineMiRNA::Programs;

# ABSTRACT: Executing external programs

use strict;
use warnings;

use File::Spec;
use File::Which;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Copy;
use Log::Message::Simple qw[msg error debug];

use PipelineMiRNA::Paths;
use PipelineMiRNA::Data;

my %programs_config = PipelineMiRNA::Paths->get_programs_config();

my $vienna_progs_dir = $programs_config{'vienna_package'};
my $rnafold_bin  = File::Spec->catfile( $vienna_progs_dir, 'RNAfold' );
my $rnalfold_bin = File::Spec->catfile( $vienna_progs_dir, 'RNALfold' );
my $rnaeval_bin  = File::Spec->catfile( $vienna_progs_dir, 'RNAeval' );
my $b2ct_bin     = File::Spec->catfile( $programs_config{'vienna_utils'}, 'b2ct' );

my $randfold_bin = $programs_config{'randfold'};
my $exonerate_dir = $programs_config{'exonerate'};

my $exonerate_bin = 'exonerate';
my $varna_bin        = $programs_config{'varna'};
my $rnastemploop_bin = $programs_config{'rnastemloop'};
my $blastx_bin = $programs_config{'blastx'};

my $miRdup_jar = $programs_config{'miRdup'};


=method list_programs

List all the binaries needed in the pipeline.

=cut

sub list_programs {
    my @args = @_;
    return (
        $rnafold_bin,  $rnalfold_bin, $rnaeval_bin,
        $randfold_bin, $rnastemploop_bin,
	$varna_bin,    
        $miRdup_jar
    );
}

=method list_unavailable_programs

List the missing binaries among those needed in the pipeline.

=cut

sub list_unavailable_programs {
    my @args       = @_;
    my @programs   = list_programs();
    my @unexisting = grep { !-f $_ && ! which($_) } @programs;
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

=method run_rnalfold_on_file

Run RNAfold on the given FASTA file
Return whether the output file exists.

=cut

sub run_rnalfold_on_file {
    my ( $input, $output ) = @_;
    my $options = '-L 400';
    my $rnalfold_cmd = "$rnalfold_bin $options < $input > $output";
    system($rnalfold_cmd);
    return ( -e $output );
}

=method run_rnalfold

Run RNAfold on the given FASTA sequence
Return the results of run_rnalfold_on_file

=cut

sub run_rnalfold {
    my ( @args ) = @_;
    my $sequence_name = shift @args;
    my $sequence      = shift @args;
    my $temp_file     = shift @args;
    my $output_file   = shift @args;
    open( my $TEMPFILE_FH, '>', $temp_file )
      or die "Error when opening tempfile -$temp_file-: $!";
    print $TEMPFILE_FH "$sequence_name\n$sequence";
    close $TEMPFILE_FH;
    my $result = run_rnalfold_on_file($temp_file, $output_file);
    unlink $TEMPFILE_FH;
    return $result;
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

=method run_slow_randfold

Run Randfold on the given file
Return whether the output file exists.

=cut

sub run_slow_randfold {
    my ( $input, $output ) = @_;
    my $nb_iterations = 7;
    my $randfold_cmd = "$randfold_bin -d $input $nb_iterations > $output";
    system($randfold_cmd);
    return ( -e $output );
}

=method run_randfold

Run Randfold on the given file
Return whether the output file exists.

=cut

sub run_randfold {
    my @args = @_;
    my $input_file = shift @args;
    my $output_file = shift @args;
    my $iterations  = shift @args;
    my ( $input, $output ) = @_;
    my $randfold_cmd = "$randfold_bin --iterations $iterations --fast --fasta $input_file > $output";
    debug($randfold_cmd, 1);
    system($randfold_cmd);
    return ( -e $output_file );
}

=method run_exonerate

Run Exonerate on the given file
Return whether the output file exists.

=cut

sub run_exonerate {
    my ( $input, $output ) = @_;
    my $matrix_file = PipelineMiRNA::Data::get_matrix_file();
    my $mirbase_file = PipelineMiRNA::Data::get_mirbase_file();
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

    my $exonerate_cmd = '';
    if( $exonerate_dir ){
        $exonerate_cmd .="cd $exonerate_dir && "
    }
    $exonerate_cmd =
        "$exonerate_bin "
      . "--exhaustive "             # Exhaustive alignment
      . "--model affine:bestfit "   # Best location alignment of the query onto the target
      . "$mirbase_file $input "
      . "--dnasubmat $matrix_file " # Substitution matrix to be used for DNA comparison
      . '--bestn 1 '                # Report the single best result for each query
      . '--score -3 '               # Overall score threshold
      . '--gapextend -1 '           # Gap extension penalty
      . '--gapopen -1 '             # Gap open penalty
      . '--showvulgar no '          # Do not show the alignments in "vulgar" format.
      . '--showalignment no '       # Do not show the human readable alignment
      . '--verbose 0 '              # Verbose to minimum
      . "--ryo '$output_fmt'"       # Using custom output format
      . "> $output  2> /dev/null";  # Stdout to file, stderr to /dev/null

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

 Usage : PipelineMiRNA::Programs::run_mirdup_prediction_on_sequence($sequence, $result_dir);
 Input : 
 Return: -

=cut

sub run_mirdup_prediction_on_sequence {
    my @args        = @_;
    my $sequence    = shift @args;
    my $output_dir  = shift @args;

    my $miRdup_model_name = PipelineMiRNA::Data::get_mirdup_model_name();
    my $miRdup_model_path = PipelineMiRNA::Data::get_mirdup_data_path();
    my $output_root = 'out';
    my $output =
        $output_root
      . '.generatedmirnas.'
      . $miRdup_model_name
      . '.miRdupOutput.txt';
    my $run_mirdup_cmd =
"java -jar $miRdup_jar -r $vienna_progs_dir/ -d $miRdup_model_name -predict -u $sequence -f $output_root > /dev/null 2>&1";
    chdir($miRdup_model_path) or die "$!";
    system($run_mirdup_cmd);

    my @old_files = glob "$miRdup_model_path/$output_root*";

    foreach my $old_file (@old_files) {
        move( $old_file, $output_dir )
          or die "Could not move $old_file to $output_dir: $!\n";
    }
    my $output_file = File::Spec->catfile( $output_dir, $output );
    return ( -e $output_file );
}

=method run_mirdup_prediction_on_sequence_file

Predict a miRNA given a pre-miRNA using MiRdup model.

 Usage : PipelineMiRNA::Programs::run_mirdup_prediction_on_sequence_file($prediction_source_file);
 Input : A prediction source file
 Return: Name of the output file

=cut

sub run_mirdup_prediction_on_sequence_file {
    my @args                   = @_;
    my $prediction_source_file = shift @args;

    my $miRdup_model_name = PipelineMiRNA::Data::get_mirdup_model_name();
    my $miRdup_model_path = PipelineMiRNA::Data::get_mirdup_data_path();

    my $run_mirdup_cmd =
"java -jar $miRdup_jar -r $vienna_progs_dir/ -d $miRdup_model_name -predict -i $prediction_source_file -f out  > /dev/null 2>&1";
    chdir($miRdup_model_path) or die "$!";
    system($run_mirdup_cmd);

    my $output_file =
        $prediction_source_file
      . '.miRdup.predictions.txt';
    return $output_file;
}

=method run_mirdup_validation_on_file

Validates a file with miRNA mature candidates using MiRdup.

 Usage : PipelineMiRNA::Programs::run_mirdup_validation_on_file($input_file);
 Input : File correctly formatted
 Return: the output file to parse

=cut

sub run_mirdup_validation_on_file {
    my @args                   = @_;
    my $validation_source_file = shift @args;

    my $miRdup_model_name = PipelineMiRNA::Data::get_mirdup_model_name();
    my $miRdup_model_path = PipelineMiRNA::Data::get_mirdup_data_path();

    my $run_mirdup_cmd =
"java -jar $miRdup_jar -r $vienna_progs_dir/ -c $miRdup_model_name -v $validation_source_file > /dev/null 2>&1";
    chdir($miRdup_model_path) or die "$!";
    system($run_mirdup_cmd);

    my $output_file =
        $validation_source_file . '.' . $miRdup_model_name . '.' . 'miRdup.tab.txt';
    return $output_file;
}

1;
