package miRkwood::Programs;

# ABSTRACT: Executing external programs

use strict;
use warnings;

use File::Spec;
use File::Which;
use File::Temp qw/ tempdir /;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Copy;
use Log::Message::Simple qw[msg error debug];

use miRkwood::Paths;
use miRkwood::Data;

our %programs_paths;

sub init_programs{
    my @args = @_;
    my %programs_config = miRkwood::PROGRAMS_CONFIG();
    my $mirkwood_local_programs = miRkwood::Paths::get_local_programs_path();
    while (my ($key, $value) = each %programs_config) {
        $value =~ s/\${local_programs}/$mirkwood_local_programs/g;
        $programs_paths{$key} = $value;
    }
    $programs_paths{'vienna_progs'} = get_Vienna_program_path($programs_config{'rnafold'});
    return 1;
}

=method get_Vienna_program_path

Infers the Vienna programs path based on the configuration value

=cut

sub get_Vienna_program_path {
    my @args = @_;
    my $rnafold_exe = shift @args;
    my $folder;
    if (-f $rnafold_exe){
        $folder = $rnafold_exe;
    } else {
        $folder = which($rnafold_exe);
    }
    if ( ! defined($folder) || $folder eq '' ){
        die( 'Cannot find required third-party software: rnafold.' );
    }
    return File::Basename::dirname($folder);
}

=method list_programs

List all the binaries needed in the pipeline.

=cut

sub list_programs {
    my @args = @_;
    my %progs = %programs_paths;
    delete $progs{'vienna_progs'};
    delete $progs{'mirdup'};
    delete $progs{'varna'};
    delete $progs{'piccolo'};

    # These programs are used only by abinitio pipelines and
    # for optional features.
    delete $progs{'tRNAscanSE'};
    delete $progs{'blastx'};
    delete $progs{'rnammer'};
    delete $progs{'rnashuffles'};
    return values %progs;
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

=method run_varna_on_ct_file

Run VARNA on the given CT file
Return whether the output file exists.

=cut

sub run_varna_on_ct_file {
    my ( $candidate_ct_file, $varna_image ) = @_;
    my $varna_cmd =
"/usr/bin/java -cp $programs_paths{'varna'} fr.orsay.lri.varna.applications.VARNAcmd -titleSize 0 -i $candidate_ct_file -o $varna_image > /dev/null 2>&1";
    system($varna_cmd);
    return ( -e $varna_image );
}

=method run_varna_on_structure

Run VARNA on the given structure
Return whether the output file exists.

=cut

sub run_varna_on_structure {
    my ( $sequence, $structure, $varna_image ) = @_;
    if ( -f $programs_paths{'varna'} ){
        my $varna_cmd =
    "/usr/bin/java -cp $programs_paths{'varna'} fr.orsay.lri.varna.applications.VARNAcmd -titleSize 0 -sequenceDBN '$sequence' -structureDBN '$structure' -o $varna_image > /dev/null 2>&1";
        system($varna_cmd);
        return ( -e $varna_image );
    }
    else{
        debug('[WARNING] VARNA is not installed. Cannot create the image with VARNA.', miRkwood->DEBUG());
    }
}

=method convert_to_ct

Convert to CT format using b2ct
Return whether the output file exists.

=cut

sub convert_to_ct {
    my ( $rnafold_out, $ct_file ) = @_;
    my $b2ct_cmd = "$programs_paths{'b2ct'} < $rnafold_out > $ct_file";
    system($b2ct_cmd);
    chmod 0777, $ct_file;
    return ( -e $ct_file );
}

=method run_rnalfold_on_file

Run RNAlfold on the given FASTA file
Return whether the output file exists.

=cut

sub run_rnalfold_on_file {
    my ( $input, $output ) = @_;
    my $options = '-L 400';
    my $rnalfold_cmd = "$programs_paths{'rnalfold'} $options < $input > $output";
    #~ debug( '     ' . $rnalfold_cmd, miRkwood->DEBUG());
    debug( "       Running RNAlfold on $input" . ' [' . localtime() . ']', miRkwood->DEBUG());
    system($rnalfold_cmd);
    return ( -e $output );
}

=method run_rnafold_on_file

Run RNAfold on the given FASTA file
Return whether the output file exists.

=cut

sub run_rnafold_on_file {
    my ( $input, $output ) = @_;
    my $rnafold_cmd = "$programs_paths{'rnafold'} --noPS < $input > $output";
    #~ debug( '        ' . $rnafold_cmd, miRkwood->DEBUG());
    system($rnafold_cmd);
    return ( -e $output );
}

=method run_rnalfold

Run RNAlfold on the given FASTA sequence
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
    print {$TEMPFILE_FH} ">$sequence_name\n$sequence"
      or die "Error when writing in $temp_file: $!";
    close $TEMPFILE_FH
      or die "Error when closing $temp_file: $!";
    my $result = run_rnalfold_on_file($temp_file, $output_file);
    unlink $temp_file;
    return $result;
}

=method run_rnafold

Run RNAfold on the given FASTA sequence
Return the results of run_rnafold_on_file

=cut

sub run_rnafold {
    my ( @args ) = @_;
    my $sequence_name = shift @args;
    my $sequence      = shift @args;
    my $temp_file     = shift @args;
    my $output_file   = shift @args;
    open( my $TEMPFILE_FH, '>', $temp_file )
      or die "Error when opening tempfile -$temp_file-: $!";
    print {$TEMPFILE_FH} ">$sequence_name\n$sequence"
      or die "Error when writing in $temp_file: $!";
    close $TEMPFILE_FH
      or die "Error when closing $temp_file: $!";
    my $result = run_rnafold_on_file($temp_file, $output_file);
    unlink $temp_file;
    return $result;
}

=method run_rnaeval

Run RNAeval on the file
Return whether the output file exists.

=cut

sub run_rnaeval {
    my ( $input, $output ) = @_;
    my $rnaeval_cmd = "$programs_paths{'rnaeval'} < $input > $output";
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
    my $randfold_cmd = "$programs_paths{'rnashuffles'} -d $input $nb_iterations > $output";
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
    my $randfold_cmd = "$programs_paths{'rnashuffles'} --iterations $iterations --fast --fasta $input_file > $output_file";
    debug( '        ' . $randfold_cmd, miRkwood->DEBUG());
    system($randfold_cmd);
    return ( -e $output_file );
}

=method run_exonerate

Run Exonerate on the given file
Return whether the output file exists.

=cut

sub run_exonerate {
    my ( $input, $output ) = @_;
    my $matrix_file = miRkwood::Data::get_matrix_file();
    my $mirbase_file = miRkwood::Data::get_mirbase_file();
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

    my $exonerate_cmd = qw{};
    $exonerate_cmd =
        "$programs_paths{'exonerate'} "
      . '--exhaustive '             # Exhaustive alignment
      . '--model affine:bestfit '   # Best location alignment of the query onto the target
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

    debug( '        ' . $exonerate_cmd, miRkwood->DEBUG());
    system($exonerate_cmd);
    return ( -e $output );
}

=method run_piccolo

  Run piccolo on the given file
  Return whether the output file exists.
  piccolo output format is the same as exonerate, so
  exonerate output parser can be used for piccolo output file.

=cut
sub run_piccolo {
    my ( $input, $output ) = @_;
    my $mirbase_file = miRkwood::Data::get_mirbase_file();
    if ( ! -f $programs_paths{'piccolo'} ) {
        debug('[WARNING] piccolo is not installed. Cannot align against miRBase.', miRkwood->DEBUG());
    }
    else {
        my $cmd = "$programs_paths{'piccolo'} -r $mirbase_file -i $input --half > $output 2> /dev/null";
        #~ debug( '          ' . $cmd, miRkwood->DEBUG());
        debug( "         Running piccolo on $input" . ' [' . localtime() . ']', miRkwood->DEBUG());
        system($cmd);
    }
    return ( -e $output );
}

=method run_rnastemloop

Run RNAstemloop on the given file
Return whether the output files exist.

=cut

sub run_rnastemloop {
    my ( $input, $output_stemloop, $output_optimal ) = @_;
    my $rnastemloop_cmd = "$programs_paths{'rnastemloop'} -i $input -s $output_stemloop -o $output_optimal";
    #~ debug( '     ' . $rnastemloop_cmd, miRkwood->DEBUG());
    debug( "       Running RNAstemloop on $input" . ' [' . localtime() . ']', miRkwood->DEBUG());
    system($rnastemloop_cmd);
    return ( -e $output_stemloop && -e $output_optimal);
}

=method run_blast

Run BLAST.

Usage:
miRkwood::Programs::run_blast($query_file, $blast_database_file,
                                   $blastx_options, $blast_output_file)

Return whether the output file exists.

=cut

sub run_blast {
    my ( $query, $database, $blastx_options, $output ) = @_;
    my $blastx_cmd =
        "$programs_paths{'blastx'} "
      . "-query $query "
      . "-db $database "
      . "$blastx_options "
      . "-out $output";
    debug($blastx_cmd, miRkwood->DEBUG());
    system($blastx_cmd);
    return ( -e $output );
}

=method train_mirdup

Train a MiRdup model using the given data

 Usage : miRkwood::Programs::train_mirdup($matures_miRNA, $hairpins_precursors);
 Input : Matures miRNAs file, in FASTA format
         Hairpins precursors file, in FASTA format
 Return: -

=cut

sub train_mirdup {
    my @args                = @_;
    my $matures_miRNA       = shift @args;
    my $hairpins_precursors = shift @args;
    my $run_mirdup_cmd =
"java -Xms500m -Xmx1500m -jar $programs_paths{'mirdup'} -r $programs_paths{'vienna_progs'} -m $matures_miRNA -h $hairpins_precursors";
    system($run_mirdup_cmd);
    return;
}

=method run_mirdup_prediction_on_sequence

Predict a miRNA given a pre-miRNA using MiRdup model.

 Usage : miRkwood::Programs::run_mirdup_prediction_on_sequence($sequence, $result_dir);
 Input : 
 Return: -

=cut

sub run_mirdup_prediction_on_sequence {
    my @args        = @_;
    my $sequence    = shift @args;
    my $output_dir  = shift @args;

    my $miRdup_model_name = miRkwood::Data::get_mirdup_model_name();
    my $miRdup_model_path = miRkwood::Data::get_mirdup_data_path();
    my $output_root = 'out';
    my $output =
        $output_root
      . '.generatedmirnas.'
      . $miRdup_model_name
      . '.miRdupOutput.txt';

    my $run_mirdup_cmd = "cd $miRdup_model_path && "
        . "java -jar $programs_paths{'mirdup'}"
        . " -r $programs_paths{'vienna_progs'}/"
        . " -d $miRdup_model_name"
        . ' -predict'
        . " -u $sequence"
        . " -f $output_root"
        . " > /dev/null 2>&1";
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

 Usage : miRkwood::Programs::run_mirdup_prediction_on_sequence_file($prediction_source_file);
 Input : A prediction source file
 Return: Name of the output file

=cut

sub run_mirdup_prediction_on_sequence_file {
    my @args                   = @_;
    my $prediction_source_file = shift @args;

    my $miRdup_model_name = miRkwood::Data::get_mirdup_model_name();
    my $miRdup_model_path = miRkwood::Data::get_mirdup_data_path();

    my $run_mirdup_cmd = "cd $miRdup_model_path && "
        . "java -jar $programs_paths{'mirdup'}"
        . " -r $programs_paths{'vienna_progs'}/"
        . " -d $miRdup_model_name"
        . ' -predict'
        . " -i $prediction_source_file"
        . ' -f out'
        . "  > /dev/null 2>&1";
    system($run_mirdup_cmd);

    my $output_file =
        $prediction_source_file
      . '.miRdup.predictions.txt';
    return $output_file;
}

=method run_mirdup_validation_on_file

Validates a file with miRNA mature candidates using MiRdup.

 Usage : miRkwood::Programs::run_mirdup_validation_on_file($input_file);
 Input : File correctly formatted
 Return: the output file to parse

=cut

sub run_mirdup_validation_on_file {
    my @args                   = @_;
    my $validation_source_file = shift @args;

    my $miRdup_model_name = miRkwood::Data::get_mirdup_model_name();
    my $miRdup_model_path = miRkwood::Data::get_mirdup_data_path();

    if ( -f $programs_paths{'mirdup'} ){

        my $run_mirdup_cmd = "cd $miRdup_model_path && "
            . "java -jar $programs_paths{'mirdup'}"
            . " -r $programs_paths{'vienna_progs'}/"
            . " -c $miRdup_model_name"
            . " -v $validation_source_file"
            . " > /dev/null";
        system($run_mirdup_cmd);

        my @useless_created_by_mirdup_files = ("$validation_source_file.arff",
                "$validation_source_file.arff",
                "$validation_source_file.folded",
                "$validation_source_file.$miRdup_model_name.miRdup.txt",
                "$validation_source_file.$miRdup_model_name.miRdupOutput.txt");

        foreach my $file ( @useless_created_by_mirdup_files ){
            if ( -e $file ){
                unlink $file;
            }
        }

    }
    else{
        debug('[WARNING] miRdup is not installed. Cannot validate the miRNA with miRdup.', miRkwood->DEBUG());
    }

    my $output_file =
        $validation_source_file . '.' . $miRdup_model_name . '.' . 'miRdup.tab.txt';
    return $output_file;
}


=method run_tRNAscanSE_on_file

Run run_tRNAscanSE on the given file
Return whether the output file exists.

=cut

sub run_tRNAscanSE_on_file {
    my ( $input, $output ) = @_;
    my $tRNAscanSE_cmd = qw{};
    my $tRNAscanSE_bin = $programs_paths{'tRNAscanSE'};
    $tRNAscanSE_cmd =
        "PERL5LIB=$tRNAscanSE_bin: "
      . "$tRNAscanSE_bin/tRNAscan-SE "
      . "$input "
      . '--quiet '
      . '--brief '
      . '--forceow '
      . "--output $output";
    debug("$tRNAscanSE_cmd", miRkwood->DEBUG());
    system($tRNAscanSE_cmd);
    return ( -e $output );
}

=method run_rnammer_on_file

Run RNAmmer on the given file
Return whether the output file exists.

=cut

sub run_rnammer_on_file {
    my ( $input, $kingdom, $output ) = @_;
    my $tmp_dir = tempdir( CLEANUP => 1 );
    my $rnammer_cmd = qw{};
    $rnammer_cmd =
      "$programs_paths{'rnammer'} "
      . "-T $tmp_dir "          # Temporary directory
      . "-S $kingdom "          # Kingdom
      . "-m lsu,ssu,tsu "       # Molecule types
      . "--gff $output "        # GFF output
      . "< $input";
    debug( "$rnammer_cmd", miRkwood->DEBUG() );
    system($rnammer_cmd);
    return ( -e $output );
}


1;
