package PipelineMiRNA::MainPipeline;

# ABSTRACT: The actual pipeline

use strict;
use warnings;

use CGI::Carp qw(fatalsToBrowser);
use File::Path 'rmtree';
use File::Basename;
use Cwd qw( abs_path );
use File::Copy;
use PipelineMiRNA;
use PipelineMiRNA::Paths;
use PipelineMiRNA::Utils;
use PipelineMiRNA::MiRdup;
use PipelineMiRNA::Parsers;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Candidate;
use PipelineMiRNA::Components;
use PipelineMiRNA::PosterioriTests;
use Log::Message::Simple qw[msg error debug];

### Data ##
my $dirData = PipelineMiRNA::Paths->get_data_path();

=method write_config

Write the run options to the job configuration file.

=cut

sub write_config {
	my ( $mfe, $randfold, $align, $run_options_file ) = @_;
	my $run_options = PipelineMiRNA->CONFIG();
	$run_options->param( "options.mfe",      $mfe );
	$run_options->param( "options.randfold", $randfold );
	$run_options->param( "options.align",    $align );
	PipelineMiRNA->CONFIG($run_options);
}

=method main_entry

Run the pipeline.

 Usage : PipelineMiRNA::MainPipeline::main_entry( $ifilterCDS, $imfei, $irandfold, $ialign, $idirJob, $iplant );
 Input : Booleans for pipeline options, job directory, plant identifier
 Return: -

=cut

sub main_entry {
	my ( $check, $mfe, $randfold, $align, $job_dir, $plant ) = @_;
	my $debug    = 1;
	my $log_file = File::Spec->catfile( $job_dir, 'log.log' );
	local $Log::Message::Simple::DEBUG_FH = PipelineMiRNA->LOGFH($log_file);

	my $run_options_file = PipelineMiRNA::Paths->get_job_config_path($job_dir);
	PipelineMiRNA->CONFIG_FILE($run_options_file);
	write_config( $mfe, $randfold, $align, $run_options_file );

	debug( 'BEGIN execute_scripts', $debug );
	my $sequences_input = File::Spec->catfile( $job_dir, 'Sequences.fas' );
	if ( $check eq 'checked' ) {
		debug( 'FilteringCDS', $debug );

		#Filtering CDS
		PipelineMiRNA::Components::filter_CDS( $dirData, $job_dir, $plant );
	}
	else {
		my $sequence_uploaded =
		  File::Spec->catfile( $job_dir, 'sequenceUpload.fas' );
		debug( "Moving file $sequence_uploaded to $sequences_input", $debug );
		File::Copy::move( $sequence_uploaded, $sequences_input );
	}

	##Passage du multifasta -> fasta et appel au script Stemloop
	debug( "Opening multifasta $sequences_input", $debug );
	open my $ENTREE_FH, '<', $sequences_input
	  or die "Error when opening sequences -$sequences_input-: $!";
	debug( "Calling parse_multi_fasta() on $sequences_input", $debug );
	my %tab = PipelineMiRNA::Utils::parse_multi_fasta($ENTREE_FH);
	close $ENTREE_FH;

	debug( 'Iterating over names', $debug );

    my $sequence_dir_name = 0;

    while ( my ($name, $sequence) = each %tab) {
        debug( "Considering sequence $sequence_dir_name: $name", $debug );
        $sequence_dir_name++;
        my $sequence_dir = File::Spec->catdir( $job_dir, $sequence_dir_name );
        mkdir $sequence_dir;
        process_sequence( $sequence_dir, $name, $sequence );
    }
    process_tests( $job_dir );
    return;
}

=method process_sequence

Process a single sequence

 Usage : process_sequence( $sequence_dir, $name, $sequence );
 Return: -

=cut

sub process_sequence {
    my @args = @_;
    my $sequence_dir = shift @args;
    my $name = shift @args;
    my $sequence = shift @args;

    ## Running RNALfold
    debug( 'Running RNALfold', 1 );
    my $rnalfold_output =
      File::Spec->catfile( $sequence_dir, 'RNALfold.out' );

	my $temp_file = File::Spec->catfile( $sequence_dir, 'tempFile.txt' );
    PipelineMiRNA::Programs::run_rnalfold( $name, $sequence, $temp_file, $rnalfold_output )
     or die("Problem when running RNALfold: $!");

	## Running RNAstemloop
	my $rnastemloop_out_optimal =
	  File::Spec->catfile( $sequence_dir, 'rnastemloop_optimal.out' );
	my $rnastemloop_out_stemloop =
	  File::Spec->catfile( $sequence_dir, 'rnastemloop_stemloop.out' );
	debug( "Running RNAstemloop on $rnalfold_output", 1 );
	PipelineMiRNA::Programs::run_rnastemloop( $rnalfold_output,
		$rnastemloop_out_optimal, $rnastemloop_out_stemloop )
	  or die("Problem when running RNAstemloop");

	run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_optimal,  'optimal' );
	my $rnaeval_out =
	   run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_stemloop, 'stemloop' );

	open( my $STEM_FH, '<', $rnastemloop_out_stemloop ) or die $!;
	open( my $EVAL_FH, '<', $rnaeval_out ) or die $!;
	my $res = process_RNAstemloop( $sequence_dir, 'stemloop', $STEM_FH, $EVAL_FH );
	my @hash = @{$res};
	close($STEM_FH);
	close($EVAL_FH);
    my %newHash = treat_candidates( \@hash );
    create_directories( \%newHash, $sequence_dir );
}

=method run_RNAeval_on_RNAstemloop_output


=cut

sub run_RNAeval_on_RNAstemloop_output {
	my ( $rnastemloop_out, $suffix ) = @_;
	my $current_sequence_dir = dirname($rnastemloop_out);
	debug( "Processing RNAstemloop output for $suffix $rnastemloop_out", 1 );
	my $rnaeval_out =
	  File::Spec->catfile( $current_sequence_dir, "rnaeval_$suffix.out" );

	debug( "Running RNAeval in $rnaeval_out", 1 );
	PipelineMiRNA::Programs::run_rnaeval( $rnastemloop_out, $rnaeval_out )
	  or die("Problem when running RNAeval");

	return $rnaeval_out;
}

=method process_RNAstemloop

Process the outputs of RNAstemloop + RNAeval
Writes the sequence on disk (seq.txt) and outRNAFold.txt
(for a given suffix)

=cut

sub process_RNAstemloop {
	my @args           = @_;
	my ($sequence_dir) = shift @args;
	my ($suffix)       = shift @args;
	my ($STEM_FH)      = shift @args;
	my ($EVAL_FH)      = shift @args;
	my $index          = 0;
	my $line2;
	my ( $nameSeq, $dna, $Vienna );
	my @hash = ();

	while ( my $line = <$STEM_FH> ) {

		if ( ( $line =~ /^>(.*)/ ) ) {    # nom sequence
			$nameSeq = $1;
		}
		elsif ( ( $line =~ /^[a-zA-Z]/ ) ) { # récupération de la sequence adn
			$dna = substr $line, 0, -1;
			$line2 = substr( <$EVAL_FH>, 0, -1 );    # the sequence as well

			if ( $dna ne $line2 ) {
				                                  # Should not happen
			}

		}
		elsif ( ( $line =~ /(.*)/ ) ) {
			$Vienna = $1;
			$line2  = <$EVAL_FH>;    # the structure as well, and the energy
			if ( my ( $structure, $energy ) =
				PipelineMiRNA::Parsers::parse_Vienna_line($line2) )
			{                     # We have a structure
				if ( $Vienna ne $structure ) {
				}

				if ( $nameSeq =~ /.*__(\d*)-(\d*)$/ ) {
					my $mfei = PipelineMiRNA::Utils::compute_mfei($dna, $energy);

					$hash[ $index++ ] = {
						"name"      => $nameSeq,
						"start"     => $1,
						"end"       => $2,
						"mfei"      => $mfei,
						"dna"       => $dna,
						"structure" => $structure,
						"energy"    => $energy
					};

				}
			}
			else {
				debug( "No structure found in $line2", 1 );
			}    # if $line2
		}
		else {

			# Should not happen
		}    #if $line1
	}    #while $line=<IN>
	return \@hash;
}

=method treat_candidates

Process the candidates and try merging them.

=cut

sub treat_candidates {

	my (@hash)   = @{ +shift };
	my %newHash  = ();
	my %tempHash = ();
	my $i        = 0;
	foreach my $key ( keys @hash ) {

		my $nb        = scalar @hash;
		my $start     = $hash[$key]{"start"};
		my $end       = $hash[$key]{"end"};
		my $mfei      = $hash[$key]{"mfei"};
		
		my $nameSeq   = $hash[$key]{"name"};
		my $structure = $hash[$key]{"structure"};
		my $dna       = $hash[$key]{"dna"};
		my $energy    = $hash[$key]{"energy"};
		if (
			(
				$end >= $hash[ $key + 1 ]{"end"}
				|| ( $hash[ $key + 1 ]{"start"} < ( $start + $end ) / 2 )
			)
			&& ( $key != $nb - 1 )
		  )
		{
			$tempHash{$nameSeq} = {
				"mfei"      => $mfei,
				"dna"       => $dna,
				"structure" => $structure,
				"energy"    => $energy
			};

		}
		else {
			$tempHash{$nameSeq} = {
				"mfei"      => $mfei,
				"dna"       => $dna,
				"structure" => $structure,
				"energy"    => $energy
			};
			my $max;
			my @keys =
			  sort { $tempHash{$a}{"mfei"} <=> $tempHash{$b}{"mfei"} }
			  keys(%tempHash);
			foreach my $key (@keys) {
				if ( $i == 0 ) {
					$max = $key;
					$newHash{$max}{'max'} = {
						"mfei"      => PipelineMiRNA::Utils::restrict_num_decimal_digits($tempHash{$key}{"mfei"},3),
						"dna"       => $tempHash{$key}{"dna"},
						"structure" => $tempHash{$key}{"structure"},
						"energy"    => $tempHash{$key}{"energy"}
					};

				}
				else {

					$newHash{$max}{$key} = {
						"mfei"      => PipelineMiRNA::Utils::restrict_num_decimal_digits($tempHash{$key}{"mfei"},3),
						"dna"       => $tempHash{$key}{"dna"},
						"structure" => $tempHash{$key}{"structure"},
						"energy"    => $tempHash{$key}{"energy"}
					};

				}
				$i++;
			}
			%tempHash = ();
			$i        = 0;
		}

	}

	return %newHash;
}

=method create_directories

Create the necessary directories.

=cut

sub create_directories {

	my (%newHash) = %{ +shift };
	my $current_sequence_dir = shift;
	
	my $candidate_counter = 0;
	
	foreach my $key ( sort keys %newHash ) {
	    $candidate_counter++;
		my $candidate_dir = File::Spec->catdir( $current_sequence_dir, $candidate_counter );
		mkdir $candidate_dir;

		#Writing seq.txt
		my $candidate_sequence =
		  File::Spec->catfile( $candidate_dir, 'seq.txt' );
		open( my $OUT, '>', $candidate_sequence )
		  or die "Error when opening $candidate_sequence: $!";
		print $OUT ">$key\n$newHash{$key}{'max'}{'dna'}\n";
		close $OUT;

		for my $name ( keys %{ $newHash{$key} } ) {
			if ( $name eq 'max' ) {
				process_outRNAFold(
					$candidate_dir,
					'optimal',
					$key,
					$newHash{$key}{'max'}{'dna'},
					$newHash{$key}{'max'}{'structure'},
					$newHash{$key}{'max'}{'energy'}
				);
				process_outRNAFold(
					$candidate_dir,
					'stemloop',
					$key,
					$newHash{$key}{'max'}{'dna'},
					$newHash{$key}{'max'}{'structure'},
					$newHash{$key}{'max'}{'energy'}
				);
			}
			else {

				#Writing alternativeCandidates.txt
				my $alternative_candidates =
				  File::Spec->catfile( $candidate_dir,
					'alternativeCandidates.txt' );
				open( my $OUT2, '>>', $alternative_candidates )
				  or die "Error when opening $alternative_candidates: $!";
				print $OUT2
">$name\t$newHash{$key}{$name}{'dna'}\t$newHash{$key}{$name}{'structure'}\t$newHash{$key}{$name}{'mfei'}\n";

			}
		}
		close $OUT;
	}

}

=method process_outRNAFold

Writing (pseudo) rnafold output

=cut

sub process_outRNAFold {
	my ( $candidate_dir, $suffix, $nameSeq, $dna, $structure, $energy ) = @_;

	my $candidate_rnafold_output =
	  File::Spec->catfile( $candidate_dir, "outRNAFold_$suffix.txt" );

	open( my $OUT2, '>', $candidate_rnafold_output )
	  or die "Error when opening $candidate_rnafold_output: $!";
	print $OUT2 ">$nameSeq\n$dna\n$structure ($energy)\n";
	close $OUT2;

}

=method process_tests

Perform the a posteriori tests for a given job

=cut

sub process_tests {
	my ( $job_dir ) = @_;
	debug( "A posteriori tests in $job_dir", 1 );
	##Traitement fichier de sortie outStemloop
	opendir DIR, $job_dir;    #ouverture répertoire job
	my @dirs;
	@dirs = readdir DIR;
	closedir DIR;
    my $candidates_dir = File::Spec->catdir($job_dir, 'candidates');
    mkdir $candidates_dir;
	foreach my $dir (@dirs)    # parcours du contenu
	{
		debug( "Considering $dir", 1 );
		my $sequence_dir = File::Spec->catdir( $job_dir, $dir );
		if (   $dir ne '.'
			&& $dir ne '..'
			&& -d $sequence_dir )    #si fichier est un répertoire
		{
			debug( "Entering sequence $sequence_dir", 1 );
			opendir DIR, $sequence_dir;    # ouverture du sous répertoire
			my @files;
			@files = readdir DIR;
			closedir DIR;
			foreach my $subDir (@files) {
				debug( "Considering $subDir", 1 );
				my $candidate_dir =
				  File::Spec->catdir( $sequence_dir, $subDir );
				if (   $subDir ne '.'
					&& $subDir ne '..'
					&& -d $candidate_dir
				  )    # si le fichier est de type repertoire
				{
					debug( "Entering candidate $subDir", 1 );
					process_tests_for_candidate( $candidate_dir, $subDir );
					debug( "Done with candidate $subDir", 1 );
					debug(
"Pseudo Serializing candidate information:\n $job_dir, $dir, $subDir",
						1
					);

					if (
						!eval {
							PipelineMiRNA::Candidate
							  ->serialize_candidate_information(
								$job_dir, $dir, $subDir, $candidates_dir );
						}
					  )
					{

						# Catching
						debug( "Serialization failed", 1 );
					}
					else {
						debug( "Done with serializing $subDir", 1 );

						# All is well
					}
				}    # foreach my $file (@files)
			}    # if directory
			debug( "Done with initial sequence $dir", 1 );
		}    # foreach my $dir (@dirs)
	}    #process_tests
	debug( "Done with all the tests", 1 );
	return 0;
}

=method process_tests_for_candidate

Perform the a posteriori tests for a given candidate

=cut

sub process_tests_for_candidate {

	my @args = @_;
	my ( $candidate_dir, $file ) = @args;

	####Traitement fichier de sortie outStemloop
	chmod 0777, $candidate_dir;

	my $seq_file = File::Spec->catfile( $candidate_dir, 'seq.txt' );
	my $candidate_rnafold_optimal_out =
	  File::Spec->catfile( $candidate_dir, 'outRNAFold_optimal.txt' );
	my $candidate_rnafold_stemploop_out =
	  File::Spec->catfile( $candidate_dir, 'outRNAFold_stemloop.txt' );

	####conversion en format CT
	my $candidate_ct_optimal_file =
	  File::Spec->catfile( $candidate_dir, 'outB2ct_optimal.ct' );
	debug( "Converting optimal to CT in $candidate_ct_optimal_file", 1 );
	PipelineMiRNA::Programs::convert_to_ct( $candidate_rnafold_optimal_out,
		$candidate_ct_optimal_file )
	  or die('Problem when converting to CT format');

	my $candidate_ct_stemloop_file =
	  File::Spec->catfile( $candidate_dir, 'outB2ct_stemloop.ct' );
	debug( "Converting stemloop to CT in $candidate_ct_stemloop_file", 1 );
	PipelineMiRNA::Programs::convert_to_ct( $candidate_rnafold_stemploop_out,
		$candidate_ct_stemloop_file )
	  or die('Problem when converting to CT format');

	my $varna_image = File::Spec->catfile( $candidate_dir, 'image.png' );
	debug( "Generating image using VARNA in $varna_image", 1 );
	PipelineMiRNA::Programs::run_varna( $candidate_ct_stemloop_file,
		$varna_image )
	  or die('Problem during image generation using VARNA');

	my $cfg = PipelineMiRNA->CONFIG();

	####calcul MFEI (appel script energie.pl)
	if ( $cfg->param('options.mfe') ) {
		debug( "Running test_mfei on $file", 1 );
		PipelineMiRNA::PosterioriTests::test_mfei( $candidate_dir,
			$candidate_ct_optimal_file, $file );
	}
	####calcul p-value randfold
	if ( $cfg->param('options.randfold') ) {
		debug( "Running test_randfold on $seq_file", 1 );
		PipelineMiRNA::PosterioriTests::test_randfold( $candidate_dir,
			$seq_file );
	}
	if ( $cfg->param('options.align') ) {
		debug( "Running test_alignment on $candidate_ct_stemloop_file", 1 );
		my $file_alignement = PipelineMiRNA::PosterioriTests::test_alignment( $candidate_dir,
			$candidate_ct_stemloop_file );
		post_process_alignments($candidate_dir, $candidate_rnafold_stemploop_out, $file_alignement);
	}    # if file

	return;
}

=method post_process_alignments


=cut

sub post_process_alignments {
    my @args = @_;
    my $candidate_dir = shift @args;
    my $candidate_rnafold_stemploop_out = shift @args;
    my $file_alignement = shift @args;

    my @res = PipelineMiRNA::Components::get_data_from_rnafold_out($candidate_rnafold_stemploop_out);
    my ($name, $position, $DNASequence, $Vienna) = @res;
    my %alignments;
    if (! eval {%alignments = PipelineMiRNA::Components::parse_custom_exonerate_output($file_alignement);}) {
        # Catching exception
    } else {
        %alignments = PipelineMiRNA::Components::merge_alignments(\%alignments);
        my $tmp_file = File::Spec->catfile($candidate_dir, "mirdup_validation.txt");
        my %mirdup_results = PipelineMiRNA::MiRdup->validate_with_mirdup($tmp_file, $name,
                                                                         $DNASequence, $Vienna,
                                                                         keys %alignments);
        my $mirdup_results_file = File::Spec->catfile($candidate_dir, 'mirdup_results.yml');
        YAML::XS::DumpFile($mirdup_results_file, %mirdup_results);

        my $alignments_results_file = File::Spec->catfile($candidate_dir, 'merged_alignments.yml');
        YAML::XS::DumpFile($alignments_results_file, %alignments);
    }
}


1;
