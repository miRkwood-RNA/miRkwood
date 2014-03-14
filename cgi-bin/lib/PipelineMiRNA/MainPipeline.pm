package PipelineMiRNA::MainPipeline;

# ABSTRACT: The actual pipeline

use strict;
use warnings;

use File::Path 'rmtree';
use File::Basename;
use Cwd qw( abs_path );
use Carp;
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
	my ( $strands, $mfe, $randfold, $align, $run_options_file ) = @_;
	my $run_options = PipelineMiRNA->CONFIG();
	$run_options->param( "options.strands",  $strands );
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
	my ( $filter, $strand, $mfe, $randfold, $align, $job_dir, $plant ) = @_;
	my $log_file = File::Spec->catfile( $job_dir, 'log.log' );
	local $Log::Message::Simple::DEBUG_FH = PipelineMiRNA->LOGFH($log_file);
    PipelineMiRNA->DEBUG(1);
	my $run_options_file = PipelineMiRNA::Paths->get_job_config_path($job_dir);
	PipelineMiRNA->CONFIG_FILE($run_options_file);
	write_config( $strand, $mfe, $randfold, $align, $run_options_file );

	debug( 'BEGIN execute_scripts', PipelineMiRNA->DEBUG() );
	my $sequences_input = File::Spec->catfile( $job_dir, 'Sequences.fas' );
	if ($filter) {
		debug( 'FilteringCDS', PipelineMiRNA->DEBUG() );
		PipelineMiRNA::Components::filter_CDS( $dirData, $job_dir, $plant );
	}
	else {
		my $sequence_uploaded =
		  File::Spec->catfile( $job_dir, 'sequenceUpload.fas' );
		debug( "Moving file $sequence_uploaded to $sequences_input", PipelineMiRNA->DEBUG() );
		File::Copy::move( $sequence_uploaded, $sequences_input );
	}

	##Passage du multifasta -> fasta et appel au script Stemloop
	debug( "Opening multifasta $sequences_input", PipelineMiRNA->DEBUG() );
	open my $ENTREE_FH, '<', $sequences_input
	  or die "Error when opening sequences -$sequences_input-: $!";
	debug( "Calling parse_multi_fasta() on $sequences_input", PipelineMiRNA->DEBUG() );
	my %tab = PipelineMiRNA::Utils::parse_multi_fasta($ENTREE_FH);
	close $ENTREE_FH;

	debug( 'Iterating over names', PipelineMiRNA->DEBUG() );

	my $sequence_dir_name = 0;

	while ( my ( $name, $sequence ) = each %tab ) {
		debug( "Considering sequence $sequence_dir_name: $name", PipelineMiRNA->DEBUG() );
		$sequence_dir_name++;
		my $sequence_dir = File::Spec->catdir( $job_dir, $sequence_dir_name );
		mkdir $sequence_dir;

		my $res = process_sequence( $sequence_dir, $name, $sequence, '+' );
		my @hash1 = @{$res};
		my @hash;
		my $cfg = PipelineMiRNA->CONFIG();
		if ( $cfg->param('options.strands') ) {
			debug( "Processing the other strand", PipelineMiRNA->DEBUG() );
			my $reversed_sequence =
			  PipelineMiRNA::Utils::reverse_complement($sequence);
			my $res2 =
			  process_sequence( $sequence_dir, $name, $reversed_sequence, '-' );
			my @hash2 = @{$res2};
			@hash = sort { $a->{start} <=> $b->{start} } ( @hash1, @hash2 );
		}
		else {
			@hash = @hash1;
		}

        if ( $cfg->param('options.mfe') ) {
            debug('Select only sequences with MFEI < -0.6', PipelineMiRNA->DEBUG() );
            @hash = grep { mfei_below_threshold($_, -0.6) } @hash;
        }

		my %candidates_hash;
		if (@hash) {
			%candidates_hash = merge_candidates( \@hash );
		}
		else {
		    %candidates_hash = ();
		}

		create_directories( \%candidates_hash, $sequence_dir );
	}
	process_tests($job_dir);
	return;
}

=method mfei_below_threshold

Return whether a given candidate has its mfei above a given threshold

=cut

sub mfei_below_threshold {
    my @args      = @_;
    my $candidate = shift @args;
    my $threshold = shift @args;
    my %current_candidate = %{ $candidate };
    my $mfei = $current_candidate{'mfei'};
    return $mfei < $threshold;
}

=method process_sequence

Process a single sequence

 Usage : process_sequence( $sequence_dir, $name, $sequence );
 Return: -

=cut

sub process_sequence {
	my @args         = @_;
	my $sequence_dir = shift @args;
	my $name         = shift @args;
	my $sequence     = shift @args;
	my $strand       = shift @args;

	## Running RNALfold
	debug( 'Running RNALfold', PipelineMiRNA->DEBUG() );
	my $rnalfold_output = File::Spec->catfile( $sequence_dir, 'RNALfold.out' );

	my $temp_file = File::Spec->catfile( $sequence_dir, 'tempFile.txt' );
	PipelineMiRNA::Programs::run_rnalfold( $name, $sequence, $temp_file,
		$rnalfold_output )
	  or die("Problem when running RNALfold: $!");

	## Running RNAstemloop
	my $rnastemloop_out_optimal =
	  File::Spec->catfile( $sequence_dir, 'rnastemloop_optimal.out' );
	my $rnastemloop_out_stemloop =
	  File::Spec->catfile( $sequence_dir, 'rnastemloop_stemloop.out' );
	debug( "Running RNAstemloop on $rnalfold_output", PipelineMiRNA->DEBUG() );
	PipelineMiRNA::Programs::run_rnastemloop( $rnalfold_output,
		$rnastemloop_out_stemloop, $rnastemloop_out_optimal )
	  or die("Problem when running RNAstemloop");

	my $rnaeval_out_optimal =
	  run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_optimal, 'optimal' );
    my $rnaeval_out_stemloop =
	run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_stemloop,  'stemloop' );
	my $seq_length = length $sequence;
	open( my $STEM_FH, '<', $rnastemloop_out_stemloop ) or die $!;
	open( my $EVAL_OPT_FH, '<', $rnaeval_out_optimal ) or die $!;
	open( my $EVAL_STEM_FH, '<', $rnaeval_out_stemloop ) or die $!;
	my $res =
	  process_RNAstemloop( $sequence_dir, $strand, $seq_length, $STEM_FH,
		$EVAL_OPT_FH, $EVAL_STEM_FH );
	close($STEM_FH);
	close($EVAL_OPT_FH);
	close($EVAL_STEM_FH);
	return $res;
}

=method run_RNAeval_on_RNAstemloop_output


=cut

sub run_RNAeval_on_RNAstemloop_output {
	my ( $rnastemloop_out, $suffix ) = @_;
	my $current_sequence_dir = dirname($rnastemloop_out);
	debug( "Processing RNAstemloop output for $suffix $rnastemloop_out", PipelineMiRNA->DEBUG() );
	my $rnaeval_out =
	  File::Spec->catfile( $current_sequence_dir, "rnaeval_$suffix.out" );

	debug( "Running RNAeval in $rnaeval_out", PipelineMiRNA->DEBUG() );
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
	my ($strand)       = shift @args;
	my ($seq_length)   = shift @args;
	my ($STEM_FH)      = shift @args;
	my ($EVAL_OPT_FH)  = shift @args;
	my ($EVAL_STEM_FH) = shift @args;
	my $index          = 0;
	my ($line_eval_opt, $line_eval_stem);
	my ( $nameSeq, $dna, $structure_stemloop );
	my @hash = ();

	while ( my $stem_line = <$STEM_FH> ) {

		if ( ( $stem_line =~ /^>(.*)/ ) ) {    # nom sequence
			$nameSeq = $1;
		}
		elsif ( ( $stem_line =~ /^[a-zA-Z]/ ) ) { # récupération de la sequence adn
			$dna = substr $stem_line, 0, -1;
			$line_eval_opt = substr( <$EVAL_OPT_FH>, 0, -1 );    # the sequence as well
			$line_eval_stem = substr( <$EVAL_STEM_FH>, 0, -1 );    # the sequence as well
			if ( $dna ne $line_eval_opt || $dna ne $line_eval_stem ) {
				warn ('The sequences differ in RNAeval and RNAstemloop output');
			}
		}
		elsif ( ( $stem_line =~ /(.*)/ ) ) {
			$structure_stemloop = $1;
			$line_eval_opt = <$EVAL_OPT_FH>;    # the structure as well, and the energy
			$line_eval_stem = <$EVAL_STEM_FH>;
			
			my ( $structure_optimal, $energy_optimal ) =
                PipelineMiRNA::Parsers::parse_Vienna_line($line_eval_opt);
			my ( $structure_stemloop, $energy_stemloop ) =
                PipelineMiRNA::Parsers::parse_Vienna_line($line_eval_stem);
			if ($structure_optimal)
			{                        # We have a structure

				if ( $nameSeq =~ /.*__(\d*)-(\d*)$/ ) {
					my $mfei =
					  PipelineMiRNA::Utils::compute_mfei( $dna, $energy_optimal );
					my ( $start, $end );
					if ( $strand eq '-' ) {
						my $res =
						  PipelineMiRNA::Utils::get_position_from_opposite_strand
						  ( $1, $2, $seq_length );
						( $start, $end ) = @{$res};
					}
					else {
						( $start, $end ) = ( $1, $2 );
					}
					$hash[ $index++ ] = {
						"name"      => $nameSeq,
						"start"     => $start,
						"end"       => $end,
						"mfei"      => $mfei,
						"dna"       => $dna,
						"structure_optimal" => $structure_optimal,
						"structure_stemloop" => $structure_stemloop,
						"strand"    => $strand,
						"energy_optimal" => $energy_optimal,
						"energy_stemloop" => $energy_stemloop
					};

				}
			}
			else {
				warn( "No structure found in $line_eval_opt" );
			}    # if $line2
		}
		else {
            warn( "Unrecognized line " );
		}    #if $line1
	}    #while $line=<IN>
	return \@hash;
}

=method is_overlapping

Test whether one candidate is overlapping with the other
(based on their positions).

=cut

sub is_overlapping {
	my @args      = @_;
	my $start     = shift @args or die('Not enough values provided');
	my $end       = shift @args or die('Not enough values provided');
	my $ref_start = shift @args or die('Not enough values provided');
	my $ref_end   = shift @args or die('Not enough values provided');
	$ref_start <= $start or die('Positions should be ordered');
	return ( $start < ( $ref_start + $ref_end ) / 2 );
}

=method is_included

Test whether one candidate is included into the other
(based on their positions).

=cut

sub is_included {
	my @args      = @_;
	my $start     = shift @args or die('Not enough values provided');
	my $end       = shift @args or die('Not enough values provided');
	my $ref_start = shift @args or die('Not enough values provided');
	my $ref_end   = shift @args or die('Not enough values provided');
	$ref_start <= $start or die('Positions should be ordered');
	return ( $end <= $ref_end );
}

=method merge_candidates

Process the candidates and try merging them.
We assume the candidates array is already sorted by growing position

=cut

sub merge_candidates {
	my (@candidates_array) = @{ +shift };
	my $nb_candidates = scalar @candidates_array;

	my @merged_candidates   = ();
	my %final_hash          = ();
	my %reference_candidate = %{ $candidates_array[0] };
	my %best_candidate      = %reference_candidate;
	my %current_candidate;
	for my $candidate_index ( 1 .. $#candidates_array ) {
		%current_candidate = %{ $candidates_array[$candidate_index] };
		my $start = $current_candidate{'start'};
		my $end   = $current_candidate{'end'};
		my ( $ref_start, $ref_end ) =
		  ( $reference_candidate{'start'}, $reference_candidate{'end'} );
		if ( is_included( $start, $end, $ref_start, $ref_end ) ) {
			if ( $best_candidate{'mfei'} <= -0.8 ) {
				push @merged_candidates, {%current_candidate};
			}
			else {
				if ( $current_candidate{'mfei'} < $best_candidate{'mfei'} ) {
					push @merged_candidates, {%best_candidate};
					%best_candidate = %current_candidate;
				}
				else {
					push @merged_candidates, {%current_candidate};
				}
			}

		}
		elsif ( is_overlapping( $start, $end, $ref_start, $ref_end ) ) {
			if ( $current_candidate{'mfei'} < $best_candidate{'mfei'} ) {
				push @merged_candidates, {%best_candidate};
				%best_candidate = %current_candidate;
			}
			else {
				push @merged_candidates, {%current_candidate};
			}
		}
		else {
			my $final_name = $best_candidate{'name'};
			$final_hash{$final_name}                 = {};
			$final_hash{$final_name}{'max'}          = {%best_candidate};
			$final_hash{$final_name}{'alternatives'} = [@merged_candidates];
			@merged_candidates                       = undef;
			@merged_candidates                       = ();
			%reference_candidate                     = %current_candidate;
			%best_candidate                          = %reference_candidate;
		}

	}    #foreach
	my $final_name = $best_candidate{'name'};
	$final_hash{$final_name}                 = {};
	$final_hash{$final_name}{'max'}          = {%best_candidate};
	$final_hash{$final_name}{'alternatives'} = [@merged_candidates];
	return %final_hash;
}

=method create_directories

Create the necessary directories.

=cut

sub create_directories {
	my @args                 = @_;
	my (%candidates_hash)    = %{ shift @args };
	my $current_sequence_dir = shift @args;
	my $candidate_counter    = 0;
	foreach my $key ( sort keys %candidates_hash ) {
		$candidate_counter++;
		my $candidate_dir =
		  File::Spec->catdir( $current_sequence_dir, $candidate_counter );
		mkdir $candidate_dir;
		my $candidate_ref = $candidates_hash{$key}{'max'};
		populate_candidate_directory( $candidate_dir, $candidate_ref,
			$candidates_hash{$key}{'alternatives'} );
	}
	return;
}

=method populate_candidate_directory

Populate a candidate directory with the sequence, strand & so on.

=cut

sub populate_candidate_directory {
	my @args               = @_;
	my $candidate_dir      = shift @args;
	my %candidate          = %{ shift @args };
	my @alternatives_array = @{ shift @args };

	#Writing seq.txt
	my $candidate_sequence = File::Spec->catfile( $candidate_dir, 'seq.txt' );
	open( my $SEQ_FH, '>', $candidate_sequence )
	  or die "Error when opening $candidate_sequence: $!";
	print $SEQ_FH ">$candidate{'name'}\n$candidate{'dna'}\n";
	close $SEQ_FH;

	#Writing strand
	my $strand_file = File::Spec->catfile( $candidate_dir, 'strand.txt' );
	open( my $STRAND_FH, '>', $strand_file )
	  or die "Error when opening $strand_file: $!";
	print $STRAND_FH $candidate{'strand'};
	close $STRAND_FH;

	process_outRNAFold( $candidate_dir, 'optimal', $candidate{'name'},
		$candidate{'dna'}, $candidate{'structure_optimal'}, $candidate{'energy_optimal'} );
	process_outRNAFold( $candidate_dir, 'stemloop', $candidate{'name'},
		$candidate{'dna'}, $candidate{'structure_stemloop'}, $candidate{'energy_stemloop'} );

	#Writing alternativeCandidates.txt
	my $alternative_candidates =
	  File::Spec->catfile( $candidate_dir, 'alternativeCandidates.txt' );
	if (@alternatives_array){
		open( my $OUT2, '>>', $alternative_candidates )
		  or die "Error when opening $alternative_candidates: $!";
		foreach my $alternative (@alternatives_array) {
			print $OUT2
">$alternative->{'name'}\t$alternative->{'dna'}\t$alternative->{'structure_optimal'}\t$alternative->{'mfei'}\n";
		}
		close $OUT2;
    }
	return;
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
	my ($job_dir) = @_;
	debug( "A posteriori tests in $job_dir", PipelineMiRNA->DEBUG() );
	##Traitement fichier de sortie outStemloop
	opendir DIR, $job_dir;    #ouverture répertoire job
	my @dirs;
	@dirs = readdir DIR;
	closedir DIR;
	my $candidates_dir = File::Spec->catdir( $job_dir, 'candidates' );
	mkdir $candidates_dir;
	foreach my $dir (@dirs)    # parcours du contenu
	{
		my $sequence_dir = File::Spec->catdir( $job_dir, $dir );
		if (   $dir ne '.'
			&& $dir ne '..'
			&& -d $sequence_dir )    #si fichier est un répertoire
		{
			debug( "Entering sequence $sequence_dir", PipelineMiRNA->DEBUG() );
			opendir DIR, $sequence_dir;    # ouverture du sous répertoire
			my @files;
			@files = readdir DIR;
			closedir DIR;
			foreach my $subDir (@files) {
				my $candidate_dir =
				  File::Spec->catdir( $sequence_dir, $subDir );
				if (   $subDir ne '.'
					&& $subDir ne '..'
					&& -d $candidate_dir
				  )    # si le fichier est de type repertoire
				{
					debug( "Entering candidate $subDir", PipelineMiRNA->DEBUG() );
					process_tests_for_candidate( $candidate_dir, $subDir );
					debug( "Done with candidate $subDir", PipelineMiRNA->DEBUG() );

					if (
						!eval {
							PipelineMiRNA::Candidate
							  ->serialize_candidate_information( $job_dir, $dir,
								$subDir, $candidates_dir );
						}
					  )
					{

						# Catching
						carp( "Serialization failed" );
					}
					else {
						debug( "Done with serializing $subDir", PipelineMiRNA->DEBUG() );

						# All is well
					}
				}    # foreach my $file (@files)
			}    # if directory
			debug( "Done with initial sequence $dir", PipelineMiRNA->DEBUG() );
		}    # foreach my $dir (@dirs)
	}    #process_tests
	debug( "Done with all the tests", PipelineMiRNA->DEBUG() );
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
	debug( "Converting optimal to CT in $candidate_ct_optimal_file", PipelineMiRNA->DEBUG() );
	PipelineMiRNA::Programs::convert_to_ct( $candidate_rnafold_optimal_out,
		$candidate_ct_optimal_file )
	  or die('Problem when converting to CT format');

	my $candidate_ct_stemloop_file =
	  File::Spec->catfile( $candidate_dir, 'outB2ct_stemloop.ct' );
	debug( "Converting stemloop to CT in $candidate_ct_stemloop_file", PipelineMiRNA->DEBUG() );
	PipelineMiRNA::Programs::convert_to_ct( $candidate_rnafold_stemploop_out,
		$candidate_ct_stemloop_file )
	  or die('Problem when converting to CT format');

	my $varna_image = File::Spec->catfile( $candidate_dir, 'image.png' );
	debug( "Generating image using VARNA in $varna_image", PipelineMiRNA->DEBUG() );
	PipelineMiRNA::Programs::run_varna( $candidate_ct_stemloop_file,
		$varna_image )
	  or carp('Problem during image generation using VARNA');

	my $cfg = PipelineMiRNA->CONFIG();

	####calcul MFEI (appel script energie.pl)

	debug( "Running test_mfei on $file", PipelineMiRNA->DEBUG() );
	PipelineMiRNA::PosterioriTests::test_mfei( $candidate_dir,
		$candidate_ct_optimal_file, $file );

	####calcul p-value randfold
	if ( $cfg->param('options.randfold') ) {
		debug( "Running test_randfold on $seq_file", PipelineMiRNA->DEBUG() );
		PipelineMiRNA::PosterioriTests::test_randfold( $candidate_dir,
			$seq_file );
	}
	if ( $cfg->param('options.align') ) {
		debug( "Running test_alignment on $candidate_ct_stemloop_file", PipelineMiRNA->DEBUG() );
		my $file_alignement =
		  PipelineMiRNA::PosterioriTests::test_alignment( $candidate_dir,
			$candidate_ct_stemloop_file );
		post_process_alignments( $candidate_dir,
			$candidate_rnafold_stemploop_out,
			$file_alignement );
	}    # if file

	return;
}

=method post_process_alignments


=cut

sub post_process_alignments {
	my @args                            = @_;
	my $candidate_dir                   = shift @args;
	my $candidate_rnafold_stemploop_out = shift @args;
	my $file_alignement                 = shift @args;

	my @res =
	  PipelineMiRNA::Components::get_data_from_rnafold_out(
		$candidate_rnafold_stemploop_out);
	my ( $name, $position, $DNASequence, $Vienna ) = @res;
	my %alignments;
	if (
		!eval {
			%alignments =
			  PipelineMiRNA::Components::parse_custom_exonerate_output(
				$file_alignement);
		}
	  )
	{

		# Catching exception
	}
	else {
		%alignments =
		  PipelineMiRNA::Components::merge_alignments( \%alignments );
		my $tmp_file =
		  File::Spec->catfile( $candidate_dir, "mirdup_validation.txt" );
		my %mirdup_results =
		  PipelineMiRNA::MiRdup->validate_with_mirdup( $tmp_file, $name,
			$DNASequence, $Vienna, keys %alignments );
		my $mirdup_results_file =
		  File::Spec->catfile( $candidate_dir, 'mirdup_results.yml' );
		YAML::XS::DumpFile( $mirdup_results_file, %mirdup_results );

		my $alignments_results_file =
		  File::Spec->catfile( $candidate_dir, 'merged_alignments.yml' );
		YAML::XS::DumpFile( $alignments_results_file, %alignments );
	}
}

1;
