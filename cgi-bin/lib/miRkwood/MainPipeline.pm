package miRkwood::MainPipeline;

# ABSTRACT: The actual pipeline

use strict;
use warnings;

use File::Path 'rmtree';
use File::Basename;
use Cwd qw( abs_path );
use Carp;
use File::Copy;
use miRkwood;
use miRkwood::Paths;
use miRkwood::Utils;
use miRkwood::MiRdup;
use miRkwood::Parsers;
use miRkwood::Programs;
use miRkwood::CandidateHandler;
use miRkwood::Components;
use miRkwood::Maskers;
use miRkwood::PosterioriTests;
use miRkwood::FileUtils;
use Log::Message::Simple qw[msg error debug];

### Data ##
my $dirData = miRkwood::Paths->get_data_path();

=method setup_logging

 Usage : setup_logging($job_dir)
 Input : The job directory
 Return: -

=cut

sub setup_logging {
    my @args = @_;
    my $job_dir = shift @args;
    my $log_file = File::Spec->catfile( $job_dir, 'log.log' );
    $Log::Message::Simple::DEBUG_FH = miRkwood->LOGFH($log_file);
    miRkwood->DEBUG(1);
    return;
}

=method init_pipeline

Initialise the pipeline setup

 Usage : init_pipeline($job_dir)
 Input : The job directory
 Return: -

=cut

sub init_pipeline {
    my @args = @_;
    my $job_dir = shift @args;
    setup_logging($job_dir);
    my $run_options_file = miRkwood::Paths->get_job_config_path($job_dir);
    miRkwood->CONFIG_FILE($run_options_file);
    miRkwood::Programs::init_programs();
    return;
}

=method fasta_pipeline

Run the pipeline in FASTA mode

 Usage : miRkwood::MainPipeline::fasta_pipeline( $idirJob );
 Input : The job directory
 Return: -

=cut

sub fasta_pipeline {
    my ( $job_dir  ) = @_;
    init_pipeline($job_dir);
    my @sequences_array = get_sequences($job_dir);
    run_pipeline_on_sequences($job_dir, \@sequences_array);
    return;
}

=method bam_pipeline

Run the pipeline in BAM mode

 Usage : miRkwood::MainPipeline::bam_pipeline( $idirJob, $bam_file, $genome_file );
 Input : The job directory
 Return: -

=cut

sub bam_pipeline {
    my @args = @_;
    my ( $job_dir, $bam_file, $genome_file ) = @args;
    init_pipeline($job_dir);
    my $clustering = miRkwood::Clusters->new($bam_file, $genome_file);
    debug( "Extracting sequences from genome using BAM clusters", miRkwood->DEBUG() );
    my @sequences_array = $clustering->get_clustered_sequences_from_bam();
    run_pipeline_on_sequences($job_dir, \@sequences_array);
    return;
}

=method run_pipeline_on_sequences

Run the pipeline on the given sequences

=cut

sub run_pipeline_on_sequences {
    my @args = @_;
    my $job_dir = shift @args;
    my $sequences_array = shift @args;
    my @sequences_array = @{$sequences_array};
    my $sequences_count = scalar @sequences_array;
    debug( "$sequences_count sequences to process", miRkwood->DEBUG() );
    my $workspace_dir = miRkwood::Paths->get_workspace_path($job_dir);
    mkdir $workspace_dir;
    my $sequence_dir_name = 0;
    foreach my $item (@sequences_array){
        my ($name, $sequence) = @{$item};
        debug( "Considering sequence $sequence_dir_name: $name", miRkwood->DEBUG() );
        $sequence_dir_name++;
        my $sequence_dir = File::Spec->catdir( $workspace_dir, $sequence_dir_name );
        mkdir $sequence_dir;
        my %candidates_hash = compute_candidates_for_sequence($sequence_dir, $item);
        create_directories( \%candidates_hash, $sequence_dir );
    }
    process_tests($job_dir);
    debug('miRkwood processing done', miRkwood->DEBUG() );
    mark_job_as_finished($job_dir);
    debug("Writing finish file", miRkwood->DEBUG() );
    return;
}

=method compute_candidates_for_sequence

Computes the candidates of the given sequence (as a couple <name, sequence>)
using the given workspace.

 Usage : my %candidates_hash = compute_candidates_for_sequence($sequence_dir, $seq_record);
 Return: A hash of candidates

=cut

sub compute_candidates_for_sequence {
    my @args = @_;
    my $sequence_dir = shift @args;
    my $sequence_tuple = shift @args;
    my ($name, $sequence) = @{$sequence_tuple};
    my $cfg = miRkwood->CONFIG();
    my $candidates = get_raw_candidates_from_sequence( $sequence_dir, $name, $sequence, '+' );
    my @candidates_hash1 = @{$candidates};
    my @candidates_hash;

    if ( $cfg->param('options.strands') ) {
        debug( "Processing the other strand", miRkwood->DEBUG() );
        my $reversed_sequence =
          miRkwood::Utils::reverse_complement($sequence);
        my $candidates2 =
          get_raw_candidates_from_sequence( $sequence_dir, $name, $reversed_sequence, '-' );
        my @candidates_hash2 = @{$candidates2};
        @candidates_hash = sort { $a->{start} <=> $b->{start} } ( @candidates_hash1, @candidates_hash2 );
    }
    else {
        @candidates_hash = sort { $a->{start} <=> $b->{start} } ( @candidates_hash1 );
    }

    if ( $cfg->param('options.mfe') ) {
        debug('Select only sequences with MFEI < -0.6', miRkwood->DEBUG() );
        @candidates_hash = grep { mfei_below_threshold($_, -0.6) } @candidates_hash;
    }

    my %candidates_hash;
    if (@candidates_hash) {
        debug("Merging candidates for sequence $name", miRkwood->DEBUG() );
        %candidates_hash = merge_candidates( \@candidates_hash );
    }
    else {
        %candidates_hash = ();
    }
    return %candidates_hash;
}


=method get_sequences

Get the sequences to process from the job directory
This includes parsing, and masking if option selected.

 Usage : my @fasta_array = get_sequences($job_dir);
 Input : The job directory
 Return: An array of couples (name, sequence)

=cut

sub get_sequences {
    my @args = @_;
    my $job_dir = shift @args;

    my %filter = get_masking_information($job_dir);
    my $sequence_uploaded =
          File::Spec->catfile( $job_dir, 'input_sequences.fas' );

    open my $ENTREE_FH, '<', $sequence_uploaded
      or die "Error when opening sequences -$sequence_uploaded-: $!";
    debug( "Calling parse_multi_fasta() on $sequence_uploaded", miRkwood->DEBUG() );
    my @sequences_array = miRkwood::Utils::parse_multi_fasta($ENTREE_FH);
    close $ENTREE_FH;
    my @sequences;
    if (%filter){
        debug( 'Masking input sequences', miRkwood->DEBUG() );
        @sequences = miRkwood::Maskers::mask_sequences(\%filter, @sequences_array);
    } else {
        @sequences = @sequences_array;
    }
    return @sequences;
}

=method get_masking_information

Get the masking information based on the job configuration

 Usage : my %filter = get_masking_information($job_dir);
 Input : The job directory
 Return: A hash (name => [ positions ])

=cut

sub get_masking_information {
    my @args = @_;
    my $job_dir = shift @args;
    my %filter;
    my $cfg = miRkwood->CONFIG();
    my $masking_folder = File::Spec->catdir($job_dir, 'masks');
    mkdir $masking_folder;
    my $sequences = File::Spec->catfile( $job_dir, 'input_sequences.fas' );

    if ($cfg->param('options.filter')) {
        debug( 'Get masking information for coding regions', miRkwood->DEBUG() );
        my $plant = $cfg->param('job.plant');
        my $blast_database = File::Spec->catfile( $dirData, "$plant.fas" );
        my %blast_mask = miRkwood::Maskers::get_coding_region_masking_information( $sequences, $masking_folder, $blast_database );
        %filter = miRkwood::Utils::merge_hashes_of_arrays(\%filter, \%blast_mask);
    }
    if ($cfg->param('options.mask-trna')) {
        debug( 'Get masking information for tRNAs', miRkwood->DEBUG() );
        my %trna_mask = miRkwood::Maskers::get_trna_masking_information( $sequences, $masking_folder );
        %filter = miRkwood::Utils::merge_hashes_of_arrays(\%filter, \%trna_mask);
    }
    if ($cfg->param('options.mask-rrna')) {
        debug( 'Get masking information for ribosomal RNAs', miRkwood->DEBUG() );
        my %rrna_mask = miRkwood::Maskers::get_rnammer_masking_information( $sequences, $masking_folder );
        %filter = miRkwood::Utils::merge_hashes_of_arrays(\%filter, \%rrna_mask);
    }
    return %filter;
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

=method get_raw_candidates_from_sequence

Process a single sequence record (as a couple <name, sequence>)
and return the raw candidates per RNALfold/RNAstemloop.

 Usage : process_sequence( $sequence_dir, $name, $sequence );
 Return: a list of candidates (as hashes)

=cut

sub get_raw_candidates_from_sequence {
	my @args         = @_;
	my $sequence_dir = shift @args;
	my $name         = shift @args;
	my $sequence     = shift @args;
	my $strand       = shift @args;

	my $rnalfold_output = run_rnalfold_on_sequence($name, $sequence, $sequence_dir);

	my ($rnastemloop_out_stemloop, $rnastemloop_out_optimal) =
	   run_RNAstemloop_on_rnalfold_output($rnalfold_output, $sequence_dir);

	my $rnaeval_out_optimal =
	  run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_optimal, 'optimal' );
    my $rnaeval_out_stemloop =
	run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_stemloop,  'stemloop' );
	my $seq_length = length $sequence;
	my $res = process_RNAstemloop_on_filenames( $sequence_dir, $strand, $seq_length,
	                                            $rnastemloop_out_stemloop,
	                                            $rnaeval_out_optimal,
	                                            $rnaeval_out_stemloop );
	return $res;
}

=method run_rnalfold_on_sequence

 Usage : run_rnalfold_on_sequence( $name, $sequence, $sequence_dir);
 Return: $rnalfold_output

=cut

sub run_rnalfold_on_sequence {
    my @args = @_;
    my $name     = shift @args;
    my $sequence = shift @args;
    my $sequence_dir = shift @args;
    debug( 'Running RNALfold', miRkwood->DEBUG() );
    my $rnalfold_output = File::Spec->catfile( $sequence_dir, 'RNALfold.out' );

    my $temp_file = File::Spec->catfile( $sequence_dir, 'tempFile.txt' );
    miRkwood::Programs::run_rnalfold( $name, $sequence, $temp_file,
        $rnalfold_output )
      or die("Problem when running RNALfold: $!");
    return $rnalfold_output;
}


=method run_RNAstemloop_on_rnalfold_output

 Usage : run_RNAstemloop_on_rnalfold_output( $rnalfold_output, $sequence_dir );
 Return: ($rnastemloop_out_stemloop, $rnastemloop_out_optimal);

=cut

sub run_RNAstemloop_on_rnalfold_output {
    my @args = @_;
    my $rnalfold_output = shift @args;
    my $sequence_dir    = shift @args;
    my $rnastemloop_out_optimal = File::Spec->catfile( $sequence_dir, 'rnastemloop_optimal.out' );
    my $rnastemloop_out_stemloop =
      File::Spec->catfile( $sequence_dir, 'rnastemloop_stemloop.out' );
    debug( "Running RNAstemloop on $rnalfold_output", miRkwood->DEBUG() );
    miRkwood::Programs::run_rnastemloop( $rnalfold_output,
        $rnastemloop_out_stemloop, $rnastemloop_out_optimal )
      or die("Problem when running RNAstemloop");
    return ($rnastemloop_out_stemloop, $rnastemloop_out_optimal);
}

=method run_RNAeval_on_RNAstemloop_output


=cut

sub run_RNAeval_on_RNAstemloop_output {
	my ( $rnastemloop_out, $suffix ) = @_;
	my $current_sequence_dir = dirname($rnastemloop_out);
	debug( "Processing RNAstemloop output for $suffix $rnastemloop_out", miRkwood->DEBUG() );
	my $rnaeval_out =
	  File::Spec->catfile( $current_sequence_dir, "rnaeval_$suffix.out" );

	debug( "Running RNAeval in $rnaeval_out", miRkwood->DEBUG() );
	miRkwood::Programs::run_rnaeval( $rnastemloop_out, $rnaeval_out )
	  or die("Problem when running RNAeval");

	return $rnaeval_out;
}

=method process_RNAstemloop_on_filenames

Pass-through method for process_RNAstemloop

=cut

sub process_RNAstemloop_on_filenames {
    my @args           = @_;
    my ($sequence_dir) = shift @args;
    my ($strand)       = shift @args;
    my ($seq_length)   = shift @args;
    my ($rnastemloop_out_stemloop)      = shift @args;
    my ($rnaeval_out_optimal)  = shift @args;
    my ($rnaeval_out_stemloop) = shift @args;

    open( my $STEM_FH, '<', $rnastemloop_out_stemloop ) or die $!;
    open( my $EVAL_OPT_FH, '<', $rnaeval_out_optimal ) or die $!;
    open( my $EVAL_STEM_FH, '<', $rnaeval_out_stemloop ) or die $!;
    my $msg = "Processing RNAstemloop ( $sequence_dir, $strand, $seq_length, $rnastemloop_out_stemloop, $rnaeval_out_optimal, $rnaeval_out_stemloop )";
    debug( $msg, miRkwood->DEBUG() );
    my $res =
      process_RNAstemloop( $sequence_dir, $strand, $seq_length, $STEM_FH,
        $EVAL_OPT_FH, $EVAL_STEM_FH );
    close($STEM_FH);
    close($EVAL_OPT_FH);
    close($EVAL_STEM_FH);
    return $res;
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

		if ( miRkwood::Utils::is_fasta_header( $stem_line )) {
			$nameSeq = substr ($stem_line, 1, -1);
		}
		elsif ( miRkwood::Utils::is_fasta_line_relaxed($stem_line ) ) {
			$dna = substr $stem_line, 0, -1;
			$line_eval_opt = substr( <$EVAL_OPT_FH>, 0, -1 );    # the sequence as well
			if ( miRkwood::Utils::is_fasta_header( $line_eval_opt ) ) {
			    $line_eval_opt = substr( <$EVAL_OPT_FH>, 0, -1 );
			}
			$line_eval_stem = substr( <$EVAL_STEM_FH>, 0, -1 );    # the sequence as well
	         if ( miRkwood::Utils::is_fasta_header( $line_eval_stem )) {
                $line_eval_stem = substr( <$EVAL_STEM_FH>, 0, -1 );
            }
			if ( $dna ne $line_eval_opt || $dna ne $line_eval_stem ) {
				warn ('The sequences differ in RNAeval and RNAstemloop output');
			}
		}
		elsif ( ( $stem_line =~ /(.*)/ ) ) {
			$structure_stemloop = $1;
			$line_eval_opt = <$EVAL_OPT_FH>;    # the structure as well, and the energy
			$line_eval_stem = <$EVAL_STEM_FH>;
			
			my ( $structure_optimal, $energy_optimal ) =
                miRkwood::Parsers::parse_Vienna_line($line_eval_opt);
			my ( $structure_stemloop, $energy_stemloop ) =
                miRkwood::Parsers::parse_Vienna_line($line_eval_stem);
			if ($structure_optimal)
			{                        # We have a structure

				if ( $nameSeq =~ /.*__(\d*)-(\d*)$/ ) {
					my ($mfei, $amfe) =
					  miRkwood::Utils::compute_mfei_and_amfe( $dna, $energy_optimal );
					my ( $start, $end );
					if ( $strand eq '-' ) {
						my $res =
						  miRkwood::Utils::get_position_from_opposite_strand
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
						"amfe"      => $amfe,
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
	my $start     = shift @args;
	my $end       = shift @args;
	my $ref_start = shift @args;
	my $ref_end   = shift @args;
	( $start ne '' && $end ne '' && $ref_start ne '' && $ref_end ne '' ) or die('Not enough values provided');
	$ref_start <= $start or die("Positions should be ordered : $ref_start <= $start");
	return ( $start < ( $ref_start + $ref_end ) / 2 );
}

=method is_included

Test whether one candidate is included into the other
(based on their positions).

=cut

sub is_included {
	my @args      = @_;
    my $start     = shift @args;
    my $end       = shift @args;
    my $ref_start = shift @args;
    my $ref_end   = shift @args;
    ( $start ne '' && $end ne '' && $ref_start ne '' && $ref_end ne '' ) or die('Not enough values provided');
	$ref_start <= $start or die("Positions should be ordered : $ref_start <= $start");
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

    debug( "Writing candidate information in $candidate_dir", miRkwood->DEBUG() );
	#Writing seq.txt
	my $candidate_sequence = File::Spec->catfile( $candidate_dir, 'seq.txt' );
	open( my $SEQ_FH, '>', $candidate_sequence )
	  or die "Error when opening $candidate_sequence: $!";
	print $SEQ_FH ">$candidate{'name'}\n$candidate{'dna'}\n";
	close $SEQ_FH;

	#Writing sequence information
	my $seq_info_file = File::Spec->catfile( $candidate_dir, 'sequence_information.txt' );
	open( my $SEQ_INFO_FH, '>', $seq_info_file )
	  or die "Error when opening $seq_info_file: $!";
	print $SEQ_INFO_FH $candidate{'strand'} . "\t" .  $candidate{'start'} . "\t" . $candidate{'end'};
	close $SEQ_INFO_FH;

	process_outRNAFold( $candidate_dir, 'optimal', $candidate{'name'},
		$candidate{'dna'}, $candidate{'structure_optimal'}, $candidate{'energy_optimal'} );
	process_outRNAFold( $candidate_dir, 'stemloop', $candidate{'name'},
		$candidate{'dna'}, $candidate{'structure_stemloop'}, $candidate{'energy_stemloop'} );

	# Writing energy file
	my $energy_file = File::Spec->catfile( $candidate_dir, 'outMFEI.txt' );
    open( my $ENERGY_FH, '>', $energy_file )
        or die "Unable to open $energy_file: $!";
    my $content = $candidate{'name'} . "\t" . $candidate{'mfei'} . "\t" . $candidate{'energy_optimal'} . "\t" . $candidate{'amfe'};
    print $ENERGY_FH $content;
    close $ENERGY_FH or die "Unable to close: $!";

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

    # Writing VARNA image
    my $cfg = miRkwood->CONFIG();

    if ( $cfg->param('options.varna') ) {
        my $varna_image = File::Spec->catfile( $candidate_dir, 'image.png' );
        debug( "Generating image using VARNA in $varna_image", miRkwood->DEBUG() );
        miRkwood::Programs::run_varna_on_structure( $candidate{'dna'}, $candidate{'structure_stemloop'}, $varna_image )
          or carp('Problem during image generation using VARNA');
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
	debug( "A posteriori tests in $job_dir", miRkwood->DEBUG() );
    my $candidates_dir = File::Spec->catdir( $job_dir, 'candidates' );
    my $workspace_dir = miRkwood::Paths->get_workspace_path($job_dir);

	my @sequence_dirs = miRkwood::FileUtils::get_dirs_from_directory($workspace_dir);

	mkdir $candidates_dir;
	foreach my $dir (@sequence_dirs)
	{
		my $sequence_dir = File::Spec->catdir( $workspace_dir, $dir );
		debug( "Entering sequence $sequence_dir", miRkwood->DEBUG() );

		my @candidate_dirs = miRkwood::FileUtils::get_dirs_from_directory($sequence_dir);
		foreach my $subDir (@candidate_dirs) {
			my $candidate_dir =
			  File::Spec->catdir( $sequence_dir, $subDir );

			debug( "Entering candidate $subDir", miRkwood->DEBUG() );
			process_tests_for_candidate( $candidate_dir, $subDir );
			debug( "Done with candidate $subDir", miRkwood->DEBUG() );
			if (
				!eval {
					miRkwood::CandidateHandler
					  ->serialize_candidate_from_run( $job_dir, $dir,
						$subDir, $candidates_dir );
				}
			  )
			{
				# Catching
				carp( "Serialization of $subDir failed" );
			}
			else {
				debug( "Done with serializing $subDir", miRkwood->DEBUG() );
				# All is well
			}
		} # foreach my $file (@files)
		debug( "Done with initial sequence $dir", miRkwood->DEBUG() );
	} # foreach my $dir (@dirs)
	debug( "Done with all the tests", miRkwood->DEBUG() );
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

    my $cfg = miRkwood->CONFIG();

	####calcul p-value randfold
	if ( $cfg->param('options.randfold') ) {
		debug( "Running test_randfold on $seq_file", miRkwood->DEBUG() );
		miRkwood::PosterioriTests::test_randfold( $candidate_dir,
			$seq_file );
	}

	if ( $cfg->param('options.align') ) {
		debug( "Running test_alignment on $candidate_rnafold_stemploop_out", miRkwood->DEBUG() );
		my $file_alignement =
		  miRkwood::PosterioriTests::test_alignment( $candidate_dir,
			$candidate_rnafold_stemploop_out );
		post_process_alignments( $candidate_dir,
			$candidate_rnafold_stemploop_out,
			$file_alignement );
	}

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
	  miRkwood::Components::get_data_from_rnafold_out(
		$candidate_rnafold_stemploop_out);
	my ( $name, $position, $DNASequence, $Vienna ) = @res;
	my %alignments;

    if (-z $file_alignement){
        return;
    }
	if (
		!eval {
			%alignments =
			  miRkwood::Components::parse_custom_exonerate_output(
				$file_alignement);
		}
	  )
	{
	    # Catching exception
        carp("Exception when parsing exonerate output $file_alignement");
        return;
	}
	else {
		%alignments =
		  miRkwood::Components::merge_alignments( \%alignments );
		my $tmp_file =
		  File::Spec->catfile( $candidate_dir, "mirdup_validation.txt" );
		my %mirdup_results =
		  miRkwood::MiRdup->validate_with_mirdup( $tmp_file, $name,
			$DNASequence, $Vienna, keys %alignments );
		my $mirdup_results_file =
		  File::Spec->catfile( $candidate_dir, 'mirdup_results.yml' );
		YAML::XS::DumpFile( $mirdup_results_file, %mirdup_results );

		my $alignments_results_file =
		  File::Spec->catfile( $candidate_dir, 'merged_alignments.yml' );
		YAML::XS::DumpFile( $alignments_results_file, %alignments );
	}
}

=method mark_job_as_finished

Mark the current job as finished

=cut

sub mark_job_as_finished {
    my ( @args ) = @_;
    my $job_dir     = shift @args;
    my $is_finished_file = File::Spec->catfile( $job_dir, 'finished' );
    open( my $finish, '>', $is_finished_file )
        or die "Error when opening $is_finished_file: $!";
    close $finish;
    return (-e $is_finished_file);
}

1;
