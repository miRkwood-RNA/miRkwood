package miRkwood::HairpinBuilder;

use strict;
use warnings;
use POSIX;

use List::Util qw(max min);

use miRkwood;
use miRkwood::Paths;
use miRkwood::Parsers;
use miRkwood::Programs;
use miRkwood::Utils;
use miRkwood::ClusterJobSebastien;
use miRkwood::CppBinarySearch;

use File::Spec;
use File::Basename;
use Log::Message::Simple qw[msg error debug];

sub new {
	my ($class, $genome_db, $workspace, $parsed_bed
	#STATS BEG
	#~ ,$stats
	#STATS END
	) = @_;
	my $this = bless {
        'genome_db'  => $genome_db,
        'workspace'  => $workspace,
        'parsed_bed' => $parsed_bed,
        # THIS DISCARDS HAIRPIN CANDIDATES THAT CONTAIN STRICTLY LESS THAN 'read_coverage_threshold' READS
        'read_coverage_threshold' => 2, # TO MODIFY IN RELEASE
    #STATS BEG
		#~ 'stats' => $stats
    #STATS END
    }, $class;

    return $this;
}

sub get_parsed_bed {
	my $this = shift;
	return $this->{'parsed_bed'};
}

# gets the sequence. start starts at 0. end is excluded
# strand is '+' or '-'
# if strand == '-', then the reverse complement is returned
sub get_sub_sequence_on_strand {
	my ($this, $chr, $start, $end, $strand) = @_;
	my $seq = $this->{'genome_db'}->seq($chr, $start+1, $end, $strand eq '+' ? 1 : -1);
	#~ $seq =~ s/[RYSWKMBDHV]/N/ig;
	return $seq;
}

sub build_hairpins {
	my $this = shift;
	my $locus = shift;

	#~ warn("HairpinBuilder: $locus->{'chr'}:", $locus->{begin}+1, '-', $locus->{end}, ' [', $locus->{strand}, "]\n");
	
	my $chr = $locus->{'chr'};
	
    my $working_chr_dir = miRkwood::Paths::get_workspace_chromosome_dir( $this->{'workspace'}, $chr );
	mkdir $working_chr_dir;

    my $working_dir = miRkwood::Paths::get_workspace_candidate_dir( $this->{'workspace'},
                                                                    $chr,
                                                                    ($locus->{begin}+1) . '-' . ($locus->{end}),
                                                                    $locus->{strand} );
	mkdir $working_dir;

	my $rnalfold_output_filename = 'rnalfold_out';

	my $genomic_seq = $this->get_sub_sequence_on_strand($chr, $locus->{begin}, $locus->{end}, $locus->{'strand'});

	my $rnalfold_output_filepath = run_rnalfold('miRnaPrecursor', $genomic_seq, $working_dir, $rnalfold_output_filename);

	my $output_stemloop_opt = File::Spec->catfile($working_dir, 'rnastemloop_optimal');
	my $output_stemloop = File::Spec->catfile($working_dir, 'rnastemloop_stemloop');

	#~ if (!-e $output_stemloop) {
		miRkwood::Programs::run_rnastemloop($rnalfold_output_filepath, $output_stemloop, $output_stemloop_opt)
		or die("Problem running RNAstemloop [Input: '$rnalfold_output_filepath', Outputs: '$output_stemloop', '$output_stemloop_opt'");
	#~ }
	#~ else {
		#~ print "\tRNAstemloop output already exists\n";
	#~ }

	my $rnaeval_out_optimal = $this->run_RNAeval_on_RNAstemloop_optimal_output($output_stemloop_opt);
	my $rnaeval_out_stemloop = $this->run_RNAeval_on_RNAstemloop_stemloop_output($output_stemloop);

	my $new_candidates = $this->process_RNAstemloop_on_filenames($output_stemloop, $rnaeval_out_optimal, $rnaeval_out_stemloop,
		$locus->{begin}, $locus->{end}-$locus->{begin}, $chr, $locus->{strand}, $locus->{miRnas});

	#~ ($new_candidates, $already_positioned) = filter_candidates_on_position($new_candidates, $already_positioned);

	my @sorted_new_candidates = sort { $a->{'start_position'} <=> $b->{'start_position'} } @{$new_candidates};

    miRkwood::Utils::display_var_sizes_in_log_file( '..... HairpinBuilder : build_hairpins' );

	return \@sorted_new_candidates;
}


=method run_RNAeval_on_RNAstemloop_optimal_output

=cut

sub run_RNAeval_on_RNAstemloop_optimal_output {
    my ($self, @args)  = @_;
    my $rnastemloop_out_optimal = shift @args;
    return $self->run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_optimal, 'optimal' );
}

=method run_RNAeval_on_RNAstemloop_stemloop_output

=cut

sub run_RNAeval_on_RNAstemloop_stemloop_output {
    my ($self, @args)  = @_;
    my $rnastemloop_out_stemloop = shift @args;
    return $self->run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_stemloop, 'stemloop' );
}

=method run_RNAeval_on_RNAstemloop_output


=cut

sub run_RNAeval_on_RNAstemloop_output {
    my ( $self, @args ) = @_;
    my ( $rnastemloop_out, $suffix ) = @args    ;
    my $current_sequence_dir = dirname($rnastemloop_out);
    debug( "     Processing RNAstemloop output for $suffix $rnastemloop_out", miRkwood->DEBUG() );
    my $rnaeval_out = File::Spec->catfile( $current_sequence_dir, "rnaeval_$suffix.out" );

    debug( "     Running RNAeval in $rnaeval_out", miRkwood->DEBUG() );
    
    #~ if (!-e $rnaeval_out) {
		miRkwood::Programs::run_rnaeval( $rnastemloop_out, $rnaeval_out ) or 
		die("Problem when running RNAeval [Input file: '$rnastemloop_out' Output file: '$rnaeval_out'");
	#~ }
    return $rnaeval_out;
}

sub process_RNAstemloop_on_filenames {
    my ($self, @args)  = @_;
    my ($rnastemloop_out_stemloop)      = shift @args;
    my ($rnaeval_out_optimal)  = shift @args;
    my ($rnaeval_out_stemloop) = shift @args;
    my ($sequence_begin) = shift @args;
    my ($seq_len) = shift @args;
    my $chr = shift @args;
    my ($strand) = shift @args;
    my ($sequence_miRnas) = shift @args;

    open( my $STEM_FH, '<', $rnastemloop_out_stemloop ) or die "Error opening $rnastemloop_out_stemloop: $!";
    open( my $EVAL_OPT_FH, '<', $rnaeval_out_optimal ) or die $!;
    open( my $EVAL_STEM_FH, '<', $rnaeval_out_stemloop ) or die $!;
    my $msg = "     Processing RNAstemloop ( $rnastemloop_out_stemloop, $rnaeval_out_optimal, $rnaeval_out_stemloop )";
    debug( $msg, miRkwood->DEBUG() );
    my $candidates = $self->process_RNAstemloop($STEM_FH, $EVAL_OPT_FH, $EVAL_STEM_FH, $sequence_begin, $seq_len, $chr, $strand,
    $sequence_miRnas);
    close($STEM_FH);
    close($EVAL_OPT_FH);
    close($EVAL_STEM_FH);
    return $candidates;
}

sub get_contained_read_coverage {
	my ($parsed_bed, $chr, $region_begin, $region_end, $strand) = @_;
	my $reads = $parsed_bed->{$chr}{$strand};

	my $low = miRkwood::CppBinarySearch::lower_bound($reads, 0, scalar @{$reads}, $region_begin);
	if ($low == scalar @{$reads}) {
		return 0;
	}
	my $high = miRkwood::CppBinarySearch::upper_bound($reads, $low, scalar @{$reads}, $region_end);

	my $cov = 0;

	for (my $i = $low; $i < $high; $i++) {
		my $read_begin = $reads->[$i]{'begin'};
		my @read_ends = keys %{$reads->[$i]{'ends'}};
		foreach my $read_end (@read_ends) {
			if ($read_end <= $region_end) {
				$cov += $reads->[$i]{'ends'}{$read_end};
			}
		}
	}

	return $cov;
}

sub get_contained_reads {
	my ($parsed_bed, $chr, $region_begin, $region_end, $strand) = @_;
	my $reads = $parsed_bed->{$chr}{$strand};

	my $low = miRkwood::CppBinarySearch::lower_bound($reads, 0, scalar @{$reads}, $region_begin);
	if ($low == scalar @{$reads}) {
		return {};
	}
	my $high = miRkwood::CppBinarySearch::upper_bound($reads, $low, scalar @{$reads}, $region_end);

	my %result = ();

	for (my $i = $low; $i < $high; $i++) {
		my $read_begin = $reads->[$i]{'begin'};
		my @read_ends = keys %{$reads->[$i]{'ends'}};
		foreach my $read_end (@read_ends) {
			if ($read_end <= $region_end) {
				$result{($read_begin+1).'-'.$read_end} = $reads->[$i]{'ends'}{$read_end};
			}
		}
	}

	return \%result;
}


=method raw_regions_intertect

Static private helper function. You shouldnt use this function.

=cut
sub raw_regions_intertect {
	my ($a, $b) = @_;
	if ($a->{begin} <= $b->{begin} && $a->{end} > $b->{begin}) {
		return 1;
	}
	elsif ($b->{begin} <= $a->{begin} && $b->{end} > $a->{begin}) {
		return 1;
	}
	return 0;
}


=method eval_stemloop_output

Static private helper function. You shouldnt use this function.

=cut
sub test_miRna_intersection {
	my $stemloop = shift;
	my $miRnaPos = shift;
	
	if (!defined ($miRnaPos)) {
		return 1;
	}
	
	foreach my $miRnas ( @{$miRnaPos} ) {
		if ($miRnas->{source} == miRkwood::ClusterJobSebastien->DUE_TO_TWO_SPIKES ) {
			if (raw_regions_intertect($stemloop, $miRnas->{first}) && raw_regions_intertect($stemloop, $miRnas->{second})) {
				return 1;
			}
		}
		else {
			if (raw_regions_intertect($stemloop, $miRnas->{from_read})) {
				return 1;
			}
		}
	}
	return 0;
}


sub eval_single_stemloop {
	my $this = shift;
	my $chr = shift;
	my $strand = shift;
	my $stemloop = shift;
	my $miRnaPos = shift;
	
	if (test_miRna_intersection($stemloop, $miRnaPos) == 0) {
		return 0;
	}
	
	my $parsed_bed = $this->get_parsed_bed;
	
	# No bed supplied
	if (!defined ($parsed_bed)) {
		return 1;
	}
	if (get_contained_read_coverage($parsed_bed, $chr, $stemloop->{'begin'}, $stemloop->{'end'}, $strand) >= 
				$this->{'read_coverage_threshold'}) {
		return 1;
	}
	else {
		return 0;
	}
}

=method process_RNAstemloop

Process the outputs of RNAstemloop + RNAeval (both optimal and stemloop)
Returns a list of candidates.

=cut

sub process_RNAstemloop {
	my ($self, @args)  = @_;
	my ($STEM_FH)      = shift @args;
	my ($EVAL_OPT_FH)  = shift @args;
	my ($EVAL_STEM_FH) = shift @args;
	my ($seq_begin) = shift @args;
	my ($seq_len) = shift @args;
	my ($chr) = shift @args;
	my ($strand) = shift @args;
	my ($sequence_miRnas) = shift @args;
	my ($parsed_bed) = $self->get_parsed_bed;

	my ($line_eval_opt, $line_eval_stem);
	my ($nameSeq, $dna, $structure_stemloop);
	my @candidates_array = ();
	
	my $DEBUG_NAME = $chr. ':' . ($seq_begin+1) . '-' . ($seq_begin+$seq_len) . ' [' . $strand . ']';

	while ( my $stem_line = <$STEM_FH> ) {
		if (miRkwood::Utils::is_fasta_header($stem_line)) {
			$nameSeq = substr ($stem_line, 1, -1);
		}
		elsif (miRkwood::Utils::is_fasta_line_relaxed($stem_line)) {
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
				warn ('The sequences differ in RNAeval and RNAstemloop output for '. $DEBUG_NAME);
			}
		}
		elsif ($stem_line =~ /(.*)/) {
			$structure_stemloop = $1;
			$line_eval_opt = <$EVAL_OPT_FH>;    # the structure as well, and the energy
			$line_eval_stem = <$EVAL_STEM_FH>;

			my ( $structure_optimal, $energy_optimal ) =
				miRkwood::Parsers::parse_Vienna_line($line_eval_opt);
			my ( $structure_stemloop, $energy_stemloop ) =
				miRkwood::Parsers::parse_Vienna_line($line_eval_stem);
			if ($structure_optimal) { # We have a structure
				if ($nameSeq =~ /.*__(\d*)-(\d*)$/) {
					my ($mfei, $amfe) = miRkwood::Utils::compute_mfei_and_amfe($dna, $energy_optimal);
					#~ my ($start, $end);
                    my $stemloop = {};
					if ($strand eq '-') {
						#~ ($start, $end) = @{ miRkwood::Utils::get_position_from_opposite_strand( $1, $2, $seq_len) };
                        $stemloop = {begin => $seq_len - $2 + $seq_begin, end => $seq_len - $1 +1 + $seq_begin};
					}
					else {
						#~ ($start, $end) = ($1, $2);
                        $stemloop = {begin => $1 + $seq_begin-1, end => $2 + $seq_begin};
					}
					if ($self->eval_single_stemloop($chr, $strand, $stemloop, $sequence_miRnas) == 1) {
						my $cluster_position = ($seq_begin+1). '-' . ($seq_begin+$seq_len);
						my $res = {
							'name' => $chr. '__' .($stemloop->{'begin'}+1).'-'.$stemloop->{'end'} . $strand,
							'strand' => $strand,
							'sequence' => $dna,
							'start_position' => $stemloop->{'begin'}+1, # 1-based
							'end_position' => $stemloop->{'end'}, # excludes the end
							'mfei' => $mfei,
							'amfe' => $amfe,
							'structure_optimal' => $structure_optimal,
							'energy_optimal' => $energy_optimal,
							'structure_stemloop' => $structure_stemloop,
							'energy_stemloop' => $energy_stemloop,
							'reads' => get_contained_reads($parsed_bed, $chr, $stemloop->{'begin'}, $stemloop->{'end'}, $strand),
							'cluster' => $cluster_position
						};
						push @candidates_array, $res;
					}
					
    #STATS BEG
					#~ else {
						#~ $self->update_stats_for_discarded_miRNA($sequence_miRnas);
					#~ }
    
    #STATS END
				}
			}
			else {
				warn( "No structure found in $line_eval_opt for $DEBUG_NAME" );
				#STATS BEG
				#~ $self->update_stats_for_discarded_miRNA($sequence_miRnas);
				#STATS END
			}    # if $line2
		}
		else {
			warn( "Unrecognized line for $DEBUG_NAME" );
		}    #if $line1
	}    #while $line=<IN>
	return \@candidates_array;
}

=method run_rnalfold

Static private helper function. You shouldnt use this function.

=cut
sub run_rnalfold {
	my $seq_name = shift;
	my $seq = shift;
	my $output_folder = shift;
	my $output_file = shift;
	my $rnalfold_output = File::Spec->catfile($output_folder, $output_file);
    my $temp_file = File::Spec->catfile($output_folder, $output_file . '_sequence.txt');
    #~ if (!-e $rnalfold_output) {
		miRkwood::Programs::run_rnalfold($seq_name, $seq, $temp_file, $rnalfold_output) or die("Problem when running RNALfold [$temp_file]: $!");
	#~ }
    return $rnalfold_output;
}


#STATS BEG
    
#~ sub update_stats_for_discarded_miRNA {
	#~ my $this = shift;
	#~ my $sequence_miRnas = shift;
	#~ 
	#~ foreach my $miRnas ( @{$sequence_miRnas} ) {
		#~ if ($miRnas->{source} == miRkwood::ClusterJobSebastien::DUE_TO_TWO_SPIKES) {
			#~ $this->{'stats'}{'COUPLE_DISCARDED'}++;
		#~ }
		#~ else {
			#~ $this->{'stats'}{'SINGLE_DISCARDED'}++;
		#~ }
	#~ }
#~ }
    
#STATS END
1;
