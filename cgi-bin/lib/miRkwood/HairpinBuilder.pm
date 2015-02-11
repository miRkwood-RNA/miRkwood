package miRkwood::HairpinBuilder;

use strict;
use warnings;
use POSIX;

use miRkwood;

use miRkwood::Utils;

use List::Util qw(max min);


use miRkwood::Parsers;
use miRkwood::Programs;
use miRkwood::Utils;
use miRkwood::ClusterJobSebastien;

use File::Spec;
use File::Basename;
use Log::Message::Simple qw[msg error debug];

sub new {
	my ($class, $genome_db, $workspace, $parsed_bed) = @_;
	return bless {
        'genome_db'  => $genome_db,
        'workspace'  => $workspace,
        'parsed_bed' => $parsed_bed
    }, $class;
}

# gets the sequence. start starts at 0. end is excluded
sub get_sub_sequence {
	my ($this, $chr, $start, $end) = @_;
	# Chr1, 1, 1
	return substr($this->{genome_db}{$chr}, $start, $end-$start);
}


# gets the sequence. start starts at 0. end is excluded
# strand is 1 or -1
# if strand == -1, then the reverse complement is returned
sub get_sub_sequence_on_strand {
	my ($this, $chr, $start, $end, $strand) = @_;
	my $seq = $this->get_sub_sequence($chr, $start, $end);
	if ($strand eq '-') {
		return miRkwood::Utils::reverse_complement($seq);
	}
	return $seq;
}

sub build_hairpins {
	my $this = shift;
	my $locus = shift;
	
	my $chr = $locus->{'chr'};
	#~ my $strand = $locus->{strand} eq '+' ? 1 : -1;
	#~ my $seq_id = $chr . '__' . $locus->{begin}+1 . '-' . ($locus->{end});
	
	my $working_dir = File::Spec->catdir($this->{'workspace'}, $chr);
	mkdir $working_dir;
	$working_dir = File::Spec->catdir($working_dir, $locus->{begin} . '-' . ($locus->{end}-1) . $locus->{'strand'});
	mkdir $working_dir;

	my $rnalfold_output_filename = 'rnalfold_out';

	my $genomic_seq = $this->get_sub_sequence_on_strand($chr, $locus->{begin}, $locus->{end}, $locus->{'strand'});

	my $rnalfold_output_filepath = run_rnalfold('miRnaPrecursor', $genomic_seq, $working_dir, $rnalfold_output_filename);

	my $output_stemloop_opt = File::Spec->catfile($working_dir, 'rnastemloop_optimal');
	my $output_stemloop = File::Spec->catfile($working_dir, 'rnastemloop_stemloop');

	run_rnastemloop($rnalfold_output_filepath, $output_stemloop, $output_stemloop_opt);

	my $rnaeval_out_optimal = $this->run_RNAeval_on_RNAstemloop_optimal_output($output_stemloop_opt);
	my $rnaeval_out_stemloop = $this->run_RNAeval_on_RNAstemloop_stemloop_output($output_stemloop);

	my $new_candidates = $this->process_RNAstemloop_on_filenames($output_stemloop, $rnaeval_out_optimal, $rnaeval_out_stemloop,
		$locus->{begin}, $locus->{end}-$locus->{begin}, $chr, $locus->{strand}, $locus->{miRnas}, 
		$this->{'parsed_bed'}
		);

	#~ ($new_candidates, $already_positioned) = filter_candidates_on_position($new_candidates, $already_positioned);

	my @sorted_new_candidates = sort { $a->{'start_position'} <=> $b->{'start_position'} } @{$new_candidates};
	
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
    miRkwood::Programs::run_rnaeval( $rnastemloop_out, $rnaeval_out ) or die('Problem when running RNAeval');

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
    my $parsed_bed = shift @args;

    open( my $STEM_FH, '<', $rnastemloop_out_stemloop ) or die "Error opening $rnastemloop_out_stemloop: $!";
    open( my $EVAL_OPT_FH, '<', $rnaeval_out_optimal ) or die $!;
    open( my $EVAL_STEM_FH, '<', $rnaeval_out_stemloop ) or die $!;
    my $msg = "     Processing RNAstemloop ( $rnastemloop_out_stemloop, $rnaeval_out_optimal, $rnaeval_out_stemloop )";
    debug( $msg, miRkwood->DEBUG() );
    my $candidates = $self->process_RNAstemloop($STEM_FH, $EVAL_OPT_FH, $EVAL_STEM_FH, $sequence_begin, $seq_len, $chr, $strand,
    $sequence_miRnas 
    ,$parsed_bed
    );
    close($STEM_FH);
    close($EVAL_OPT_FH);
    close($EVAL_STEM_FH);
    return $candidates;
}

# See std::lower_bound in the C++ standard library
sub lowerbound_binsearch {
	my $arr = shift;
	my $first = shift;
	my $last = shift; # last is excluded, past the last item
	my $val = shift;
	
	my $len = $last - $first;
	while ($len > 0) {
		my $half = $len >> 1;
		my $middle = $first + $half;
		if ($arr->[$middle]{'begin'} < $val) {
			$first = $middle;
			$first++;
			$len = $len - $half - 1;
		}
		else {
			$len = $half;
		}
	}
	return $first;
}

# See std::upper_bound in the C++ standard library
sub upperbound_binsearch {
	my $arr = shift;
	my $first = shift;
	my $last = shift; # last is excluded, past the last item
	my $val = shift;
	
	my $it;
	my ($count, $step);
	$count = $last - $first;
	while ($count > 0)	{
		$it = $first; $step=$count/2; $it += $step;
		if (!($val < $arr->[$it]{'begin'})) {
			$it++;
			$first = $it;
			$count -= $step + 1;
		}
		else {
			$count= $step;
		}
	}
	return $first;
}

sub get_contained_reads {
	my ($parsed_bed, $chr, $region_begin, $region_end, $strand) = @_;
	my $reads = $parsed_bed->{$chr}{$strand};

	my $low = lowerbound_binsearch($reads, 0, scalar @{$reads}, $region_begin);
	my $high = upperbound_binsearch($reads, $low, scalar @{$reads}, $region_end);

	my %result = ();
	if ($low == scalar @{$reads}) {
		return \%result;
	}
	#~ if ($high == scalar @{$reads}) {
		#~ $high--;
	#~ }

	for (my $i = $low; $i < $high; $i++) {
		my $read_begin = $reads->[$i]{'begin'}-1;
		my @read_ends = keys %{$reads->[$i]{'ends'}};
		foreach my $read_end (@read_ends) {
			if ($read_end <= $region_end) {
				$result{"$read_begin-".($read_end-1)} = $reads->[$i]{'ends'}{$read_end};
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
sub eval_single_stemloop {
	my $stemloop = shift;
	my $miRnaPos = shift;
	
	if (!defined ($miRnaPos)) {
		return 1;
	}

	foreach my $miRnas ( @{$miRnaPos} ) {
		if ($miRnas->{source} == miRkwood::ClusterJobSebastien::DUE_TO_TWO_SPIKES) {
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
	my ($parsed_bed) = shift @args;

	my ($line_eval_opt, $line_eval_stem);
	my ($nameSeq, $dna, $structure_stemloop);
	my @candidates_array = ();

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
				warn ('The sequences differ in RNAeval and RNAstemloop output');
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
					if (eval_single_stemloop($stemloop, $sequence_miRnas) == 1) {
						my $cluster_position = $seq_begin. '-' . ($seq_begin+$seq_len-1);
						my $res = {
							'name' => $chr. '__' .($stemloop->{'begin'}+1).'-'.$stemloop->{'end'} . $strand,
							'strand' => $strand,
							'sequence' => $dna,
							'start_position' => $stemloop->{'begin'}+1, # starts at 1
							'end_position' => $stemloop->{'end'}, # includes the end
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
    miRkwood::Programs::run_rnalfold($seq_name, $seq, $temp_file, $rnalfold_output) or die("Problem when running RNALfold: $!");
    return $rnalfold_output;
}


=method run_rnastemloop

Static private helper function. You shouldnt use this function.

=cut
sub run_rnastemloop {
	my ( $input, $output_stemloop, $output_optimal ) = @_;
	return miRkwood::Programs::run_rnastemloop($input, $output_stemloop, $output_optimal);
}

1;
