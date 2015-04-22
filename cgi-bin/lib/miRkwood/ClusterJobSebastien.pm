package miRkwood::ClusterJobSebastien;

use strict;
use warnings;
use POSIX;

use miRkwood;
use miRkwood::MiRnaDuplexDetector;
use miRkwood::HairpinBuilder;
use miRkwood::Utils;


#~ use List::BinarySearch qw( binsearch  binsearch_pos  binsearch_range );
use List::Util qw(max min);

use constant {
	DUE_TO_SINGLE_SPIKE => 1,
	DUE_TO_TWO_SPIKES => 2,

	ISOLATED_SPIKE => 3,
	COUPLE_OF_SPIKES => 4,
	MULTIPLE_SPIKES => 5
};

use miRkwood::Parsers;
use miRkwood::Programs;
use miRkwood::Utils;
use miRkwood::SequenceJob;

use File::Spec;
use File::Basename;
use Log::Message::Simple qw[msg error debug];

=method new

Constructor

 Usage : my $cluster_job = ClusterJob->new($workspace_dir)
 Input : The genome file
 Return: The cunstructed instance

=cut

sub new {
    my ( $class,
    $genome_db ) = @_;
    my $self = bless {
        genome_db => $genome_db
    }, $class;
    return $self;
}


sub init_from_clustering {
	my ($this, $clustering) = @_;
	$this->{chr_info} = $clustering->{'chr_info'};
	$this->{accepting_time} = $clustering->{'accepting_time'};
    return;
}


# gets the sequence. start starts at 0. end is excluded
sub get_sub_sequence {
	my ($this, $chr, $start, $end) = @_;
	return $this->{'genome_db'}->seq($chr, $start+1, $end, 1);
}


sub display_chr_coordinatees {
	my ($chr, $begin, $end) = @_;
	return $chr. ':' . ($begin+1) . '-' . $end;
}


sub extract_spike_train_per_chr {
	my $this = shift;
	my $trains_for_chr = shift;

	my @spike_array = ();
	foreach my $train (@{$trains_for_chr}) {
		foreach my $spike (@{$train->{spikes}}) {
			push @spike_array, { 'spike' => $spike };
		}
	}

	my $accepting_time = $this->{'accepting_time'};
	my $min_length = 21;
	my $min_length_for_neighbor = 40;

	for (my $i = 0, my $e = scalar @spike_array; $i < $e; $i++) {
		my $current_spike = $spike_array[$i]{'spike'};
		my @neighbors = ();
		my $unknown_strand = $current_spike->{'strand'} eq '?';

		for (my $j = $i+1; $j < $e; $j++) {
			my $neighbor = $spike_array[$j]{'spike'};
			if ($neighbor->{'begin'} - $current_spike->{'end'} < $accepting_time) {
				if ($unknown_strand || $neighbor->{'strand'} eq $current_spike->{'strand'} || $neighbor->{'strand'} eq '?') {
					push @neighbors, $neighbor;
				}
			}
			else {
				last;
			}
		}
		$spike_array[$i]{'neighbors'} = \@neighbors;
		$spike_array[$i]{'class'} = (scalar(@neighbors) == 0 ? ISOLATED_SPIKE : scalar(@neighbors) == 1 ? COUPLE_OF_SPIKES : MULTIPLE_SPIKES);
	}

	return \@spike_array;
}

=method process_window_spikes

Retrieve candidates miRNA from windows spike analysis.

 Usage : my $miRnaPos = $self->process_window_spikes($windows),
 Input : The windows returned by get_windows
 Return: A hash ref {
				chr => [array of miRNA candidates (hash ref) {
					strand => '+' or '-',
					source => DUE_TO_SINGLE_SPIKE or DUE_TO_TWO_SPIKES,
					first => First miRNA candidate (lower chr coordinates). Hash ref {
						begin => start position (1-based)
						end => end position (excluded)
						strand => '+' or '-' #This parameter may be inaccurate, you should rely on the top-level strand key instead
						}
					second => Second miRNA candidate (duplex) (larger chr coordinates). Hash ref {
						begin => start position (1-based)
						end => end position (excluded)
						strand => '+' or '-' #This parameter may be inaccurate, you should rely on the top-level strand key instead
						}
					from_read => miRNA candidate that was deduced from a read spike. Hash ref on 'first' or 'second'. This key exists only if 'source' is DUE_TO_SINGLE_SPIKE
					detected => miRNA candidate that was deduced by folding the sequence. Hash ref on 'first' or 'second'. Doesnt point on the same miRNA than from_read.
						This key exists only if 'source' is DUE_TO_SINGLE_SPIKE
				}]
			}

=cut
sub process_spikes {
	my $this = shift;
	my $spikes_per_chr = shift;
	my %miRnaPos = ();

	#STATS BEG
	#~ my %stats = (
	#~ 'ISOLATED_SPIKE_DISCARDED' => 0,
	#~ 'COUPLE_OF_SPIKES_DISCARDED' => 0,
	#~ 'COUPLE_OF_SPIKES_TO_ISOLATED_SPIKE' => 0
	#~ );
	#STATS END

	$this->{detector} = MiRnaDuplexDetector::MiRnaDetector->new(5000);

	foreach my $chr (keys %{ $this->{chr_info} }) {
		$miRnaPos{$chr} = $this->process_spikes_for_chr($chr, $spikes_per_chr->{$chr}
		#STATS BEG
		#~ , \%stats
		#STATS END
		);
	}
	return \%miRnaPos
	#STATS BEG
	#~ , \%stats
	#STATS END
	;
}

=method __enlarged_spike

Static private helper function. You shouldnt use this function.

=cut
sub __enlarged_spike {
	my ($spike, $min_length, $chr_length) = @_;
	my $spike_length = $spike->{'end'} - $spike->{'begin'};
	if ($spike_length < $min_length) {
		my $new_spike = {begin => $spike->{'begin'}, end => $spike->{'end'}};
		my $length_by2 = int(($min_length - $spike_length)/2);
		$new_spike->{'begin'} = max (0, ($new_spike->{'begin'} - $length_by2));
		$new_spike->{'end'} = min($chr_length, $new_spike->{'begin'} + $min_length);
		$new_spike->{'strand'} = $spike->{'strand'};
		#~ $new_spike->{read_count} = $spike->{read_count};
		#~ $new_spike->{forward_read_count} = $spike->{forward_read_count};
		return $new_spike;
	}
	return $spike;
}


=method __force_spike_size

Static private helper function. You shouldnt use this function.

=cut
sub __force_spike_size {
	my ($spike, $length, $chr_length) = @_;
	my $spike_length = $spike->{'end'} - $spike->{'begin'};
	if ($spike_length < $length) {
		my $new_spike = {begin => $spike->{'begin'}, end => $spike->{'end'}};
		my $length_by2 = int(($length - $spike_length)/2);
		$new_spike->{'begin'} = max (0, ($new_spike->{'begin'} - $length_by2));
		$new_spike->{'end'} = min ($chr_length, $new_spike->{'begin'} + $length);
		$new_spike->{'strand'} = $spike->{'strand'};
		#~ $new_spike->{read_count} = $spike->{read_count};
		#~ $new_spike->{forward_read_count} = $spike->{forward_read_count};
		return $new_spike;
	}
	elsif ($spike_length > $length) {
		my $new_spike = {begin => $spike->{'begin'}, end => $spike->{'end'}};
		my $length_by2 = int(($spike_length - $length)/2);
		$new_spike->{'begin'} = $new_spike->{'begin'} + $length_by2;
		$new_spike->{'end'} = $new_spike->{'begin'} + $length;
		$new_spike->{'strand'} = $spike->{'strand'};
		#~ $new_spike->{read_count} = $spike->{read_count};
		#~ $new_spike->{forward_read_count} = $spike->{forward_read_count};
		return $new_spike;
	}
	return $spike;
}


=method __look_forward

Static private helper function. You shouldnt use this function.

=cut
sub __look_forward {
	my ($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $enlarged_spike, $window_length) = @_;
	my $genome_seq_length = $genome_seq_end-$genome_seq_beg;
	my @candidate_region = (min ($enlarged_spike->{'end'}+15, $chr_length), min ($enlarged_spike->{'end'}+$window_length+15, $chr_length));
	if ($candidate_region[1] - $candidate_region[0] < $enlarged_spike->{'end'} - $enlarged_spike->{'begin'}) {
		# At the right end of the chromosome
		return 0;
	}

	if ($enlarged_spike->{'begin'} < $genome_seq_beg || $candidate_region[0] < $genome_seq_beg ||
		$enlarged_spike->{'end'} > $genome_seq_end || $candidate_region[1] > $genome_seq_end) {
		debug('Looking forward... Genome seq: '. display_chr_coordinatees($chr, $genome_seq_beg, $genome_seq_end) . "\n", 1);
		debug('Too long to test ' . display_chr_coordinatees($chr, $enlarged_spike->{'begin'}, $enlarged_spike->{'end'}) .
			' against ' . display_chr_coordinatees($chr, $candidate_region[0], $candidate_region[1]) . "\n", 1);
		return 0;
	}
	my $results;
	my $strand = $enlarged_spike->{'strand'};
	if ($enlarged_spike->{'strand'} eq '?') {
		$results = $detector->detect_on_strand('+', $genome_seq, $genome_seq_length, $enlarged_spike->{'begin'}-$genome_seq_beg, $enlarged_spike->{'end'}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
		$strand = '+';
		if ( scalar @{$results} ) {
			my @farthest = @{$results->[0]};
			foreach my $r ( @{$results} ) {
				@farthest = @{$r} if $r->[0] > $farthest[0]; # if begin_sequence of 'r' > begin_sequence of 'farthest'
			}
			my $other = {begin => $candidate_region[0]+$farthest[0], end => $candidate_region[0]+$farthest[1], strand => $enlarged_spike->{'strand'}};
			push @{$miRnaPos}, {strand => $strand, from_read => $enlarged_spike, detected => $other, first => $enlarged_spike, second => $other, source => DUE_TO_SINGLE_SPIKE};
		}
		$results = $detector->detect_on_strand('-', $genome_seq, $genome_seq_length, $enlarged_spike->{'begin'}-$genome_seq_beg, $enlarged_spike->{'end'}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
		$strand = '-';
	}
	else {
		$results = $detector->detect_on_strand($enlarged_spike->{'strand'}, $genome_seq, $genome_seq_length, $enlarged_spike->{'begin'}-$genome_seq_beg, $enlarged_spike->{'end'}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
	}
	# $results = [ [[begin_read, end_read], [begin_sequence, end_sequence]], [[begin_read, end_read], [begin_sequence, end_sequence]], ...]
	if ( scalar @{$results} ) {
		my @farthest = @{$results->[0]};
		foreach my $r ( @{$results} ) {
			@farthest = @{$r} if $r->[0] > $farthest[0]; # if begin_sequence of 'r' > begin_sequence of 'farthest'
		}
		my $other = {begin => $candidate_region[0]+$farthest[0], end => $candidate_region[0]+$farthest[1], strand => $enlarged_spike->{'strand'}};
		push @{$miRnaPos}, {strand => $strand, from_read => $enlarged_spike, detected => $other, first => $enlarged_spike, second => $other, source => DUE_TO_SINGLE_SPIKE};
		return 1;
	}
	return 0;
}


=method __look_backward

Static private helper function. You shouldnt use this function.

=cut
sub __look_backward {
	my ($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $enlarged_spike, $window_length) = @_;
	my $genome_seq_length = $genome_seq_end-$genome_seq_beg;
	my @candidate_region = (max ($enlarged_spike->{'begin'}-$window_length-15, 0), max ($enlarged_spike->{'begin'}-15, 0));
	if ($candidate_region[1] - $candidate_region[0] < $enlarged_spike->{'end'} - $enlarged_spike->{'begin'}) {
		return 0;
	}
	if ($enlarged_spike->{'begin'} < $genome_seq_beg || $candidate_region[0] < $genome_seq_beg ||
		$enlarged_spike->{'end'} > $genome_seq_end || $candidate_region[1] > $genome_seq_end) {
 		debug('Looking backward... Genome seq: '.display_chr_coordinatees($chr, $genome_seq_beg, $genome_seq_end) . "\n", 1);
		debug('Too long to test ' . display_chr_coordinatees($chr, $enlarged_spike->{'begin'}, $enlarged_spike->{'end'}) .
			' against ' . display_chr_coordinatees($chr, $candidate_region[0], $candidate_region[1]) . "\n", 1);
		return 0;
	}

	my $results;
	my $strand = $enlarged_spike->{'strand'};
	if ($enlarged_spike->{'strand'} eq '?') {
		$results = $detector->detect_on_strand('+', $genome_seq, $genome_seq_length, $enlarged_spike->{'begin'}-$genome_seq_beg, $enlarged_spike->{'end'}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
		$strand = '+';
		if (scalar @{$results}) {
			my @farthest = @{$results->[0]};
			foreach my $r ( @{$results} ) {
				@farthest = @{$r} if $r->[0] < $farthest[0]; # if begin_sequence of 'r' < begin_sequence of 'farthest'
			}
			my $other = {begin => $candidate_region[0]+$farthest[0], end => $candidate_region[0]+$farthest[1], strand => $enlarged_spike->{'strand'}};
			push @{$miRnaPos}, {strand => $strand, from_read => $enlarged_spike, detected => $other, first => $other, second => $enlarged_spike, source => DUE_TO_SINGLE_SPIKE};
		}
		$results = $detector->detect_on_strand('-', $genome_seq, $genome_seq_length, $enlarged_spike->{'begin'}-$genome_seq_beg, $enlarged_spike->{'end'}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
		$strand = '-';
	}
	else {
		$results = $detector->detect_on_strand($enlarged_spike->{'strand'}, $genome_seq, $genome_seq_length, $enlarged_spike->{'begin'}-$genome_seq_beg, $enlarged_spike->{'end'}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
	}
	return 0 if scalar @{$results} == 0;
	my @farthest = @{$results->[0]};
	foreach my $r ( @{$results} ) {
		@farthest = @{$r} if $r->[0] < $farthest[0]; # if begin_sequence of 'r' < begin_sequence of 'farthest'
	}
	my $other = {begin => $candidate_region[0]+$farthest[0], end => $candidate_region[0]+$farthest[1], strand => $enlarged_spike->{'strand'}};
	push @{$miRnaPos}, {strand => $strand, from_read => $enlarged_spike, detected => $other, first => $other, second => $enlarged_spike, source => DUE_TO_SINGLE_SPIKE};
	return 1;
}


=method __look_both_ways

Static private helper function. You shouldnt use this function.

=cut
sub __look_both_ways {
	my ($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $min_length, $window_length) = @_;
	my $result_back = __look_backward($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $window_length);
	my $result_forw = __look_forward($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $window_length);
	return $result_back || $result_forw;
}


=method __match_with_neighbor

Static private helper function. You shouldnt use this function.

=cut
sub __match_with_neighbor {
	my ($detector, $miRnaPos, $strand, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $neighbor, $min_length) = @_;
	my $genome_seq_length = $genome_seq_end-$genome_seq_beg;
	my $enlarged_neighbor = __enlarged_spike($neighbor, $min_length, $chr_length);
	my $spike_len = $current_spike->{'end'} - $current_spike->{'begin'};
	my $neighbor_len = $enlarged_neighbor->{'end'} - $enlarged_neighbor->{'begin'};
	my $results;
	if ($current_spike->{'begin'} < $genome_seq_beg || $enlarged_neighbor->{'begin'} < $genome_seq_beg ||
		$current_spike->{'end'} > $genome_seq_end || $enlarged_neighbor->{'end'} > $genome_seq_end) {
		debug('Matching neighbors... Genome seq: '.display_chr_coordinatees($chr, $genome_seq_beg, $genome_seq_end) . "\n", 1);
		debug('Too long to test ' . display_chr_coordinatees($chr, $current_spike->{'begin'}, $current_spike->{'end'}) .
			' against ' . display_chr_coordinatees($chr, $enlarged_neighbor->{'begin'}, $enlarged_neighbor->{'end'}) . "\n", 1);
		return 0;
	}

	if ($spike_len < $neighbor_len) {
		$results = $detector->detect_on_strand($strand, $genome_seq, $genome_seq_length, $current_spike->{'begin'}-$genome_seq_beg,
		$current_spike->{'end'}-$genome_seq_beg, $enlarged_neighbor->{'begin'}-$genome_seq_beg, $enlarged_neighbor->{'end'}-$genome_seq_beg);
		if ( scalar @{$results} ) {
			push @{$miRnaPos}, {strand => $strand, first => $current_spike, second => __force_spike_size($neighbor, $spike_len, $chr_length), source => DUE_TO_TWO_SPIKES};
			return 1;
		}
	}
	else {
		$results = $detector->detect_on_strand($strand, $genome_seq, $genome_seq_length, $enlarged_neighbor->{'begin'}-$genome_seq_beg, $enlarged_neighbor->{'end'}-$genome_seq_beg,
		$current_spike->{'begin'}-$genome_seq_beg, $current_spike->{'end'}-$genome_seq_beg);
		if ( scalar @{$results} ) {
			push @{$miRnaPos}, {strand => $strand, first => $current_spike, second => __force_spike_size($neighbor, $spike_len, $chr_length), source => DUE_TO_TWO_SPIKES};
			return 1;
		}
	}
	return 0;
}


=method process_window_spikes_for_chr

Private helper function. You shouldnt use this function.

=cut
sub process_spikes_for_chr {
	my ($this, $chr, $spikes_for_chr
	#~ STAT BEG
	#~ , $stats
	#~ STAT END
	) = @_;
	$this->{detector} = MiRnaDuplexDetector::MiRnaDetector->new(5000);
	my $accepting_time = $this->{accepting_time};
	my @miRnaPos = (); # Array of [5p, 3p]
	my $genome = $this->{'genome_db'};
	my $detector = $this->{detector};
	my $min_length = 21;
	my $min_length_for_neighbor = 40;
	my $chr_length = $this->{chr_info}{$chr};
	my $min_length_on_isolated_spike = 21;
	my $window_length = 400;
	my $enlarging = max($window_length, $accepting_time);
	#~ my $min_error_for_single_spike = $detector->miRnaHigherScoreThreshold();
	#~ my $min_error_for_two_spikes = $detector->miRnaLowerScoreThreshold();

	foreach my $current_spike_entry (@{$spikes_for_chr}) {
		my $current_spike = $current_spike_entry->{'spike'};
		my $unknown_strand = $current_spike->{'strand'} eq '?';
		my $enlarged_spike = __force_spike_size($current_spike, $min_length, $chr_length);

		my $genome_seq_beg = max($enlarged_spike->{'begin'}-50, 0);
		my $genome_seq_end = min($enlarged_spike->{'end'}+$accepting_time*2+100, $chr_length);
		if ($genome_seq_end-$genome_seq_beg > $detector->admissibleTextLength()) {
			$genome_seq_end = $genome_seq_beg+$detector->admissibleTextLength();
		}
		my $genome_seq = $this->get_sub_sequence($chr, $genome_seq_beg, $genome_seq_end);
		my $has_results = 0;

		#~ STATS BEG
		#~ my $stats__couple_spike = $current_spike_entry->{'class'} != ISOLATED_SPIKE ? 1 : 0;
		#~ STATS END
		if ($current_spike_entry->{'class'} != ISOLATED_SPIKE) {
			#~ $detector->setMiRnaMinErrorsThreshold($min_error_for_two_spikes);
			foreach my $neighbor (@{$current_spike_entry->{'neighbors'}}) {
				if ($unknown_strand || $neighbor->{'strand'} eq '?') {
					$has_results = __match_with_neighbor($detector, \@miRnaPos, '+', $genome_seq_beg, $genome_seq_end, $genome_seq,
									$chr, $chr_length, $enlarged_spike, $neighbor, $min_length_for_neighbor);
					my $has_results_2 = __match_with_neighbor($detector, \@miRnaPos, '-', $genome_seq_beg, $genome_seq_end, $genome_seq,
								$chr, $chr_length, $enlarged_spike, $neighbor, $min_length_for_neighbor);
					$has_results = $has_results || $has_results_2;
				}
				else {
					$has_results = __match_with_neighbor($detector, \@miRnaPos, $enlarged_spike->{'strand'}, $genome_seq_beg, $genome_seq_end, $genome_seq,
									$chr, $chr_length, $enlarged_spike, $neighbor, $min_length_for_neighbor);
				}
			}
		}
		if ($has_results == 0) {
			#~ $detector->setMiRnaMinErrorsThreshold($min_error_for_single_spike);

			$genome_seq_beg = max($enlarged_spike->{'begin'}-$window_length-15, 0);
			$genome_seq_end = min($enlarged_spike->{'end'}+$window_length+15, $chr_length);
			$genome_seq = $this->get_sub_sequence($chr, $genome_seq_beg, $genome_seq_end);

			$has_results = __look_both_ways($detector, \@miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $enlarged_spike,
			$min_length_on_isolated_spike, $window_length);

			#~ STATS BEG
			#~ if ($stats__couple_spike == 1 && $has_results != 0) {
				#~ $stats->{'COUPLE_OF_SPIKES_TO_ISOLATED_SPIKE'}++;
			#~ }
			#~ STATS END
		}

		#~ STATS BEG
		#~ if ($has_results == 0) {
			#~ if ($stats__couple_spike == 1) {
				#~ $stats->{'COUPLE_OF_SPIKES_DISCARDED'}++;
			#~ }
			#~ else {
				#~ $stats->{'ISOLATED_SPIKE_DISCARDED'}++;
			#~ }
		#~ }
		#~ STATS END
	}
	#~ $detector->setMiRnaMinErrorsThreshold($min_error_for_single_spike);
    miRkwood::Utils::display_var_sizes_in_log_file( '..... ClusterJobSebastien : process_spikes_for_chr' );
	return \@miRnaPos;
}


=method region_intertects

Static private helper function. You shouldnt use this function.

=cut
sub region_intertects {
	my ($a, $b, $ratio) = @_;
	if ($a->{'begin'} <= $b->{'begin'} && $a->{'end'} > $b->{'begin'}) {
		if ($b->{'end'} <= $a->{'end'}) {
			return 1;
		}
		return ($a->{'end'}-$b->{'begin'})/($b->{'end'}-$a->{'begin'}) >= $ratio ? 1 : 0;
	}
	elsif ($b->{'begin'} <= $a->{'begin'} && $b->{'end'} > $a->{'begin'}) {
		if ($a->{'end'} <= $b->{'end'}) {
			return 1;
		}
		return ($b->{'end'}-$a->{'begin'})/($a->{'end'}-$b->{'begin'}) >= $ratio ? 1 : 0;
	}
	return 0;
}


=method raw_regions_intertect

Static private helper function. You shouldnt use this function.

=cut
sub raw_regions_intertect {
	my ($a, $b) = @_;
	if ($a->{'begin'} <= $b->{'begin'} && $a->{'end'} > $b->{'begin'}) {
		return 1;
	}
	elsif ($b->{'begin'} <= $a->{'begin'} && $b->{'end'} > $a->{'begin'}) {
		return 1;
	}
	return 0;
}


=method merge_regions

Static private helper function. You shouldnt use this function.

=cut
# $a and $b must intersect
# the result is merged in $a
sub merge_regions {
	my ($a, $b) = @_;
	if ($a->{'begin'} <= $b->{'begin'}) {
		if ($b->{'end'} > $a->{'end'}) {
			$a->{'end'} = $b->{'end'};
		}
	}
	else {
		if ($b->{'end'} <= $a->{'end'}) {
			$a->{'begin'} = $b->{'begin'};
		}
		else {
			$a->{'begin'} = $b->{'begin'};
			$a->{'end'} = $b->{'end'};
		}
	}
	push @{$a->{'miRnas'}}, @{$b->{'miRnas'}};
    return;
}

sub merge_overlapping_regions {
	my $regions = shift;
	my @sorted_regions = sort {$a->{'begin'} <=> $b->{'begin'}} @{$regions};

	if (scalar @sorted_regions == 0) {
		return [];
	}

	my @stack = ($sorted_regions[0]);
	for (my $i = 1; $i < scalar @sorted_regions; $i++) {
		my $top = $stack[0];
		if ($top->{'end'} < $sorted_regions[$i]{'begin'}) {
			unshift @stack, $sorted_regions[$i];
		}
		elsif ($top->{'end'} < $sorted_regions[$i]{'end'}) {
			$top->{'end'} = $sorted_regions[$i]{'end'};
			push @{$top->{'miRnas'}}, @{$sorted_regions[$i]{'miRnas'}};
			shift @stack;
			unshift @stack, $top;
		}
	}

	return \@stack;
}


=method compute_candidate_precursors_from_miRnaPos_for_chr

Static private helper function. You shouldnt use this function.

Computes candidate precursors based miRNA candidates. This is made by adding 100 nt on both sides of the miRna couple.
Regions that overlap by more than 60% are merged together.

 Usage : my $candidate_precursors = $self->compute_candidate_precursors_from_miRnaPos($miRnas);
 Input : The miRNAs returned by process_window_spikes
 Return: A hash ref {
			chr => [array of precursors {
				begin => start position (1-based)
				end => end position (excluded)
				strand => '+' or '-'
				miRnas => [] Array reference of miRNA candidates contained in the precursor. Picked from $miRnas.
				}]
			}

=cut
sub compute_candidate_precursors_from_miRnaPos_for_chr {
	my $this = shift;
	my $chr = shift;
	my $miRnaPos = shift;
	my $chr_length = shift;
    my $average_coverage = shift;
	my $read_coverage_threshold = shift;
	my $peak_padding = shift;
	my $parsed_bed = shift;

	if (scalar @{$miRnaPos} == 0) {
		return [];
	}
	my %regions = ('+' => [], '-' => []);

	for (my $i = 0; $i < scalar @{$miRnaPos}; $i++) {
		my $current_miRna = $miRnaPos->[$i];
		my $region_begin = max(0, $current_miRna->{first}{'begin'}-$peak_padding);
		my $region_end = min($chr_length, $current_miRna->{second}{'end'}+$peak_padding);
		my $threshold = ( $region_end - $region_begin - 19 ) / $average_coverage;
		my $final_threshold = max( $read_coverage_threshold, $threshold );
		if (miRkwood::HairpinBuilder::get_contained_read_coverage($parsed_bed, $chr, $region_begin, $region_end, $current_miRna->{'strand'})
		< $final_threshold) {
			next;
		}
		my %region = (
			'begin'  => $region_begin,
			'end'    => $region_end,
			'strand' => $current_miRna->{'strand'},
			'miRnas' => [$current_miRna],
			'chr' => $chr);

		push @{$regions{$current_miRna->{'strand'}}}, \%region;
	}

	$regions{'+'} = merge_overlapping_regions($regions{'+'});
	$regions{'-'} = merge_overlapping_regions($regions{'-'});

	my @final_regions = sort {$a->{'begin'} <=> $b->{'begin'}} (@{$regions{'+'}}, @{$regions{'-'}});
	return \@final_regions;
}

1;
