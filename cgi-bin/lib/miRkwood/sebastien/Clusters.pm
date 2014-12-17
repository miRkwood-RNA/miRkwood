# ABSTRACT: Handling the cluster making from a BAM file.

package Clusters;

use strict;
use warnings;
use POSIX;
use MiRnaDuplexDetector;
use Bio::DB::Fasta;
use KMean;

use List::Util qw(max min);

use constant {
	DUE_TO_SINGLE_SPIKE => 1,
	DUE_TO_TWO_SPIKES => 2,
};

use lib '/home/sebastien/mirkwood/cgi-bin/lib';
use miRkwood::Programs;


=method new

Constructor

 Usage : my $clustering = Clusters->new('genome.fa'),
 Input : The genome file
 Return: The cunstructed instance

=cut

sub new {
    my ( $class, @args ) = @_;
    my $genome_file = shift @args;
    my $genome_db = Bio::DB::Fasta->new($genome_file);
    my $self = bless {
        genome_file => $genome_file,
        genome_db => $genome_db,
        accepting_time => 350
    }, $class;
    my %chr_info = $self->get_chromosomes_info_from_genome_file();
    $self->{chr_info} = \%chr_info;
    return $self;
}

=method get_read_distribution_from_bam

Retrieve the reads from a bam file.

 Usage : my $reads = $self->get_read_distribution_from_bam('reads.bam'),
 Input : The bam file
 Return: A hash reference {chr => {
					begin_pos => {
						read_count => read depth,
						end => maximum end coordinates of all reads starting at begin_pos,
						forward_read_count => read depth of reads mapped onto forward strand
						}
					}
				}

=cut
sub get_read_distribution_from_bam {
    my ( $this, $bam_file ) = @_;
    # go chr by chr, using the .fai index file to get the chr names
    my %reads = ();
    my @chrs = keys %{ $this->{chr_info} };

    foreach my $chr (@chrs) {
        my $samtools_cmd = "samtools view $bam_file $chr";
        open( my $DEPTH, '-|', $samtools_cmd ) or die "Can't open '$samtools_cmd': $!";
        $reads{$chr} = __get_read_distribution_from_bam_for_chr($DEPTH);
        close $DEPTH;
    }
    return \%reads;
}

=method __get_read_distribution_from_bam_for_chr

Static private helper function. You shouldnt use this function.

 Usage : my $reads = __get_read_distribution_from_bam_for_chr($sam_file_handle),
 Input : The bam file path
 Return: A hash reference {begin_pos => {
					read_count => read depth,
					end => maximum end coordinates of all reads starting at begin_pos,
					forward_read_count => read depth of reads mapped onto forward strand
					}
				}

=cut
sub __get_read_distribution_from_bam_for_chr {
    my $HANDLE = shift;
    my %chr_reads = ();
    while (<$HANDLE>) {
        chomp;
        my @fields = split("\t");
        my $pos = $fields[3];
        my $strand = '+';
        $strand = '-' if $fields[1] & 0x10;
        my $end = $fields[3] + length $fields[9];
        if (defined $chr_reads{$pos}) {
			$chr_reads{$pos}->{read_count}++;
			$chr_reads{$pos}->{end} = $end if $end > $chr_reads{$pos}->{end};
			$chr_reads{$pos}->{forward_read_count}++ if $strand eq '+';
        }
        else {
			$chr_reads{$pos} = {read_count => 1, end => $end, forward_read_count => ($strand eq '+') ? 1 : 0};
        }
    }
    return \%chr_reads;
}


=method get_read_distribution_from_bed

Retrieve the reads from a miRkwood-normalized bed file. The bed file must have the following columns (tab-separated):
# 0: chromosom name,
# 1: starting position (0-based, starting at 0),
# 2: end positing (excluded)
# 3: sequence name
# 4: depth (integer)
# 5: strand ('+' or '-')

 Usage : my $reads = $self->get_read_distribution_from_bed('reads.bed'),
 Input : The bed file path
 Return: A hash reference {chr => {
					begin_pos => {
						read_count => read depth,
						end => maximum end coordinates of all reads starting at begin_pos,
						forward_read_count => read depth of reads mapped onto forward strand
						}
					}
				}

=cut
sub get_read_distribution_from_bed {
    my ($this, $bed_file) = @_;
    my %reads = ();
    foreach my $chr (keys %{$this->{chr_info}}) {
		$reads{$chr} = {};
    }
    open( my $HANDLE, '<', $bed_file) or die "Can't open '$bed_file': $!";
    while (<$HANDLE>) {
        chomp;
        my @fields = split("\t");
        if (scalar @fields != 6) {
			next;
        }
        my $pos = $fields[1]+1;
        my $strand = $fields[5];
        my $end = $fields[2]+1;
        my $chr = $fields[0];
        if (defined $reads{$chr}{$pos}) {
			$reads{$chr}{$pos}{read_count} += $fields[4];
			$reads{$chr}{$pos}{end} = $end if $end > $reads{$chr}{$pos}{end};
			$reads{$chr}{$pos}{forward_read_count}+= $fields[4] if $strand eq '+';
        }
        else {
			$reads{$chr}{$pos} = {read_count => $fields[4], end => $end, forward_read_count => ($strand eq '+') ? $fields[4] : 0};
        }
    }
    close $HANDLE;
    return \%reads;
}

=method get_faidx_file


=cut
sub get_faidx_file {
    my ( $self, @args ) = @_;
    my $expected_faidx = $self->{genome_file} . '.fai';
    if ( !-e $expected_faidx ) {
        my $samtools_cmd = "samtools faidx $self->{genome_file}";
        system $samtools_cmd;
    }
    return $expected_faidx;
}


=method get_chromosomes_info_from_genome_file

Retrieve chromosomes name and length and from FAI file.

 Usage : my %chr_info_output = $self->get_chromosomes_info_from_genome_file(),
 Input : The genome file
 Return: A hash {name => length}

=cut

sub get_chromosomes_info_from_genome_file {
    my ($self, @args) = @_;
    my $genome_file = $self->{genome_file};
    my $fai_file = $self->get_faidx_file();
    my %chr_lengths;
    open( my $FAI, '<', "$fai_file" )
      or die "Error when opening $fai_file: $!";
    while (<$FAI>) {
        chomp;
        my @fields = split( "\t", $_ );
        my $chr_name = $fields[0];
        my $chr_length = $fields[1];
        $chr_lengths{$chr_name} = $chr_length;
    }
    return %chr_lengths;
}


=method get_windows

Retrieve regions of interest by reading through the reads. The threshold defines the minimum read depth for the read to be caught in a window.
A window is made of read trains (a read train is a region of overlapping reads). A window gather a cluster of read trains (no two consecutive read trains are separated by more than
$self->{accepting_time} (= 350 nt by default)).

 Usage : my $windows = $self->get_windows($read_distribution, $detection_threshold),
 Input : The read distribution returned from get_read_distribution_from_bed or get_read_distribution_from_bam
 Return: A hash reference {
				chr => [array of hash refs {
							begin => start position (1-based)
							end => end position (excluded)
							read_count => The number of reads contained in the window (considering the read depth)
							forward_read_count => The number of reads contained in the window mapped onto the forward strand (considering the read depth)
							trains => [array of hash refs {
								begin => start position (1-based)
								end => end position (excluded)
								last_read_begin => start position of the last read of the train (1-based)
								read_count => The number of reads contained in the window (considering the read depth)
								forward_read_count => The number of reads contained in the window mapped onto the forward strand (considering the read depth)
								spikes => [array of hash refs {
									begin => start position (1-based)
									end => end position (excluded)
									read_count => The number of reads contained in the window (considering the read depth)
									forward_read_count => The number of reads contained in the window mapped onto the forward strand (considering the read depth)
									strand => '+' or '-' or '?'
									trigger => the local read coverage that triggered the spike
									}]
								}]
							}
						]
				}

=cut
sub get_windows {
	my $this = shift;
	my $read_distribution_per_chr = shift;
	my $train_detection_threshold = shift;
	my %windows_per_chr = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
		print "\t", $chr, "\n";
		$windows_per_chr{$chr} = $this->__get_windows_for_chr($read_distribution_per_chr->{$chr}, $train_detection_threshold);
	}
	return \%windows_per_chr;
}


=method get_windows

Static private helper function that returns the strand of a region based on read statistics.

 Usage : my $strand = get_strand($forward_read_count, $read_count),
 Input : The number of reads mapped onto the forward strand, and mapped on both strand
 Return: '+' (if more than 70% of reads are on the forward strand) or '-' (if more than 70% of reads are on reverse strand) or '?' (otherwise)

=cut
sub get_strand {
	my ($forward_read_count, $read_count) = @_;
	return $forward_read_count >= $read_count*.7 ? '+' : $forward_read_count <= $read_count*.3 ? '-' : '?';
}


=method __get_windows_for_chr

Private helper function that returns the windows for a given chr. You shouldnt use this function.

=cut
sub __get_windows_for_chr {
	my $this = shift;
	my $read_distribution = shift;
	my $train_detection_threshold = shift;

	my @positions = sort {$a <=> $b} keys %$read_distribution;

	my @windows = ();
	if (scalar(@positions) == 0) {
		return \@windows;
	}

	my %current_train = (begin => $positions[0], end => $read_distribution->{$positions[0]}{end}, read_count => 0, forward_read_count => 0,
	spikes => [], classifier => KMean->new());
	$current_train{classifier}->add_point(0);
	my %window = (begin => 0, read_count => 0, forward_read_count => 0, end => 0, trains => []);
	my $last_read_count = 0;
	my $spike_detection = 2;
	my $total_read_count = 0;
	my @end_reads = ();

	#DEGUB
	my $seeking_pos = 0;
	# END DEBUG

	foreach my $position (@positions) {
		my $read_locus = $read_distribution->{$position};
		if ($current_train{begin} <= $position && $position < $current_train{end}) {
			$current_train{end} = max $current_train{end}, $read_locus->{end};
			$current_train{begin_offset_count} += $position - $current_train{begin};
			$current_train{length_count} += $read_locus->{end} - $current_train{begin};
			$current_train{read_count} += $read_locus->{read_count};
			$current_train{forward_read_count} += $read_locus->{forward_read_count};
			$current_train{last_read_begin} = $position;
			__get_windows_process_train_spikes(\%current_train, $position, $read_locus, \$total_read_count, \@end_reads, \$last_read_count, $spike_detection, $seeking_pos);
		}
		else { # the read train ended
			__get_windows_maintain_read_count(\%current_train, $position, \$total_read_count, \$last_read_count, \@end_reads, $spike_detection);
			__get_windows_finish_train_spikes(\%current_train);
			$this->__get_windows_process_train_from_distribution(\@windows, \%current_train, \%window, $train_detection_threshold);
			%current_train = (begin => $position, end => $read_locus->{end}, read_count => $read_locus->{read_count}, forward_read_count => $read_locus->{forward_read_count}, spikes => [],
			classifier => $current_train{classifier});
			$total_read_count = 0;
			@end_reads = ();
			$last_read_count = 0;
			__get_windows_process_train_spikes(\%current_train, $position, $read_locus, \$total_read_count, \@end_reads, \$last_read_count, $spike_detection, $seeking_pos);
		}
		$last_read_count = $total_read_count;
	}
	__get_windows_finish_train_spikes(\%current_train);
	$this->__get_windows_process_train_from_distribution(\@windows, \%current_train, \%window, $train_detection_threshold);
	if ($window{begin} != $window{end}) {
		$this->__add_candidate_window_from_train(\@windows, \%window);
	}
	return \@windows;
}


=method __get_windows_maintain_read_count

Static private helper function. You shouldnt use this function.

=cut
sub __get_windows_maintain_read_count {
	my ($current_train, $position, $total_read_count, $last_read_count, $end_reads, $spike_detection) = @_;
	my $spikes = $current_train->{spikes};

	# Maintaining read count
	if (scalar @$spikes && $spikes->[-1]{end} == -1) {
		my $last_spike = $spikes->[-1];
		while (scalar @$end_reads && $end_reads->[0]{end} <= $position) {
			$$last_read_count = $$total_read_count;
			$$total_read_count -= $end_reads->[0]{read_count};
			my $trigger = $end_reads->[0]{end};
			shift @$end_reads;
			my $added = $current_train->{classifier}->add_point($$total_read_count);
			if ($added == KMean::ASSIGNED_FIRST_CLASS) {
				$last_spike->{end} = $trigger;
				$last_spike->{strand} = get_strand $last_spike->{forward_read_count}, $last_spike->{read_count};
				$current_train->{classifier}->clear_points();
				$current_train->{classifier}->add_point($$total_read_count);
				last;
			}
		}
	}
	while (scalar @$end_reads && $end_reads->[0]{end} <= $position) {
		$$last_read_count = $$total_read_count;
		$$total_read_count -= $end_reads->[0]{read_count};
		shift @$end_reads;
	}
}


=method __get_windows_finish_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __get_windows_finish_train_spikes {
	my ($current_train) = @_;
	my $spikes = $current_train->{spikes};
	if (scalar @$spikes && $spikes->[-1]{end} == -1) {
		my $last_spike = $spikes->[-1];
		$last_spike->{end} = $current_train->{end};
		$last_spike->{strand} = get_strand $last_spike->{forward_read_count}, $last_spike->{read_count};
	}
	elsif (scalar @$spikes == 0 && $current_train->{end} - $current_train->{begin} <= 40) {
		my $strand = get_strand $current_train->{forward_read_count}, $current_train->{read_count};
		push @$spikes, {begin => $current_train->{begin}, end => $current_train->{end}, trigger => 0,
		read_count => $current_train->{read_count}, forward_read_count => $current_train->{forward_read_count}, strand => $strand};
	}
	$current_train->{classifier}->clear_points();
	$current_train->{classifier}->add_point(0);
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __get_windows_process_train_spikes {
	my ($current_train, $position, $read_locus, $total_read_count, $end_reads, $last_read_count, $spike_detection, $seeking_pos) = @_;
	my $spikes = $current_train->{spikes};

	# Maintaining read count
	__get_windows_maintain_read_count($current_train, $position, $total_read_count, $last_read_count, $end_reads, $spike_detection);
	$$total_read_count += $read_locus->{read_count};
	push @$end_reads, {end => $read_locus->{end}, read_count => $read_locus->{read_count}};
	@$end_reads = sort {$a->{end} <=> $b->{end}} @$end_reads;
	# Looking for spikes
	my $added = $current_train->{classifier}->add_point($$total_read_count);
	if ($added == KMean::ASSIGNED_SECOND_CLASS) {
		if (scalar @$spikes) {
			my $last_spike = $spikes->[-1];
			if ($last_spike->{end} == -1) {
				if ($current_train->{classifier}->class_of($last_spike->{trigger}) == KMean::ASSIGNED_FIRST_CLASS) {
					$last_spike->{begin} = $position;
					$last_spike->{trigger} = $$total_read_count;
					$last_spike->{read_count} = $read_locus->{read_count};
					$last_spike->{forward_read_count} = $read_locus->{forward_read_count};
				}
				return;
			}
		}
		push @$spikes, {begin => $position, end => -1, trigger => $$total_read_count, read_count => $read_locus->{read_count},
		forward_read_count => $read_locus->{forward_read_count}};
	}
	elsif ($added == KMean::ASSIGNED_FIRST_CLASS) {
		my $last_spike = $spikes->[-1];
		if ($last_spike->{end} == -1) {
			$last_spike->{end} = $position;
			$last_spike->{strand} = get_strand $last_spike->{forward_read_count}, $last_spike->{read_count};
		}
		$current_train->{classifier}->clear_points();
		$current_train->{classifier}->add_point($$total_read_count);
	}
	elsif (scalar @$spikes && $spikes->[-1]{end} == -1) {
		my $last_spike = $spikes->[-1];
		$last_spike->{read_count} += $read_locus->{read_count};
		$last_spike->{forward_read_count} += $read_locus->{forward_read_count};
	}
}


=method __get_windows_process_train_spikes

Private helper function. You shouldnt use this function.

=cut
sub __get_windows_process_train_from_distribution {
	my ($this, $windows, $current_train, $window, $train_detection_threshold) = @_;
	my $accepting_time = $this->{accepting_time};
	my $reverse_read_count = $current_train->{read_count} - $current_train->{forward_read_count};
	if ($current_train->{forward_read_count} >= $train_detection_threshold || $reverse_read_count >= $train_detection_threshold) { # This train should be in a window
		if ($window->{begin} == $window->{end}) { # There were no current window
			$window->{begin} = $current_train->{begin};
			$window->{end} = $current_train->{end};
			$window->{read_count} = $current_train->{read_count};
			$window->{forward_read_count} = $current_train->{forward_read_count};
			$window->{trains} = [{begin => $current_train->{begin}, end => $current_train->{end}, last_read_begin => $current_train->{last_read_begin},
			read_count => $current_train->{read_count}, forward_read_count => $current_train->{forward_read_count}, spikes => $current_train->{spikes}}];
		}
		elsif ($current_train->{begin} - $window->{end} < $accepting_time) { # We accept this train as it his not far from the last window
			$window->{end} = $current_train->{end};
			$window->{read_count} += $current_train->{read_count};
			$window->{forward_read_count} += $current_train->{forward_read_count};
			push @{$window->{trains}}, {begin => $current_train->{begin}, end => $current_train->{end}, last_read_begin => $current_train->{last_read_begin},
			read_count => $current_train->{read_count}, forward_read_count => $current_train->{forward_read_count}, spikes => $current_train->{spikes}};
		}
		else { # This train should go on a separate window
			$this->__add_candidate_window_from_train($windows, $window);
			$window->{begin} = $current_train->{begin};
			$window->{end} = $current_train->{end};
			$window->{read_count} = $current_train->{read_count};
			$window->{forward_read_count} = $current_train->{forward_read_count};
			$window->{trains} = [{begin => $current_train->{begin}, end => $current_train->{end}, last_read_begin => $current_train->{last_read_begin},
			read_count => $current_train->{read_count}, forward_read_count => $current_train->{forward_read_count}, spikes => $current_train->{spikes}}];
		}
	}
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
sub process_window_spikes {
	my $this = shift;
	my $windows_per_chr = shift;
	my %miRnaPos = ();

	$this->{detector} = MiRnaDuplexDetector::MiRnaDetector->new(5000);
	
	foreach my $chr (keys %{ $this->{chr_info} }) {
		print "\t", $chr, "\n";
		$miRnaPos{$chr} = $this->process_window_spikes_for_chr($chr, $windows_per_chr->{$chr});
	}
	return \%miRnaPos;
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __enlarged_spike {
	my ($spike, $min_length, $chr_length) = @_;
	my $spike_length = $spike->{end} - $spike->{begin};
	if ($spike_length < $min_length) {
		my $new_spike = {begin => $spike->{begin}, end => $spike->{end}};
		my $length_by2 = int(($min_length - $spike_length)/2);
		$new_spike->{begin} = max (1, ($new_spike->{begin} - $length_by2));
		$new_spike->{end} = $new_spike->{begin} + $min_length;#min ($chr_length, ($new_spike->{end} + $min_length - $spike_length - $length_by2));
		$new_spike->{strand} = $spike->{strand};
		return $new_spike;
	}
	return $spike;
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __shrink_spike {
	my ($spike, $max_length) = @_;
	my $spike_length = $spike->{end} - $spike->{begin};
	if ($spike_length > $max_length) {
		my $new_spike = {begin => $spike->{begin}, end => $spike->{end}};
		my $length_by2 = int(($max_length - $spike_length)/2);
		$new_spike->{begin} = $new_spike->{begin} + $length_by2;
		$new_spike->{end} = $new_spike->{begin} + $max_length;
		$new_spike->{strand} = $spike->{strand};
		return $new_spike;
	}
	return $spike;
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __force_spike_size {
	my ($spike, $length, $chr_length) = @_;
	my $spike_length = $spike->{end} - $spike->{begin};
	if ($spike_length < $length) {
		my $new_spike = {begin => $spike->{begin}, end => $spike->{end}};
		my $length_by2 = int(($length - $spike_length)/2);
		$new_spike->{begin} = max (1, ($new_spike->{begin} - $length_by2));
		$new_spike->{end} = $new_spike->{begin} + $length;#min ($chr_length, ($new_spike->{end} + $min_length - $spike_length - $length_by2));
		$new_spike->{strand} = $spike->{strand};
		return $new_spike;
	}
	elsif ($spike_length > $length) {
		my $new_spike = {begin => $spike->{begin}, end => $spike->{end}};
		my $length_by2 = int(($length - $spike_length)/2);
		$new_spike->{begin} = $new_spike->{begin} + $length_by2;
		$new_spike->{end} = $new_spike->{begin} + $length;
		$new_spike->{strand} = $spike->{strand};
		return $new_spike;
	}
	return $spike;
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __look_backward_only {
	my ($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $min_length, $window_length) = @_; # min_length = 40, window_length = 300
# 	my $enlarged_spike = __force_spike_size($current_spike, $min_length, $chr_length);
# 	my $mismatches_on_isolated_spike = int($mismatches*($enlarged_spike->{end} - $enlarged_spike->{begin})/$min_length);
	return __look_backward($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $window_length);
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __look_forward {
	my ($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $enlarged_spike, $window_length) = @_;
	my @candidate_region = (min ($enlarged_spike->{end}+15, $chr_length), min ($enlarged_spike->{end}+$window_length, $chr_length));
	if ($candidate_region[1] - $candidate_region[0] < $enlarged_spike->{end} - $enlarged_spike->{begin}) {
		# At the right end of the chromosome
		return 0;
	}

	if ($enlarged_spike->{begin} < $genome_seq_beg || $candidate_region[0] < $genome_seq_beg ||
		$enlarged_spike->{end} > $genome_seq_end || $candidate_region[1] > $genome_seq_end) {
		print "Looking forward... Genome seq: $chr:$genome_seq_beg-", $genome_seq_end-1, "\n";
		print "Too long to test $chr:$enlarged_spike->{begin}-", $enlarged_spike->{end}-1, " against $chr:$candidate_region[0]-", $candidate_region[1]-1, "\n";
		return 0;
	}
	my $results;
	my $strand = $enlarged_spike->{strand};
	if ($enlarged_spike->{strand} eq '?') {
		$results = $detector->detect_on_strand('+', $genome_seq, $enlarged_spike->{begin}-$genome_seq_beg, $enlarged_spike->{end}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
		$strand = '+';
		if (scalar @$results) {
			my @farthest = @{$results->[0]};
			foreach my $r (@$results) {
				@farthest = @$r if $r->[0] > $farthest[0]; # if begin_sequence of 'r' > begin_sequence of 'farthest'
			}
			my $other = {begin => $candidate_region[0]+$farthest[0], end => $candidate_region[0]+$farthest[1], strand => $enlarged_spike->{strand}};
			push @$miRnaPos, {strand => $strand, from_read => $enlarged_spike, detected => $other, first => $enlarged_spike, second => $other, source => DUE_TO_SINGLE_SPIKE};
		}
		$results = $detector->detect_on_strand('-', $genome_seq, $enlarged_spike->{begin}-$genome_seq_beg, $enlarged_spike->{end}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
		$strand = '-';
	}
	else {
		$results = $detector->detect_on_strand($enlarged_spike->{strand}, $genome_seq, $enlarged_spike->{begin}-$genome_seq_beg, $enlarged_spike->{end}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
	}
	# $results = [ [[begin_read, end_read], [begin_sequence, end_sequence]], [[begin_read, end_read], [begin_sequence, end_sequence]], ...]
	if (scalar @$results) {
		my @farthest = @{$results->[0]};
		foreach my $r (@$results) {
			@farthest = @$r if $r->[0] > $farthest[0]; # if begin_sequence of 'r' > begin_sequence of 'farthest'
		}
		my $other = {begin => $candidate_region[0]+$farthest[0], end => $candidate_region[0]+$farthest[1], strand => $enlarged_spike->{strand}};
		push @$miRnaPos, {strand => $strand, from_read => $enlarged_spike, detected => $other, first => $enlarged_spike, second => $other, source => DUE_TO_SINGLE_SPIKE};
		return 1;
	}
	return 0;
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __look_backward {
	my ($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $enlarged_spike, $window_length) = @_;
	my @candidate_region = (max ($enlarged_spike->{begin}-$window_length-15, 1), max ($enlarged_spike->{begin}-16, 1));
	if ($candidate_region[1] - $candidate_region[0] < $enlarged_spike->{end} - $enlarged_spike->{begin}) {
		return 0;
	}
	if ($enlarged_spike->{begin} < $genome_seq_beg || $candidate_region[0] < $genome_seq_beg ||
		$enlarged_spike->{end} > $genome_seq_end || $candidate_region[1] > $genome_seq_end) {
		print "Looking backward... Genome seq: $chr:$genome_seq_beg-", $genome_seq_end-1, "\n";
		print "Too long to test $chr:$enlarged_spike->{begin}-", $enlarged_spike->{end}-1, " against $chr:$candidate_region[0]-", $candidate_region[1]-1, "\n";
		return 0;
	}
	
	my $results;
	my $strand = $enlarged_spike->{strand};
	if ($enlarged_spike->{strand} eq '?') {
		$results = $detector->detect_on_strand('+', $genome_seq, $enlarged_spike->{begin}-$genome_seq_beg, $enlarged_spike->{end}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
		$strand = '+';
		if (scalar @$results) {
			my @farthest = @{$results->[0]};
			foreach my $r (@$results) {
				@farthest = @$r if $r->[0] < $farthest[0]; # if begin_sequence of 'r' < begin_sequence of 'farthest'
			}
			my $other = {begin => $candidate_region[0]+$farthest[0], end => $candidate_region[0]+$farthest[1], strand => $enlarged_spike->{strand}};
			push @$miRnaPos, {strand => $strand, from_read => $enlarged_spike, detected => $other, first => $other, second => $enlarged_spike, source => DUE_TO_SINGLE_SPIKE};
		}
		$results = $detector->detect_on_strand('-', $genome_seq, $enlarged_spike->{begin}-$genome_seq_beg, $enlarged_spike->{end}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
		$strand = '-';
	}
	else {	
		$results = $detector->detect_on_strand($enlarged_spike->{strand}, $genome_seq, $enlarged_spike->{begin}-$genome_seq_beg, $enlarged_spike->{end}-$genome_seq_beg,
		$candidate_region[0]-$genome_seq_beg, $candidate_region[1]-$genome_seq_beg);
	}
	return 0 if scalar @$results == 0;
	my @farthest = @{$results->[0]};
	foreach my $r (@$results) {
		@farthest = @$r if $r->[0] < $farthest[0]; # if begin_sequence of 'r' < begin_sequence of 'farthest'
	}
	my $other = {begin => $candidate_region[0]+$farthest[0], end => $candidate_region[0]+$farthest[1], strand => $enlarged_spike->{strand}};
	push @$miRnaPos, {strand => $strand, from_read => $enlarged_spike, detected => $other, first => $other, second => $enlarged_spike, source => DUE_TO_SINGLE_SPIKE};
	return 1;
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __look_both_ways {
	my ($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $min_length, $window_length) = @_;
	my $result_back = __look_backward($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $window_length);
	my $result_forw = __look_forward($detector, $miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $window_length);
	return $result_back || $result_forw;
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub __match_with_neighbor {
	my ($detector, $miRnaPos, $strand, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $current_spike, $neighbor, $min_length) = @_;
	$neighbor->{visited} = 1;
	my $enlarged_neighbor = __enlarged_spike($neighbor, $min_length, $chr_length);
	my $spike_len = $current_spike->{end} - $current_spike->{begin};
	my $neighbor_len = $enlarged_neighbor->{end} - $enlarged_neighbor->{begin};
	my $results;
	if ($current_spike->{begin} < $genome_seq_beg || $enlarged_neighbor->{begin} < $genome_seq_beg ||
		$current_spike->{end} > $genome_seq_end || $enlarged_neighbor->{end} > $genome_seq_end) {
		print "Matching neighbors... Genome seq: $chr:$genome_seq_beg-", $genome_seq_end-1, "\n";
		print "Too long to test $chr:$current_spike->{begin}-", $current_spike->{end}-1, " against $chr:$enlarged_neighbor->{begin}-", $enlarged_neighbor->{end}-1, "\n";
		return 0;
	}
	
	if ($spike_len < $neighbor_len) {
		$results = $detector->detect_on_strand($strand, $genome_seq, $current_spike->{begin}-$genome_seq_beg,
		$current_spike->{end}-$genome_seq_beg, $enlarged_neighbor->{begin}-$genome_seq_beg, $enlarged_neighbor->{end}-$genome_seq_beg);
		if (scalar @$results) {
			push @$miRnaPos, {strand => $strand, first => $current_spike, second => __force_spike_size($neighbor, $spike_len, $chr_length), source => DUE_TO_TWO_SPIKES};
			return 1;
		}
	}
	else {
		$results = $detector->detect_on_strand($strand, $genome_seq, $enlarged_neighbor->{begin}-$genome_seq_beg, $enlarged_neighbor->{end}-$genome_seq_beg,
		$current_spike->{begin}-$genome_seq_beg, $current_spike->{end}-$genome_seq_beg);
		if (scalar @$results) {
			push @$miRnaPos, {strand => $strand, first => $current_spike, second => __force_spike_size($neighbor, $spike_len, $chr_length), source => DUE_TO_TWO_SPIKES};
			return 1;
		}
	}
	return 0;
}


=method __get_windows_process_train_spikes

Private helper function. You shouldnt use this function.

=cut
sub process_window_spikes_for_chr {
	my ($this, $chr, $windows) = @_;

	my $accepting_time = $this->{accepting_time};
	my @miRnaPos = (); # Array of [5p, 3p]
	my $genome = $this->{genome_db};
	my $detector = $this->{detector};
	my $min_length = 21;
	my $min_length_for_neighbor = 40;
	my $chr_length = $this->{chr_info}{$chr}+1;
	my $min_length_on_isolated_spike = 21;
	my $window_length = 400;
	my $enlarging = max($window_length, $accepting_time);
	my $min_error_for_single_spike = $detector->miRnaMinErrorsThreshold();
	my $min_error_for_two_spikes = $detector->miRnaMaxErrorsThreshold();
	foreach my $window (@$windows) {
		foreach my $train (@{$window->{trains}}) {
			my $spikes = $train->{spikes};
			for (my $i = 0, my $e = scalar @$spikes; $i < $e; $i++) {
				my $current_spike = $spikes->[$i];
				my $unknown_strand = $current_spike->{strand} eq '?';
				my $enlarged_spike = __force_spike_size($current_spike, $min_length, $chr_length);
				my $genome_seq_beg = max($enlarged_spike->{begin}-50, 1);
				my $genome_seq_end = min($enlarged_spike->{end}+$accepting_time*2+100, $chr_length);
				if ($genome_seq_end-$genome_seq_beg > $detector->admissibleTextLength()) {
					$genome_seq_end = $genome_seq_beg+$detector->admissibleTextLength();
				}
				my $genome_seq = $genome->seq($chr, $genome_seq_beg, $genome_seq_end-1);
				my $has_neighbor = 0; # 0 == no neighbor, -1 = neighbor of other strand, 1 : neighbor on same strand
				my $has_results = 0;
				$detector->setMiRnaMinErrorsThreshold($min_error_for_two_spikes);
				for (my $j = $i+1; $j < $e; $j++) {
					my $neighbor = $spikes->[$j];
					if ($neighbor->{begin} - $enlarged_spike->{end} < $accepting_time) {
						if ($unknown_strand || $neighbor->{strand} eq $current_spike->{strand}) {
							$has_neighbor = 1;
						}
						else {
							$has_neighbor = -1 if $has_neighbor == 0;
							next;
						}
					}
					else {
						last;
					}
					if ($has_neighbor == 1) {
						if ($unknown_strand) {
							$has_results = __match_with_neighbor($detector, \@miRnaPos, '+', $genome_seq_beg, $genome_seq_end, $genome_seq,
											$chr, $chr_length, $enlarged_spike, $neighbor, $min_length_for_neighbor);
# 							if (!$has_results) {
								my $has_results_2 = __match_with_neighbor($detector, \@miRnaPos, '-', $genome_seq_beg, $genome_seq_end, $genome_seq,
											$chr, $chr_length, $enlarged_spike, $neighbor, $min_length_for_neighbor);
								$has_results = $has_results || $has_results_2;
# 							}
						}
						else {
							$has_results = __match_with_neighbor($detector, \@miRnaPos, $enlarged_spike->{strand}, $genome_seq_beg, $genome_seq_end, $genome_seq,
											$chr, $chr_length, $enlarged_spike, $neighbor, $min_length_for_neighbor);
						}
					}
				}
				next if $has_results;
				$detector->setMiRnaMinErrorsThreshold($min_error_for_single_spike);

				$genome_seq_beg = max($enlarged_spike->{begin}-$window_length-15, 1);
				$genome_seq_end = min($enlarged_spike->{end}+$window_length+15, $chr_length);
				$genome_seq = $genome->seq($chr, $genome_seq_beg, $genome_seq_end-1);
				
				__look_both_ways($detector, \@miRnaPos, $genome_seq_beg, $genome_seq_end, $genome_seq, $chr, $chr_length, $enlarged_spike,
				$min_length_on_isolated_spike, $window_length);
					
			}
		}
	}
	$detector->setMiRnaMinErrorsThreshold($min_error_for_single_spike);
	return \@miRnaPos;
}


=method __get_windows_process_train_spikes

Private helper function. You shouldnt use this function.

=cut
# args = (this, windows_ref, begin, end, read_count, forward_read_count
sub __add_candidate_window_from_train {
	my ($this, $windows, $window) = @_;
	if ($window->{forward_read_count} >= $window->{read_count}*.7) {
		push @$windows, {begin => $window->{begin}, end => $window->{end}, strand => '+', read_count => $window->{read_count}, trains => $window->{trains}};
	}
	elsif ($window->{forward_read_count} <= $window->{read_count}*.3) {
		push @$windows, {begin => $window->{begin}, end => $window->{end}, strand => '-', read_count => $window->{read_count}, trains => $window->{trains}};
	}
	else {
		push @$windows, {begin => $window->{begin}, end => $window->{end}, strand => '?', read_count => $window->{read_count}, trains => $window->{trains}};
	}
}


=method __get_windows_process_train_spikes

Private helper function. You shouldnt use this function.

=cut
sub __get_windows_process_train {
	my ($this, $windows, $current_train, $window) = @_;
	my $accepting_time = $this->{accepting_time};
	if ($current_train->{read_count} >= $this->{threshold}) { # This train should be in a window
		my $avg_begin = floor($current_train->{begin} + $current_train->{begin_offset_count}/$current_train->{read_count});
		my $avg_end = $current_train->{begin} + ceil($current_train->{length_count}/$current_train->{read_count});
		if ($window->{begin} == $window->{end}) { # There were no current window
			$window->{begin} = $current_train->{begin};
			$window->{end} = $current_train->{end};
			$window->{read_count} = $current_train->{read_count};
			$window->{forward_read_count} = $current_train->{forward_read_count};
			$window->{trains} = [{begin => $current_train->{begin}, end => $current_train->{end}, last_read_begin => $current_train->{last_read_begin},
			read_count => $current_train->{read_count}, forward_read_count => $current_train->{forward_read_count}, avg_begin => $avg_begin, avg_end => $avg_end}];
		}
		elsif ($current_train->{begin} - $window->{end} < $accepting_time) { # We accept this train as it his not far from the last window
			$window->{end} = $current_train->{end};
			$window->{read_count} += $current_train->{read_count};
			$window->{forward_read_count} += $current_train->{forward_read_count};
			push @{$window->{trains}}, {begin => $current_train->{begin}, end => $current_train->{end}, last_read_begin => $current_train->{last_read_begin},
			read_count => $current_train->{read_count}, forward_read_count => $current_train->{forward_read_count}, avg_begin => $avg_begin, avg_end => $avg_end};
		}
		else { # This train should go on a separate window
			$this->__add_candidate_window_from_train($windows, $window);
			$window->{begin} = $current_train->{begin};
			$window->{end} = $current_train->{end};
			$window->{read_count} = $current_train->{read_count};
			$window->{forward_read_count} = $current_train->{forward_read_count};
			$window->{trains} = [{begin => $current_train->{begin}, end => $current_train->{end}, last_read_begin => $current_train->{last_read_begin},
			read_count => $current_train->{read_count}, forward_read_count => $current_train->{forward_read_count}, avg_begin => $avg_begin, avg_end => $avg_end}];
		}
	}
}


=method compute_candidate_precursors_from_miRnaPos

Computes candidate precursors based miRNA candidates. This is made by adding 100 nt on both sides of the miRna couple.
Regions that overlap by more than 60% are merged together.

 Usage : my $candidate_precursors = $self->compute_candidate_precursors_from_miRnaPos($miRnas),
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
sub compute_candidate_precursors_from_miRnaPos {
	my $this = shift;
	my $miRnaPosPerChr = shift;
	my %candidate_region = ();

	foreach my $chr (keys %{ $this->{chr_info} }) {
		print "\t", $chr, "\n";
		my $miRnaPos = $miRnaPosPerChr->{$chr};
		$candidate_region{$chr} = compute_candidate_precursors_from_miRnaPos_for_chr($miRnaPos, $this->{chr_info}{$chr});
	}
	return \%candidate_region;
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub region_intertects {
	my ($a, $b, $ratio) = @_;
	if ($a->{begin} <= $b->{begin} && $a->{end} > $b->{begin}) {
		if ($b->{end} <= $a->{end}) {
			return 1;
		}
		return ($a->{end}-$b->{begin})/($b->{end}-$a->{begin}) >= $ratio ? 1 : 0;
	}
	elsif ($b->{begin} <= $a->{begin} && $b->{end} > $a->{begin}) {
		if ($a->{end} <= $b->{end}) {
			return 1;
		}
		return ($b->{end}-$a->{begin})/($a->{end}-$b->{begin}) >= $ratio ? 1 : 0;
	}
	return 0;
}


=method __get_windows_process_train_spikes

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


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
# $a and $b must intersect
# the result is merged in $a
sub merge_regions {
	my ($a, $b) = @_;
	if ($a->{begin} <= $b->{begin}) {
		if ($b->{end} > $a->{end}) {
			$a->{end} = $b->{end};
		}
	}
	else { # ($b->{begin} <= $a->{begin})
		if ($b->{end} <= $a->{end}) {
			$a->{begin} = $b->{begin};
		}
		else {
			$a->{begin} = $b->{begin};
			$a->{end} = $b->{end};
		}
	}
	push @{$a->{miRnas}}, $b->{miRnas}[0];
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub compute_candidate_precursors_from_miRnaPos_for_chr {

	sub binary_insert {
		my $inter_ratio = 0.6;
		my ($ary, $val) = @_;
		my ($min, $max, $middle, $last) = (0, scalar @{$ary}, 0, -1);
		if ($max > 0) {
			while (1) {
				$middle = int(($min+$max)/2);
				if ($min >= $max) {
					last;
				}
				if ($ary->[$middle]{begin} < $val->{begin}) {
					$min = $middle+1;
				}
				elsif ($ary->[$middle]{begin} > $val->{begin}) {
					$max = $middle;
				}
				else {
					last;
				}
			}
		}
		if ($middle < scalar @{$ary}) {
			if (region_intertects($ary->[$middle], $val, $inter_ratio)) {
				merge_regions($ary->[$middle], $val);
				return;
			}
		}
		if ($middle > 0) {
			if (region_intertects($ary->[$middle-1], $val, $inter_ratio)) {
				merge_regions($ary->[$middle-1], $val);
				return;
			}
		}
		splice @$ary, $middle, 0, $val;
	}

	my $miRnaPos = shift;
	my $chr_length = shift;

	my $padding = 100;

	my %regions = ('+' => [], '-' => []);
	if (scalar @$miRnaPos == 0) {
		return $regions{'+'};
	}

	my $current_miRna = $miRnaPos->[0];
	my %last_region = ('+' => undef, '-' => undef);

	push @{$regions{$current_miRna->{strand}}}, {begin => max(1, $current_miRna->{first}{begin}-$padding), end => min($chr_length, $current_miRna->{second}{end}+$padding),
	strand => $current_miRna->{strand}, miRnas => [$current_miRna]};


	for (my $i = 1, my $e = scalar @$miRnaPos; $i < $e; $i++) {
		my %region;
		my $current_miRna = $miRnaPos->[$i];
		%region = (begin => max(1, $current_miRna->{first}{begin}-$padding), end => min($chr_length, $current_miRna->{second}{end}+$padding), strand => $current_miRna->{strand},
		miRnas => [$current_miRna]);
		binary_insert $regions{$current_miRna->{strand}}, \%region;
	}

	my @final_regions = sort {$a->{begin} <=> $b->{begin}} (@{$regions{'+'}}, @{$regions{'-'}});
	return \@final_regions;
}


=method __get_windows_process_train_spikes

Private helper function. You shouldnt use this function.

=cut
# args = (this, windows_ref, begin, end, read_count, forward_read_count
sub __add_candidate_window {
	my ($this, $windows, $begin, $end, $read_count, $forward_read_count) = @_;
	if ($forward_read_count >= $read_count*.7) {
		push @$windows, {begin => $begin, end => $end, strand => '+'};
	}
	elsif ($forward_read_count <= $read_count*.3) {
		push @$windows, {begin => $begin, end => $end, strand => '-'};
	}
	else {
		push @$windows, {begin => $begin, end => $end, strand => '?'};
	}
}



sub export_windows_to_gff {
	my ($this, $windows_per_chr, $output_folder, $gff_annotation_file) = @_;

	my $filename = "$output_folder/window_list.gff3";
	open(my $fh, '>', $filename);
	foreach my $chr (keys %{ $this->{chr_info} }) {
		my $windows = $windows_per_chr->{$chr};
		foreach my $window (@$windows) { # Already sorted
			my $trains = $window->{trains};
			my $spikes = '';
			foreach my $train (@$trains) {
				foreach my $spike (@{$train->{spikes}}) {
					$spikes .= $spike->{begin}. '-' . ($spike->{end}-1) . ';';
				}
			}
			print $fh $chr,"\t.\twindow\t",$window->{begin},"\t",$window->{end}-1,"\t.\t", $window->{strand} ,"\t.\tReads=" .$window->{read_count}. ";Trains=" . scalar(@$trains) .
			";" . $spikes . "\n";
# 					  1		  2	 3			4					  5					 6  7  8  9
		}
	}
	close $fh;

# 	if ($gff_annotation_file) {
# 		system("intersectBed -a $gff_annotation_file -b $filename -v > $output_folder/non_detected_windows_$this->{threshold}.csv");
# 		system("intersectBed -b $gff_annotation_file -a $filename -v > $output_folder/false_positive_windows_$this->{threshold}.csv");
# 	}
}


sub export_miRnaPos_to_gff {
	my ($this, $miRnaPosPerChr, $output_folder, $gff_annotation_file) = @_;

	my $filename = "$output_folder/miRna_list.gff3";
	open(my $fh, '>', $filename);
	my $counter = 0, my $counter_single = 0;
	foreach my $chr (keys %{ $this->{chr_info} }) {
		my $miRnaPos = $miRnaPosPerChr->{$chr};
		foreach my $pair (@$miRnaPos) { # Already sorted
			my $miRNA_5p = $pair->{strand} eq '+' ? $pair->{first} : $pair->{second};
			print $fh $chr,"\t.\tmiRna-5p\t",$miRNA_5p->{begin},"\t",$miRNA_5p->{end}-1,"\t.\t", $pair->{strand} ,"\t.\t",
			$pair->{source} == DUE_TO_SINGLE_SPIKE ? 'Origin=Single;' : 'Origin=Couple;', "\n";
			my $miRNA_3p = $pair->{strand} eq '+' ? $pair->{second} : $pair->{first};
			print $fh $chr,"\t.\tmiRna-3p\t",$miRNA_3p->{begin},"\t",$miRNA_3p->{end}-1,"\t.\t", $pair->{strand} ,"\t.\t",
			$pair->{source} == DUE_TO_SINGLE_SPIKE ? 'Origin=Single;' : 'Origin=Couple;', "\n";
			$counter++;
			$counter_single++ if $pair->{source} == DUE_TO_SINGLE_SPIKE;
		}
	}
	close $fh;

# 	print $counter_single/$counter, "% of miRna are detected from a single spike.";

	if ($gff_annotation_file) {
		system("intersectBed -a $gff_annotation_file -b $filename -v -s > $output_folder/non_detected_miRna.csv");
		system("intersectBed -a $gff_annotation_file -b $filename -wa -u -s > $output_folder/detected_miRna.csv");
		system("intersectBed -b $gff_annotation_file -a $filename -v -s > $output_folder/false_positive_miRna.csv");
	}
}

sub export_precursors_to_gff {
	my ($this, $regionsPerChr, $output_folder, $output_file, $gff_annotation_file, $bed_file) = @_;

	my $filename = "$output_folder/$output_file";
	open(my $fh, '>', $filename.'.csv');
	my $counter = 0, my $counter_single = 0;
	foreach my $chr (keys %{ $this->{chr_info} }) {
		my $regions = $regionsPerChr->{$chr};
# 		print $chr, "\n";
		foreach my $region (@$regions) { # Already sorted
			print $fh $chr,"\t.\tmiRna-precursor\t",$region->{begin},"\t",$region->{end}-1,"\t.\t", $region->{strand} ,"\t.\t";
			foreach my $miRna (@{$region->{miRnas}}) {
				if ($miRna->{strand} eq '+') {
					print $fh "{5p=$chr:", $miRna->{first}{begin}, "-", $miRna->{first}{end}-1, ",";
					print $fh "3p=$chr:", $miRna->{second}{begin}, "-", $miRna->{second}{end}-1, "};";
				}
				else {
					print $fh "{5p=$chr:", $miRna->{second}{begin}, "-", $miRna->{second}{end}-1, ",";
					print $fh "3p=$chr:", $miRna->{first}{begin}, "-", $miRna->{first}{end}-1, "};";
				}
			}
			print $fh "\n";
		}
	}
	close $fh;

	my @files = ();
	open (my $fh2, '<', $filename.".csv");
	open (my $fh_out, '>', $filename."_.csv");
	while (<$fh2>) {
		chomp;
		my @f = split("\t");
		splice @f, 8, 1;
		print $fh_out join("\t", @f), "\n";
	}
	close $fh2;
	close $fh_out;
	__intersectBed($filename.'_.csv', $bed_file, $filename.'_count.csv');
	if ($gff_annotation_file) {
		system("intersectBed -a $gff_annotation_file -b $filename.csv -v -s > $output_folder/non_detected_precursors__against__$output_file.csv");
		system("intersectBed -a $gff_annotation_file -b $filename.csv -wa -u -s > $output_folder/detected_precursors__against__$output_file.csv");
		system("intersectBed -b $gff_annotation_file -a ".$filename."_.csv -v -s > $output_folder/false_positive_precursor__against__$output_file.csv");
		@files = ("$output_folder/false_positive_precursor__against__$output_file",
		"$output_folder/non_detected_precursors__against__$output_file", "$output_folder/detected_precursors__against__$output_file");
	}

	foreach my $file (@files) {
		__intersectBed($file.'.csv', $bed_file, $file.'_count.csv');
	}
}

sub __intersectBed {
# Same order of chromosomes
	my $a = shift;
	my $b = shift;
	my $out = shift;
	my $chunk_size = 30 << 20;
	my $fileCount = 1;

	if (!-e $b. '_0.bed') {
		open (my $fh_b, '<', $b);
		open (my $fh_bout, '>', $b. '_0.bed');
		my $size = 0;

		while (<$fh_b>) {
			$size += length $_;
			if ($size > $chunk_size) {
				close $fh_bout;
				open ($fh_bout, '>', $b. '_' .$fileCount. '.bed');
				$fileCount++;
				$size = 0;
			}
			print $fh_bout $_;
		}
		close $fh_b;
		close $fh_bout;
	}
	else {
		while (-e $b. '_' .$fileCount. '.bed') {
			$fileCount++;
		}
	}

	for (my $i = 0; $i < $fileCount; $i++) {
		system("intersectBed -a $a -b $b" . '_' .$i. '.bed -c -s > '.$a.'_bed_result_' . $i . '.bed');
	}
	my $cat_cmd = 'cat';
	for (my $i = 0; $i < $fileCount; $i++) {
		$cat_cmd .= ' '.$a . '_bed_result_' .$i. '.bed';
	}
	$cat_cmd .= ' | sort > ' .$out. '_';
	system($cat_cmd);
	
	for (my $i = 0; $i < $fileCount; $i++) {
		unlink $a . '_bed_result_' .$i. '.bed';
	}

	open (my $fh_out_, '<', $out. '_');
	open (my $fh_out, '>', $out);
	my @last_line = ();
	while (<$fh_out_>) {
		chomp;
		my @fields = split("\t");
		if (scalar(@fields) == 0) {
			next;
		}
		if (scalar(@last_line) == 0) {
			@last_line = @fields;
			$last_line[-1] = 0;
		}
		my $all_equal = 1;
		for (my $i = 0, my $e = scalar(@fields)-1; $i < $e; $i++) {
			if ($fields[$i] ne $last_line[$i]) {
				$all_equal = 0;
				next;
			}
		}
		if ($all_equal == 1) {
			$last_line[-1] += $fields[-1];
		}
		else {
			print $fh_out join("\t", @last_line), "\n";
			@last_line = @fields;
		}
	}
	if (scalar(@last_line) > 0) {
		print $fh_out join("\t", @last_line), "\n";
	}
	close $fh_out_;
	close $fh_out;
	unlink $out. '_';
}


=method __get_windows_process_train_spikes

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


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub run_rnastemloop {
	my ( $input, $output_stemloop, $output_optimal ) = @_;
	return miRkwood::Programs::run_rnastemloop($input, $output_stemloop, $output_optimal);
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
# returns an array of {begin, end, sequence, secondary_structure}
sub parse_stemloop_output {
	my $seq_name = shift;
	my $region = shift;
	my $safe_seq_name = quotemeta $seq_name;
	my $input = shift;
	my @results = ();
	open (my $fh, '<', $input) or die "Can't open $input: $!";
	if ($region->{strand} eq '+') {
		my $seq_begin = $region->{begin};
		while (my $line = <$fh>) {
			chomp $line;
			# if $line starts with ">$seqname" . "__DIGITS-DIGITS"
			if ($line =~ /^>${safe_seq_name}__(\d+)-(\d+)/) {
				my %entry = (begin => $1 + $seq_begin-1, end => $2 + $seq_begin);
				$line = <$fh>;
				chomp $line;
				$entry{sequence} = $line;
				$line = <$fh>;
				chomp $line;
				$entry{secondary_structure} = $line;
				push @results, \%entry;
			}
		}
	}
	else {
		my $seq_begin = $region->{begin};
		my $seq_len = $region->{end} - $seq_begin;
		while (my $line = <$fh>) {
			chomp $line;
			# if $line starts with ">$seqname" . "__DIGITS-DIGITS"
			if ($line =~ /^>${safe_seq_name}__(\d+)-(\d+)/) {
				my %entry = (begin => ($seq_len - $2 + $seq_begin), end => ($seq_len - $1 +1 + $seq_begin));
				$line = <$fh>;
				chomp $line;
				$entry{sequence} = $line;
				$line = <$fh>;
				chomp $line;
				$entry{secondary_structure} = $line;
				push @results, \%entry;
			}
		}
	}
	return \@results;
}


=method __get_windows_process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub eval_stemloop_output {
	my $stemloop_results = shift;
	my $miRnaPos = shift;

	my @retained_stemloops = ();

	foreach my $stemloop (@$stemloop_results) {
		foreach my $miRnas (@$miRnaPos) {
			if ($miRnas->{source} == DUE_TO_TWO_SPIKES) {
				if (raw_regions_intertect($stemloop, $miRnas->{first}) && raw_regions_intertect($stemloop, $miRnas->{second})) {
					push @retained_stemloops, $stemloop;
					last;
				}
			}
			else {
				if (raw_regions_intertect($stemloop, $miRnas->{from_read})) {
					push @retained_stemloops, $stemloop;
					last;
				}
			}
		}
	}
	return \@retained_stemloops;
}


=method apply_structure_criterion

Runs rnalfold and rnastemloop on the candidate precursors returned by compute_candidate_precursors_from_miRnaPos.
This elimates candidates where no stem loops are possible, or where the detected miRnas fall outside the rnastemloop-predicted precursors.

 Usage : my $retained_regions = $self->apply_structure_criterion($candidate_precursors),
 Input : The candidate precursors returned by compute_candidate_precursors_from_miRnaPos.
 Return: A hash ref {
			chr => [array of precursors {
				begin => start position (1-based). This is left unchanged from the passed $candidate_precursors.
				end => end position (excluded). This is left unchanged from the passed $candidate_precursors.
				strand => '+' or '-'. This is left unchanged from the passed $candidate_precursors.
				miRnas => [] Array reference of miRNA candidates contained in the precursor. This is left unchanged from the passed $candidate_precursors.
				stemloop => [array of valid secondary structures {
					begin => start position (1-based). This is left unchanged from the passed $candidate_precursors.
					end => end position (excluded). This is left unchanged from the passed $candidate_precursors.
					sequence => genomic sequence.
					secondary_structure => The secondary structure returned from rnastemloop
					}]
				}]
			}

=cut
sub apply_structure_criterion {
	my ($this, $regionsPerChr, $cache_folder) = @_;
	my %retained_regions_per_chr = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
		print "\t", $chr, "\n";
		my $regions = $regionsPerChr->{$chr};
		$retained_regions_per_chr{$chr} = $this->apply_structure_criterion_per_chr($chr, $regions, $cache_folder);
	}
	return \%retained_regions_per_chr;
}


=method __get_windows_process_train_spikes

Private helper function. You shouldnt use this function.

=cut
sub apply_structure_criterion_per_chr {
	my ($this, $chr, $regions, $cache_folder) = @_;
	my $genome = $this->{genome_db};
	my @retained_regions = ();
	my $eliminated_by_stemloop = 0;
	foreach my $region (@$regions) {
		my $strand = $region->{strand} eq '+' ? 1 : -1;
		my $seq_filename = $region->{strand}. '_' . $chr . '_' . $region->{begin} . '-' . ($region->{end}-1);
		my $rnalfold_output_filename = "rnalfold_out_" . $seq_filename;
		my $rnalfold_output_filepath = $cache_folder.'/rnalfold_cache/'.$rnalfold_output_filename;
		if (!-e $rnalfold_output_filepath) {
			$rnalfold_output_filepath = run_rnalfold('miRnaPrecursor', $genome->seq($chr, $region->{begin}, $region->{end}-1, $strand),
			$cache_folder.'/rnalfold_cache', $rnalfold_output_filename);
		}
		my $output_stemloop_opt = $cache_folder.'/rnastemloop_cache/rnastemloop_' .$seq_filename. '_opt';
		my $output_stemloop = $cache_folder.'/rnastemloop_cache/rnastemloop_' .$seq_filename;

		if (!-e $output_stemloop) {
			run_rnastemloop($rnalfold_output_filepath, $output_stemloop, $output_stemloop_opt);
		}
		my $results = parse_stemloop_output('miRnaPrecursor', $region, $output_stemloop);
		$eliminated_by_stemloop++ if scalar @$results == 0;
		my $retained_stemloops = eval_stemloop_output($results, $region->{miRnas});
		if (scalar(@$retained_stemloops) != 0) {
			$region->{stemloop} = $retained_stemloops;
			push @retained_regions, $region;
		}
	}
	print "\tStemloop eliminated $eliminated_by_stemloop regions.\n";
	return \@retained_regions;
}

1;
