# ABSTRACT: Handling the cluster making from a BAM file.

package miRkwood::ClustersSebastien;

use strict;
use warnings;
use POSIX;

use miRkwood::KMeanSebastien;

use List::Util qw(max min);


=method new

Constructor

 Usage : my $clustering = Clusters->new('genome.fa'),
 Input : The genome file
 Return: The cunstructed instance

=cut

sub new {
    my ( $class, $genome_file ) = @_;
    my $self = bless {
        genome_file => $genome_file,
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
    my %parsed_reads = ();
    my @chrs = keys %{ $this->{chr_info} };

    foreach my $chr (@chrs) {
        my $samtools_cmd = "samtools view $bam_file $chr";
        open( my $DEPTH, '-|', $samtools_cmd ) or die "Can't open '$samtools_cmd': $!";
        $reads{$chr} = __get_read_distribution_from_bam_for_chr($DEPTH);
        close $DEPTH;
    }
    return (\%reads, \%parsed_reads);
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
    my %parsed_reads = ('+' => {}, '-' => {});
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
        if (defined $parsed_reads{$strand}{$pos}) {
			$parsed_reads{$strand}{$pos}{'depth'}++;
			if (defined $parsed_reads{$strand}{$pos}{'ends'}{$end}) {
				$parsed_reads{$strand}{$pos}{'ends'}{$end}++;
			}
			else {
				$parsed_reads{$strand}{$pos}{'ends'}{$end} = 1;
			}
        }
        else {
			$parsed_reads{$strand}{$pos} = {'depth' => 1, 'ends' => {$end => 1}};
        }
    }
    return (\%chr_reads, \%parsed_reads);
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
    my %parsed_reads = ();
    foreach my $chr (keys %{$this->{chr_info}}) {
		$reads{$chr} = {};
		$parsed_reads{$chr} = {'+' => {}, '-' => {}};
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
        if (defined $parsed_reads{$chr}{$strand}{$pos}) {
			$parsed_reads{$chr}{$strand}{$pos}{'depth'} += $fields[4];
			if (defined $parsed_reads{$chr}{$strand}{$pos}{'ends'}{$end}) {
				$parsed_reads{$chr}{$strand}{$pos}{'ends'}{$end}+=$fields[4];
			}
			else {
				$parsed_reads{$chr}{$strand}{$pos}{'ends'}{$end} = $fields[4];
			}
        }
        else {
			$parsed_reads{$chr}{$strand}{$pos} = {'depth' => $fields[4], 'ends' => {$end => $fields[4]}};
        }
    }
    close $HANDLE;
    return (\%reads, \%parsed_reads);
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
	spikes => [], classifier => KMeanSebastien->new());
	$current_train{classifier}->add_point(0);
	my %window = (begin => 0, read_count => 0, forward_read_count => 0, end => 0, trains => []);
	my $last_read_count = 0;
	my $spike_detection = 2;
	my $total_read_count = 0;
	my @end_reads = ();

	foreach my $position (@positions) {
		my $read_locus = $read_distribution->{$position};
		if ($current_train{begin} <= $position && $position < $current_train{end}) {
			$current_train{end} = max $current_train{end}, $read_locus->{end};
			$current_train{begin_offset_count} += $position - $current_train{begin};
			$current_train{length_count} += $read_locus->{end} - $current_train{begin};
			$current_train{read_count} += $read_locus->{read_count};
			$current_train{forward_read_count} += $read_locus->{forward_read_count};
			$current_train{last_read_begin} = $position;
			__get_windows_process_train_spikes(\%current_train, $position, $read_locus, \$total_read_count, \@end_reads, \$last_read_count, $spike_detection);
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
			__get_windows_process_train_spikes(\%current_train, $position, $read_locus, \$total_read_count, \@end_reads, \$last_read_count, $spike_detection);
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


=method __add_candidate_window_from_train

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
			if ($added == KMeanSebastien::ASSIGNED_FIRST_CLASS) {
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
	my ($current_train, $position, $read_locus, $total_read_count, $end_reads, $last_read_count, $spike_detection) = @_;
	my $spikes = $current_train->{spikes};

	# Maintaining read count
	__get_windows_maintain_read_count($current_train, $position, $total_read_count, $last_read_count, $end_reads, $spike_detection);
	$$total_read_count += $read_locus->{read_count};
	push @$end_reads, {end => $read_locus->{end}, read_count => $read_locus->{read_count}};
	@$end_reads = sort {$a->{end} <=> $b->{end}} @$end_reads;
	# Looking for spikes
	my $added = $current_train->{classifier}->add_point($$total_read_count);
	if ($added == KMeanSebastien::ASSIGNED_SECOND_CLASS) {
		if (scalar @$spikes) {
			my $last_spike = $spikes->[-1];
			if ($last_spike->{end} == -1) {
				if ($current_train->{classifier}->class_of($last_spike->{trigger}) == KMeanSebastien::ASSIGNED_FIRST_CLASS) {
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
	elsif ($added == KMeanSebastien::ASSIGNED_FIRST_CLASS) {
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


=method __get_windows_process_train_from_distribution

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

1;
