# ABSTRACT: Handling the cluster making from a BAM file.

package Clusters;

use strict;
use warnings;
use POSIX;
use MiRnaDuplexDetector;
use Bio::DB::Fasta;
use KMean;

use List::Util qw(max min);


=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my $bam_file    = shift @args;
    my $genome_file = shift @args;
    my $genome_db = Bio::DB::Fasta->new($genome_file);
    my $self = bless {
        bam_file => $bam_file,
        genome_file => $genome_file,
        genome_db => $genome_db
    }, $class;
    my %chr_info = $self->get_chromosomes_info_from_genome_file();
    $self->{chr_info} = \%chr_info;
    return $self;
}


# Returns a hash (chr_name => read positions)
sub get_read_loci_per_chr {
    my ( $self, @args ) = @_;
    # go chr by chr, using the .fai index file to get the chr names
    my %reads = ();
    my @chrs = keys %{ $self->{chr_info} };

    foreach my $chr (@chrs) {
        my $samtools_cmd = "samtools view -F 0x4 $self->{bam_file} $chr";
        open( my $DEPTH, '-|', $samtools_cmd ) or die "Can't open '$samtools_cmd': $!";
        my $chr_reads = $self->process_samtools_view( $DEPTH, $chr);
        $reads{$chr} = $chr_reads;
        close $DEPTH;
    }
    return \%reads;
}


sub compute_read_distribution_per_chr_from_read_loci {
	my ($this, $reads_per_chr) = @_;
	my @chrs = keys %{ $this->{chr_info} };
    my %read_distribution_per_chr = ();

    foreach my $chr (@chrs) {
    # Hash of $read_position -> (Amount of reads starting at that position, Max end pos of reads)
		my %read_distribution = ();
		my $reads = $reads_per_chr->{$chr};
        foreach my $read (@$reads) {
			my $read_beg = $read->[0];
			my $read_end = $read->[1];
			if (defined $read_distribution{$read_beg}) {
				$read_distribution{$read_beg}->[0]++;
				$read_distribution{$read_beg}->[1] = max $read_distribution{$read_beg}->[1], $read_end;
			}
			else {
				$read_distribution{$read_beg} = [1, $read_end];
			}
        }
        $read_distribution_per_chr{$chr} = \%read_distribution;
    }
    return \%read_distribution_per_chr;
}

# Returns a hash (chr_name => read positions)
sub get_read_distribution_per_chr {
    my ( $self, @args ) = @_;
    # go chr by chr, using the .fai index file to get the chr names
    my %reads = ();
    my @chrs = keys %{ $self->{chr_info} };

    foreach my $chr (@chrs) {
        my $samtools_cmd = "samtools view -F 0x4 $self->{bam_file} $chr";
        open( my $DEPTH, '-|', $samtools_cmd ) or die "Can't open '$samtools_cmd': $!";
        my $chr_reads = $self->process_samtools_view_for_distribution( $DEPTH, $chr);
        $reads{$chr} = $chr_reads;
        close $DEPTH;
    }
    return \%reads;
}


sub process_samtools_view {
    my ($self, $RESULT, $chr) = @_;
    my @chr_reads = ();
    while (<$RESULT>) {
        chomp;
        my @fields = split( "\t", $_ );
        my $strand = '+';
        $strand = '-' if $fields[1] & 0x10;
        push @chr_reads, {begin => $fields[3], end => $fields[3] + length $fields[9], strand => $strand};
    }
    my @chr_reads_sorted = sort {$a->{begin} <=> $b->{begin}} @chr_reads;
    return \@chr_reads_sorted;
}


sub process_samtools_view_for_distribution {
    my ($self, $RESULT, $chr) = @_;
    my %chr_reads = ();
    while (<$RESULT>) {
        chomp;
        my @fields = split( "\t", $_ );
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


sub compute_poisson_parameters {
    my $this = shift;
    my $reads = shift;
    my $spaces = $this->compute_read_spaces($reads);
    return $this->__estimate_poisson_parameters_per_chr($spaces);
}


sub compute_poisson_parameters_chr_free {
    my $this = shift;
    my $reads = shift;
    my $spaces = $this->compute_read_spaces_chr_free($reads);
    return $this->__estimate_poisson_parameters_chr_free($spaces);
}


sub compute_read_spaces {
    my ($self, $reads_per_chr) = @_;
    my %spaces_per_chr = ();
    foreach my $chr (keys %{ $self->{chr_info} }) {
        my $reads = $reads_per_chr->{$chr};
        my @spaces = ();
        for (my $i = 0, my $e = scalar(@$reads)-1; $i < $e; $i++) {
            push @spaces, $reads->[$i+1] - $reads->[$i];
        }
        $spaces_per_chr{$chr} = \@spaces;
    }
    return \%spaces_per_chr;
}


sub compute_read_spaces_chr_free {
    my ($self, $reads_per_chr) = @_;
    my @spaces = ();
    foreach my $chr (keys %{ $self->{chr_info} }) {
        my $reads = $reads_per_chr->{$chr};
        for (my $i = 0, my $e = scalar(@$reads)-1; $i < $e; $i++) {
            push @spaces, $reads->[$i+1] - $reads->[$i];
        }
    }
    return \@spaces;
}

# args: reads_distrib_per_chr, sliding window size
sub compute_read_density_distribution {
	my ($this, $reads_distrib_per_chr, $sliding_window_size) = @_;
	my %densities = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
        my $reads = $reads_distrib_per_chr->{$chr};
        for (my $i = 0, my $e = $this->{chr_info}->{$chr}-$sliding_window_size; $i < $e; $i++) {
			my $read_count = 0;
			for (my $j = $i, my $e2 = $i + $sliding_window_size; $j < $e2; $j++) {
				if (defined $reads->{$j}) {
					$read_count += $reads->{$j};
				}
			}
# 			if ($read_count >= 1000000) {
# 				print "Huge density detected on chr $chr, at $i\n";
# 			}
			if (defined $densities{$read_count}) {
				$densities{$read_count}++;
			}
			else {
				$densities{$read_count} = 1;
			}
        }
    }
    return \%densities;
}


sub compute_read_train_per_chr {
	my ($this, $read_loci_per_chr) = @_;
	my %trains_per_chr = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
		my $params = $this->compute_read_train_chr_free($read_loci_per_chr->{$chr}, $chr);
		$trains_per_chr{$chr} = $params;
	}
	return \%trains_per_chr;
}


sub compute_read_train_chr_free {
	my ($this, $read_loci, $chr_name) = @_;
	my @trains = ();
	my $__current_train_ref = $read_loci->[0];
	my %current_train = %$__current_train_ref; # Start pos, end pos
	$current_train{count} = 1;
	foreach my $read_locus (@$read_loci) {
		if ($current_train{begin} <= $read_locus->{begin} && $read_locus->{begin} < $current_train{end}) {
			$current_train{end} = max $current_train{end}, $read_locus->{end};
			$current_train{count}++;
		}
		else {
			push @trains, [$current_train{begin}, $current_train{end}, $current_train{count}];
			%current_train = %$read_locus;
			$current_train{count} = 1;
		}
	}
	if ($current_train{begin} < $current_train{end}) {
		push @trains, [$current_train{begin}, $current_train{end}, $current_train{count}];
	}
	return \@trains;
}

sub compute_train_distributions_per_chr {
	my ($this, $read_trains_per_chr) = @_;
	my %length_distribution = ();
	my %delay_distribution = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
		my $read_trains = $read_trains_per_chr->{$chr};
		my $e = scalar(@$read_trains)-1;
		for (my $i = 0; $i < $e; $i++) {
			my $length = $read_trains->[$i]->[1] - $read_trains->[$i]->[0];
			if (defined $length_distribution{$length}) {
				$length_distribution{$length}->[0]++;
				if ($length_distribution{$length}->[0] < 16) {
					$length_distribution{$length}->[1] .= ",$chr:$read_trains->[$i]->[0]-". ($read_trains->[$i]->[1]-1);
				}
			}
			else {
				$length_distribution{$length} = [1, ",$chr:$read_trains->[$i]->[0]-". ($read_trains->[$i]->[1]-1)];
			}
			my $delay = $read_trains->[$i+1]->[0] - $read_trains->[$i]->[1];
			if (defined $delay_distribution{$delay}) {
				$delay_distribution{$delay}->[0]++;
				if ($delay_distribution{$delay}->[0] < 16) {
					$delay_distribution{$delay}->[1] .= ",$chr:$read_trains->[$i]->[1]-". ($read_trains->[$i+1]->[0]-1);
				}
			}
			else {
				$delay_distribution{$delay} = [1, ",$chr:$read_trains->[$i]->[1]-". ($read_trains->[$i+1]->[0]-1)];
			}
		}
		# we add the last one
		my $length = $read_trains->[$e]->[1] - $read_trains->[$e]->[0];
		if (defined $length_distribution{$length}) {
			$length_distribution{$length}->[0]++;
			if ($length_distribution{$length}->[0] < 16) {
				$length_distribution{$length}->[1] .= ",$chr:$read_trains->[$e]->[0]-". ($read_trains->[$e]->[1]-1);
			}
		}
		else {
			$length_distribution{$length} = [1, ",$chr:$read_trains->[$e]->[0]-". ($read_trains->[$e]->[1]-1)];
		}
	}
	return [\%length_distribution, \%delay_distribution];
}


sub plot_read_train_distribution {
	my ($this, $args, $output_folder) = @_;
	my $length_distribution = $args->[0];
	my $delay_distribution = $args->[1];
	
	my $filename = "$output_folder/train_length_distribution.csv";
	open(my $fh, '>', $filename);
	print $fh "Train length,Amount of trains\n";
	foreach my $length (sort {$a <=> $b} keys %$length_distribution) {
		print $fh "$length,$length_distribution->{$length}->[0]$length_distribution->{$length}->[1]\n";
	}
	close $fh;
	
	$filename = "$output_folder/train_delay_distribution.csv";
	open($fh, '>', $filename);
	print $fh "Train delay,Amount of trains\n";
	foreach my $delay (sort {$a <=> $b} keys %$delay_distribution) {
		print $fh "$delay,$delay_distribution->{$delay}\n";
	}
	close $fh;
}

sub plot_read_density_distribution {
	my ($self, $density_distribution, $sliding_window_size) = @_;
	my $filename = "density_distribution.csv";
	open(my $fh, '>', $filename);
	print $fh "Density,Amount of windows (window size: $sliding_window_size)\n";
	foreach my $density (sort {$a <=> $b} keys %$density_distribution) {
		print $fh $density/$sliding_window_size,",$density_distribution->{$density}->[0]$density_distribution->{$density}->[1]\n";
	}
	close $fh;
}


sub get_read_spaces_distribution {
	my ($self, $spaces_per_chr) = @_;
	my %distribution_per_chr = ();
	foreach my $chr (keys %{ $self->{chr_info} }) {
        my $spaces = $spaces_per_chr->{$chr};
        my %space_distrib = ();
        foreach my $space (@$spaces) {
			if (defined $space_distrib{$space}) {
				$space_distrib{$space}++;# = $space_distrib{$space}+1 if defined $space_distrib{$space};
			}
			else {
				$space_distrib{$space} = 1;
			}
        }
        $distribution_per_chr{$chr} = \%space_distrib;
    }
    return \%distribution_per_chr;
}


sub get_read_spaces_distribution_chr_free {
	my ($self, $spaces_chr_free) = @_;
	my %space_distrib = ();
	foreach my $space (@$spaces_chr_free) {
		if (defined $space_distrib{$space}) {
			$space_distrib{$space}++;# = $space_distrib{$space}+1 if defined $space_distrib{$space};
		}
		else {
			$space_distrib{$space} = 1;
		}
	}
    return \%space_distrib;
}


sub plot_space_distribution {
	my ($self, $distribution_per_chr) = @_;
	foreach my $chr (keys %{ $self->{chr_info} }) {
		my $cluster_count = 0;
		my $cluster_total = 0;
		my $noise_count = 0;
		my $noise_total = 0;
		my $filename = "delays_$chr.csv";
		open(my $fh, '>', $filename);
		my $distribution = $distribution_per_chr->{$chr};
		print $fh "Delays,Number\n";
		my @spaces = keys %$distribution;
		@spaces = sort {$a <=> $b} @spaces;
		for (my $i = 0, my $e = $spaces[-1]; $i <= $e; $i++) {
			if (defined $distribution->{$i}) {
				print $fh "$i,$distribution->{$i}\n";
				if ($i > 21) {
					$noise_count+=$i;
					$noise_total+=$i*$distribution->{$i};
				}
				else {
					$cluster_count+=$i;
					$cluster_total+=$i*$distribution->{$i};
				}
			}
			else {
				print $fh "$i,0\n";
				if ($i > 21) {
					$noise_count+=$i;
				}
				else {
					$cluster_count+=$i;
				}
			}
		}
		close $fh;
	}
}


sub plot_space_distribution_chr_free {
	my ($self, $distribution_chr_free) = @_;
	my $cluster_count = 0;
	my $cluster_total = 0;
	my $noise_count = 0;
	my $noise_total = 0;
	my $filename = "delay_distribution.csv";
	open(my $fh, '>', $filename);
	my $distribution = $distribution_chr_free;
	print $fh "Delay,Amount of reads\n";
	my @spaces = keys %$distribution;
	@spaces = sort {$a <=> $b} @spaces;
	for (my $i = 0, my $e = $spaces[-1]; $i <= $e; $i++) {
		if (defined $distribution->{$i}) {
			print $fh "$i,$distribution->{$i}\n";
			if ($i > 21) {
				$noise_count+=$i;
				$noise_total+=$i*$distribution->{$i};
			}
			else {
				$cluster_count+=$i;
				$cluster_total+=$i*$distribution->{$i};
			}
		}
		else {
			print $fh "$i,0\n";
			if ($i > 21) {
				$noise_count+=$i;
			}
			else {
				$cluster_count+=$i;
			}
		}
	}
	close $fh;
}


sub __estimate_poisson_parameters_per_chr {
    my ($self, $spaces_per_chr) = @_;
    my %parameters_per_chr = ();
    foreach my $chr (keys %{ $self->{chr_info} }) {
        # [Lambda_cluster, Lambda_noise]
        my $params = $self->__estimate_poisson_parameters_chr_free($spaces_per_chr->{$chr});
        $parameters_per_chr{$chr} = $params;
    }
    return \%parameters_per_chr;
}


sub __estimate_poisson_parameters_chr_free {
    my ( $this, $spaces_non_sorted ) = @_;
    my @spaces = sort {$a <=> $b} @$spaces_non_sorted;
    # Computes an estimate of the two clusters. Only look for the largest gap between two spaces, and consider it's the future margin
    my $max_delta = $spaces[1]-$spaces[0];
    my $max_delta_pos = 0;
    for (my $i = 0, my $e = scalar(@spaces)-1; $i < $e; $i++) {
        my $delta = $spaces[$i+1] - $spaces[$i];
        if ($delta > $max_delta) {
            $max_delta = $delta;
            $max_delta_pos = $i;
        }
    }
    # Compute the averages of the clusters
    my $cluster_avg = 0, my $cluster_total = 0;
    my $cluster_count = $max_delta_pos; #cluster_count defines the margin between the two clusters
    my $noise_avg = 0, my $noise_total = 0;
    my $noise_count = scalar(@spaces) - $cluster_count;
    for (my $i = 0; $i < $cluster_count; $i++) {
        $cluster_total += $spaces[$i];
    }
    for (my $i = $cluster_count, my $e = scalar(@spaces); $i < $e; $i++) {
        $noise_total += $spaces[$i];
    }
    $cluster_avg = $cluster_total/$cluster_count;
    $noise_avg = $noise_total/$noise_count;
    while (1) {
        # Assignment step
        # In this step, for every point (i.e. blank nucleotides) we check wether it is closer to the cluster avg or the noise avg, and we reassign it if needed.
        # We try to do it smartly. As we are trying to find two clusters in 1D, the two clusters are speratated by only one point/margin.
        # This means that only the points near the margin are likely to be reassigned.
        # We start to check the points near the margin upstream. If nothing gets changed, we check downstream. Note that it's impossible that we must reassign on both directions.
        # That's because the margin will either move left or right.
        
        # Check near the border, in the noise cluster (i.e. try to move the margin right)
        my $assigned_upstream = 0;
        for (my $i = $cluster_count, my $e = scalar(@spaces); $i < $e; $i++) {
            if (abs($spaces[$i] - $cluster_avg) < abs($spaces[$i] - $noise_avg)) {
                $assigned_upstream++;
                $cluster_total += $spaces[$i];
                $noise_total -= $spaces[$i];
            }
            else {
                last;
            }
        }
        # If we couldn't move the margin right, we try to move it left
        if ($assigned_upstream == 0) {
            for (my $i = $cluster_count-1; $i >= 0; $i--) {
                if (abs($spaces[$i] - $noise_avg) < abs($spaces[$i] - $cluster_avg)) {
                    $assigned_upstream--;
                    $cluster_total -= $spaces[$i];
                    $noise_total += $spaces[$i];
                }
                else {
					last;
				}
            }
        }
        # If nothing has changed (nothing was misclassified), the algorithm converged
        if ($assigned_upstream == 0) {
            last;
        }
        $cluster_count += $assigned_upstream;
        $noise_count -= $assigned_upstream;
        # Recompute the averages
        $cluster_avg = $cluster_total/$cluster_count;
        $noise_avg = $noise_total/$noise_count;
    }
    return [1/$cluster_avg, 1/$noise_avg];
}


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


sub get_windows {
	my $this = shift;
	my $params_per_chr = shift;
	my $reads_per_chr = shift;
	my %windows_per_chr = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
# 		$windows_per_chr{$chr} = $this->__get_windows_for_chr($params_per_chr->{$chr}, $reads_per_chr->{$chr});
# 		print 'Computing for ', $chr, "\n";
		$windows_per_chr{$chr} = $this->__get_windows_for_chr2($reads_per_chr->{$chr});
	}
	return \%windows_per_chr;
}

sub get_windows_from_train_analysis {
	my $this = shift;
	my $read_loci_per_chr = shift;
	my %windows_per_chr = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
		$windows_per_chr{$chr} = $this->__get_windows_from_trains_for_chr($read_loci_per_chr->{$chr});
	}
	return \%windows_per_chr;
}

sub get_windows_from_train_analysis_with_read_distribution {
	my $this = shift;
	my $read_loci_per_chr = shift;
	my $train_detection_threshold = shift;
	my %windows_per_chr = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
# 		print $chr, "\n";
		$windows_per_chr{$chr} = $this->__get_windows_from_trains_with_read_distrib_for_chr($read_loci_per_chr->{$chr}, $train_detection_threshold);
	}
	return \%windows_per_chr;
}


sub get_strand {
	my ($forward_read_count, $read_count) = @_;
	return $forward_read_count >= $read_count*.7 ? '+' : $forward_read_count <= $read_count*.3 ? '-' : '?';
}


sub __get_windows_from_trains_with_read_distrib_for_chr {
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

sub __get_windows_maintain_read_count_local_max {
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
			if ($$total_read_count - $$last_read_count <= -$spike_detection) {
				$last_spike->{end} = $trigger;
				$last_spike->{strand} = get_strand $last_spike->{forward_read_count}, $last_spike->{read_count};
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
				return;
			}
		}
		push @$spikes, {begin => $position, end => -1, trigger => $read_locus->{read_count}, read_count => $read_locus->{read_count},
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

sub __get_windows_process_train_spikes_local_max {
	my ($current_train, $position, $read_locus, $total_read_count, $end_reads, $last_read_count, $spike_detection, $seeking_pos) = @_;
	my $spikes = $current_train->{spikes};

	# Maintaining read count
	__get_windows_maintain_read_count($current_train, $position, $total_read_count, $last_read_count, $end_reads, $spike_detection);
	$$total_read_count += $read_locus->{read_count};
	push @$end_reads, {end => $read_locus->{end}, read_count => $read_locus->{read_count}};
	@$end_reads = sort {$a->{end} <=> $b->{end}} @$end_reads;
	# Looking for spikes
	if ($$total_read_count - $$last_read_count >= $spike_detection) {
		if (scalar @$spikes) {
			my $last_spike = $spikes->[-1];
			if (($position - $last_spike->{begin} < 20 && $last_spike->{trigger} < $$total_read_count)) {
				pop @$spikes;
			}
			elsif ($last_spike->{end} == -1) {
				return;
			}
		}
		push @$spikes, {begin => $position, end => -1, trigger => $read_locus->{read_count}, read_count => $read_locus->{read_count},
		forward_read_count => $read_locus->{forward_read_count}};
	}
	elsif ($$total_read_count - $$last_read_count <= -$spike_detection && scalar @$spikes) {
		my $last_spike = $spikes->[-1];
		if ($last_spike->{end} == -1) {
			$last_spike->{end} = $position;
			$last_spike->{strand} = get_strand $last_spike->{forward_read_count}, $last_spike->{read_count};
		}
	}
	elsif (scalar @$spikes && $spikes->[-1]{end} == -1) {
		my $last_spike = $spikes->[-1];
		$last_spike->{read_count} += $read_locus->{read_count};
		$last_spike->{forward_read_count} += $read_locus->{forward_read_count};
	}
}


sub __get_windows_process_train_from_distribution {
	my ($this, $windows, $current_train, $window, $train_detection_threshold) = @_;
	my $accepting_time = 300;
	if ($current_train->{read_count} >= $train_detection_threshold) { # This train should be in a window
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


sub process_window_spikes {
	my $this = shift;
	my $windows_per_chr = shift;
	my $w_len = shift;
	my $mismatches = shift;
	$w_len = 12 if !defined $w_len;
	$mismatches = 2 if !defined $mismatches;
	my %miRnaPos = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
		$miRnaPos{$chr} = $this->process_window_spikes_for_chr($chr, $windows_per_chr->{$chr}, $w_len, $mismatches);
	}
	return \%miRnaPos;
}


sub __enlarged_spike {
	my ($spike, $min_length, $chr_length) = @_;
	my $spike_length = $spike->{end} - $spike->{begin};
	if ($spike_length < $min_length) {
		my $new_spike = {begin => $spike->{begin}, end => $spike->{end}};
		my $length_by2 = int(($min_length - $spike_length)/2);
		$new_spike->{begin} = max (1, ($new_spike->{begin} - $length_by2));
		$new_spike->{end} = min ($chr_length, ($new_spike->{end} + $min_length - $spike_length - $length_by2));
		return $new_spike;
	}
	return $spike;
}


sub __look_backward_only {
	my ($miRnaPos, $genome, $chr, $chr_length, $current_spike, $min_length, $window_length, $mismatches) = @_; # min_length = 40, window_length = 300
	my $enlarged_spike = __enlarged_spike($current_spike, $min_length, $chr_length);
	__look_backward($miRnaPos, $genome, $chr, $chr_length, $current_spike, $enlarged_spike, $genome->seq($chr, $enlarged_spike->{begin}, $enlarged_spike->{end}), $window_length, $mismatches);
}

sub __look_forward {
	my ($miRnaPos, $genome, $chr, $chr_length, $current_spike, $enlarged_spike, $spike_seq, $window_length, $mismatches) = @_;
	my @candidate_region = (min ($current_spike->{begin}+15, $chr_length), min ($current_spike->{begin}+$window_length-1, $chr_length));
	if ($candidate_region[1] - $candidate_region[0] < $enlarged_spike->{end} - $enlarged_spike->{begin}) {
		return;
	}
	my $results = MiRnaDuplexDetector::edit_distance_on_strand($current_spike->{strand}, $spike_seq, $genome->seq($chr, $candidate_region[0], $candidate_region[1]), $mismatches);
	# $results = [ [[begin_read, end_read], [begin_sequence, end_sequence]], [[begin_read, end_read], [begin_sequence, end_sequence]], ...]
	if (scalar @$results) {
		my @farthest = @{$results->[0]};
		foreach my $r (@$results) {
			@farthest = @$r if $r->[1][0] > $farthest[1][0]; # if begin_sequence of 'r' > begin_sequence of 'farthest'
		}
		push @$miRnaPos, [$current_spike, {begin => $candidate_region[0]+$farthest[1][0], end => $candidate_region[0]+$farthest[1][1], strand => $current_spike->{strand}}];
	}
}

sub __look_backward {
	my ($miRnaPos, $genome, $chr, $chr_length, $current_spike, $enlarged_spike, $spike_seq, $window_length, $mismatches) = @_;
	my @candidate_region = (max ($current_spike->{begin}-$window_length-15, 1), max ($current_spike->{begin}-16, 1));
	if ($candidate_region[1] - $candidate_region[0] < $enlarged_spike->{end} - $enlarged_spike->{begin}) {
		return;
	}
	my $results = MiRnaDuplexDetector::edit_distance_on_strand($current_spike->{strand}, $spike_seq, $genome->seq($chr, $candidate_region[0], $candidate_region[1]), $mismatches);
	return if scalar @$results == 0;
	my @farthest = @{$results->[0]};
	foreach my $r (@$results) {
		@farthest = @$r if $r->[1][0] < $farthest[1][0]; # if begin_sequence of 'r' < begin_sequence of 'farthest'
	}
	push @$miRnaPos, [{begin => $candidate_region[0]+$farthest[1][0], end => $candidate_region[0]+$farthest[1][1], strand => $current_spike->{strand}}, $current_spike];
}

sub __look_both_ways {
	my ($miRnaPos, $genome, $chr, $chr_length, $current_spike, $min_length, $window_length, $mismatches) = @_;
	my $enlarged_spike = __enlarged_spike($current_spike, 40, $chr_length);
	my $spike_seq = $genome->seq($chr, $enlarged_spike->{begin}, $enlarged_spike->{end});
	__look_backward($miRnaPos, $genome, $chr, $chr_length, $current_spike, $enlarged_spike, $genome->seq($chr, $enlarged_spike->{begin}, $enlarged_spike->{end}), $window_length, $mismatches);
	__look_forward($miRnaPos, $genome, $chr, $chr_length, $current_spike, $enlarged_spike, $genome->seq($chr, $enlarged_spike->{begin}, $enlarged_spike->{end}), $window_length, $mismatches);
}

sub __match_with_neighbor {
	my ($miRnaPos, $genome, $chr, $chr_length, $current_spike, $current_spike_seq, $neighbor, $min_length, $mismatches) = @_;
	$neighbor->{visited} = 1;
	my $enlarged_neighbor = __enlarged_spike($neighbor, $min_length, $chr_length);
	my $neighbor_seq = $genome->seq($chr, $enlarged_neighbor->{begin}, $enlarged_neighbor->{end}-1);
# 	print "Current spike: Chr1:", $current_spike->{begin}, "-", $current_spike->{end}-1, ":", $current_spike_seq, "\n" if $current_spike->{begin} == 234017 && $chr eq 'Chr1';
# 	print "Neighbor spike: Chr1:", $neighbor->{begin}, "-", $neighbor->{end}-1, ":", $genome->seq($chr, $neighbor->{begin}, $neighbor->{end}-1), "\n" if $current_spike->{begin} == 234017 && $chr eq 'Chr1';
	my $spike_len = $current_spike->{end} - $current_spike->{begin};
	my $neighbor_len = $enlarged_neighbor->{end} - $enlarged_neighbor->{begin};
	my $results;
# 	if ($spike_len < $neighbor_len) {
		$results = MiRnaDuplexDetector::edit_distance_on_strand($current_spike->{strand}, $current_spike_seq, $neighbor_seq, $mismatches);
	# 	print scalar @$results, "\n" if $current_spike->{begin} == 234017 && $chr eq 'Chr1';
		if (scalar @$results) {
			push @$miRnaPos, [$current_spike, $neighbor];
			return 1;
		}
# 	}
# 	else {
# 		$results = MiRnaDuplexDetector::edit_distance_on_strand($current_spike->{strand}, $neighbor_seq, $current_spike_seq, $mismatches);
# 	# 	print scalar @$results, "\n" if $current_spike->{begin} == 234017 && $chr eq 'Chr1';
# 		if (scalar @$results) {
# 			push @$miRnaPos, [$current_spike, $neighbor];
# 			return 1;
# 		}
# 	}
	return 0;
}


sub process_window_spikes_for_chr {
	my ($this, $chr, $windows, $mismatches) = @_;
	
	my $accepting_time = 400;
	my @miRnaPos = (); # Array of [5p, 3p]
	my $genome = $this->{genome_db};
	my $min_length = 21;
	my $min_length_for_neighbor = 40;
	my $chr_length = $this->{chr_info}{$chr};
	my $min_length_on_isolated_spike = 40;
	my $window_length = 300;
	my $mismatches_on_isolated_spike = int($mismatches*$min_length_on_isolated_spike/$min_length);
	foreach my $window (@$windows) {
		my @spikes = ();
		foreach my $train (@{$window->{trains}}) {
			@spikes = (@spikes, @{$train->{spikes}});
		}
		for (my $i = 0, my $e = @spikes; $i <$e; $i++) {
			my $current_spike = $spikes[$i];
			if ($current_spike->{strand} eq '?') {
				next;
			}
			my $enlarged_spike = __enlarged_spike($current_spike, $min_length, $chr_length);
			my $current_spike_seq = $genome->seq($chr, $enlarged_spike->{begin}, $enlarged_spike->{end}-1);
			my $has_neighbor = 0; # 0 == no neighbor, -1 = neighbor of other strand, 1 : neighbor on same strand
			my $has_results = 0;
			for (my $j = $i+1; $j < $e; $j++) {
				my $neighbor = $spikes[$j];
				if ($neighbor->{begin} - $enlarged_spike->{end} < $accepting_time) {
					if ($neighbor->{strand} eq $current_spike->{strand}) {
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
					print "\ncoucou\n" if $current_spike->{begin} == 4182143;
					$has_results = __match_with_neighbor(\@miRnaPos, $genome, $chr, $chr_length, $current_spike, $current_spike_seq, $neighbor, $min_length_for_neighbor,
					int($mismatches*($enlarged_spike->{end}-$enlarged_spike->{begin})/$min_length));
					my $enlarged_neighbor = __enlarged_spike($neighbor, $min_length_for_neighbor, $chr_length);
					print "has results = $has_results\n" if $current_spike->{begin} == 4182143;
					print "enlarged spike = Chr1:", $enlarged_spike->{begin}, '-', $enlarged_spike->{end}-1, " => $current_spike_seq\n" if $current_spike->{begin} == 4182143;
					print "enlarged neighbor = Chr1:", $enlarged_neighbor->{begin}, '-', $enlarged_neighbor->{end}-1, " => ",
					$genome->seq($chr, $enlarged_neighbor->{begin}, $enlarged_neighbor->{end}-1), "\n" if $current_spike->{begin} == 4182143;
				}
			}
			next if $has_results;
			if ($has_neighbor == 1) {
				next if $current_spike->{visited};
				__look_backward_only(\@miRnaPos, $genome, $chr, $chr_length, $current_spike, $min_length_on_isolated_spike, $window_length, $mismatches_on_isolated_spike);
			}
			else {
				__look_both_ways(\@miRnaPos, $genome, $chr, $chr_length, $current_spike, $min_length_on_isolated_spike, $window_length, $mismatches_on_isolated_spike);
			}
		}
	}
	return \@miRnaPos;
}

sub __get_windows_from_trains_for_chr {
	my $this = shift;
	my $read_loci = shift;

	my @windows = ();
	if (scalar(@$read_loci) == 0) {
		return \@windows;
	}

	my %current_train = (begin => $read_loci->[0]->{begin}, end => $read_loci->[0]->{end}, read_count => 1, forward_read_count => $read_loci->[0]->{strand} eq '+' ? 1 : 0,
	last_read_begin => $read_loci->[0]->{begin}, begin_offset_count => 0, length_count => $read_loci->[0]->{end} - $read_loci->[0]->{begin}, spikes => []);
	my %window = (begin => 0, read_count => 0, forward_read_count => 0, end => 0, trains => []);

	# score at $i = (number_of_reads_at_$i)*$alpha + $beta;
	foreach my $read_locus (@$read_loci) {
		if ($current_train{begin} <= $read_locus->{begin} && $read_locus->{begin} < $current_train{end}) {
			$current_train{end} = max $current_train{end}, $read_locus->{end};
			$current_train{begin_offset_count} += $read_locus->{begin} - $current_train{begin};
			$current_train{length_count} += $read_locus->{end} - $current_train{begin};
			$current_train{read_count}++;
			$current_train{forward_read_count}++ if $read_locus->{strand} eq '+';
			$current_train{last_read_begin} = $read_locus->{begin};
		}
		else { # the read train ended
			$this->__get_windows_process_train(\@windows, \%current_train, \%window);
			%current_train = (begin => $read_locus->{begin}, end => $read_locus->{end}, read_count => 1, forward_read_count => $read_locus->{strand} eq '+' ? 1 : 0,
			last_read_begin => $read_locus->{begin}, begin_offset_count => 0, length_count => $read_locus->{end} - $read_locus->{begin});
		}
	}
	$this->__get_windows_process_train(\@windows, \%current_train, \%window);
	if ($window{begin} != $window{end}) {
		$this->__add_candidate_window_from_train(\@windows, \%window);
	}
	return \@windows;
}

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




sub __get_windows_process_train {
	my ($this, $windows, $current_train, $window) = @_;
	my $accepting_time = 300;
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

sub score_window {
	my $window = shift;
	my $trains = $window->{trains};
	my @sorted_trains = sort {$b->{read_count} <=> $a->{read_count}} @$trains;
	my $most_covered_train = $sorted_trains[0];
	my $second_most_covered_train = scalar @sorted_trains > 1 ? $sorted_trains[0] : undef;
	my $score = 0;
	my $best_length = 24;
	$score = ($best_length - ($most_covered_train->{avg_end} - $most_covered_train->{avg_begin}))**2;
	if ($second_most_covered_train) {
		$score += ($best_length - ($second_most_covered_train->{avg_end} - $second_most_covered_train->{avg_begin}))**2;
	}
	else {
		$score += $best_length*$best_length;
	}
	my $cov = $most_covered_train->{read_count};
	for (my $i = 2, my $e = @sorted_trains; $i < $e; $i++) {
		$score += ($sorted_trains[$i]->{read_count}/$cov)*(($sorted_trains[$i]{avg_end} - $sorted_trains[$i]{avg_begin})**2);
	}
	return $score;
}

sub __get_windows_for_chr2 {
	my $this = shift;
	my $read_distribution = shift;

	my %w_start = ('+' => -1, '-' => -1);
	my %score = ('+' => 0, '-' => 0);
	my %w_last_score_pos = ('+' => -1, '-' => -1);
	my %w_last_read_endpos = ('+' => -1, '-' => -1);
	my %w_threshold_pos = ('+' => -1, '-' => -1);
	my %w_read_count = ('+' => 0, '-' => 0);
	my %w_forward_read_count = ('+' => 0, '-' => 0);

	my @windows = ();

	# score at $i = (number_of_reads_at_$i)*$alpha + $beta;
	foreach my $current_position (sort {$a <=> $b} keys %$read_distribution) {
		my $current_reads = $read_distribution->{$current_position};
		my $read_amount = $current_reads->{read_count};
		my $forward_strand_read_amount = $current_reads->{forward_read_count};
		my $reverse_strand_read_amount = $read_amount - $forward_strand_read_amount;

# 		=========================
#	 		Forward strand
# 		=========================
		if ($forward_strand_read_amount) {
			$this->__get_windows_process_read(\@windows, \%w_start, \%w_last_read_endpos, \%score, \%w_last_score_pos, \%w_threshold_pos, \%w_read_count, \%w_forward_read_count, $current_position,
			$forward_strand_read_amount, '+');
		}
		if ($reverse_strand_read_amount) {
			$this->__get_windows_process_read(\@windows, \%w_start, \%w_last_read_endpos, \%score, \%w_last_score_pos, \%w_threshold_pos, \%w_read_count, \%w_forward_read_count, $current_position,
			$reverse_strand_read_amount, '-');
		}

		$w_read_count{'+'} += $read_amount;
		$w_read_count{'-'} += $read_amount;
		$w_forward_read_count{'+'} += $forward_strand_read_amount;
		$w_forward_read_count{'-'} += $forward_strand_read_amount;

		$w_last_score_pos{'+'} = $current_position if $forward_strand_read_amount;
		$w_last_score_pos{'-'} = $current_position if $reverse_strand_read_amount;

		$w_last_read_endpos{'+'} = $current_reads->{end} if $forward_strand_read_amount && $current_reads->{end} > $w_last_read_endpos{'+'};
		$w_last_read_endpos{'-'} = $current_reads->{end} if $reverse_strand_read_amount && $current_reads->{end} > $w_last_read_endpos{'-'};
	}
	if ($w_start{'+'} != -1 && $w_threshold_pos{'+'} != -1) {
		$this->__add_candidate_window(\@windows, $w_start{'+'}, $w_last_read_endpos{'+'}, $w_read_count{'+'}, $w_forward_read_count{'+'});
	}
	if ($w_start{'-'} != -1 && $w_threshold_pos{'-'} != -1) {
		$this->__add_candidate_window(\@windows, $w_start{'-'}, $w_last_read_endpos{'-'}, $w_read_count{'-'}, $w_forward_read_count{'-'});
	}
	return \@windows;
}

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


sub __get_windows_process_read {
	my ($this, $windows, $w_start, $w_last_read_endpos, $score, $w_last_score_pos, $w_threshold_pos, $w_read_count, $w_forward_read_count, $current_position, $strand_read_amount, $strand) = @_;

	my $L_cluster = 1;
	my $L_noise = 0.1;
	my $alpha = log($L_cluster/$L_noise);
	my $beta = $L_noise - $L_cluster;
	
	my $accepting_time = 300;
	
	my $got_nonpositive = 0;
	if ($w_start->{$strand} == -1) {
		# We start a new window
		$w_start->{$strand} = $current_position;
		$w_last_score_pos->{$strand} = $current_position if $strand_read_amount;

		$score->{$strand} = $alpha*$strand_read_amount + $beta;
		$w_threshold_pos->{$strand} = -1;
		$w_last_read_endpos->{$strand} = -1;
		$w_read_count->{$strand} = 0;
		$w_forward_read_count->{$strand} = 0;
	}
	else {
		$score->{$strand} += $beta*($current_position-1-$w_last_score_pos->{$strand}); # Decrease the score on every nucleotide with no reads
		if ($score->{$strand} <= 0) {
			$score->{$strand} = 0;
			$got_nonpositive = 1;
		}
		$score->{$strand} += $alpha*$strand_read_amount + $beta;
	}

# 	if ($current_position == 7423585 || $w_start->{$strand} == 7423585) {
# 		print "\n===============\nEnjoy\n===============\n";
# 		print 'current_position = ', $current_position, "\n";
# 		print 'strand = ', $strand, "\n";
# 		print 'strand_read_amount = ', $strand_read_amount, "\n";
# 		print 'score = ', %$score, "\n";
# 		print 'w_start = ', %$w_start, "\n";
# 		print 'w_last_read_endpos = ', %$w_last_read_endpos, "\n";
# 		print 'w_threshold_pos = ', %$w_threshold_pos, "\n";
# 		print 'w_last_score_pos = ', %$w_last_score_pos, "\n";
# 		print 'w_read_count = ', %$w_read_count, "\n";
# 		print 'w_forward_read_count = ', %$w_forward_read_count, "\n";
# 		print 'got_nonpositive = ', $got_nonpositive, "\n";
# 	}
# 	
	if (($got_nonpositive || $score->{$strand} <= 0) && !($w_threshold_pos->{$strand} != -1 && $current_position - $w_threshold_pos->{$strand} < $accepting_time)) {
		if ($w_threshold_pos->{$strand} != -1) {
			if ($w_start->{$strand} == 7423585) {
				print "HURRAY!!!\n";
			}
			$this->__add_candidate_window($windows, $w_start->{$strand}, $w_last_read_endpos->{$strand}, $w_read_count->{$strand}, $w_forward_read_count->{$strand});
		}
		$w_start->{$strand} = ($score->{$strand} > 0) ? $current_position : -1;
		$w_threshold_pos->{$strand} = -1;
		$w_last_read_endpos->{$strand} = -1;
		$w_read_count->{$strand} = 0;
		$w_forward_read_count->{$strand} = 0;
	}
	if ($score->{$strand} > $this->{threshold}) {
		$score->{$strand} = 0;
		$w_threshold_pos->{$strand} = $current_position;
	}
}


sub compute_window_length_distribution {
	my ($this, $windows_per_chr) = @_;
	my %length_distribution = ();
	my $max_chr = 15;
	foreach my $chr (keys %{ $this->{chr_info} }) {
		my $windows = $windows_per_chr->{$chr};
		my $e = scalar(@$windows);
		for (my $i = 0; $i < $e; $i++) {
			my $length = $windows->[$i]->{end} - $windows->[$i]->{begin};
			if (defined $length_distribution{$length}) {
				$length_distribution{$length}->[0]++;
				if ($length_distribution{$length}->[0] <= $max_chr) {
					$length_distribution{$length}->[1] .= ','.$chr . ':' . $windows->[$i]->{begin} . '-' . ($windows->[$i]->{end}-1);
				}
			}
			else {
				$length_distribution{$length} = [1, ','.$chr . ':' . $windows->[$i]->{begin} . '-' . ($windows->[$i]->{end}-1)];
			}
		}
	}
	return \%length_distribution;
}


sub plot_window_distribution {
	my ($this, $length_distribution, $output_folder) = @_;

	my $filename = "$output_folder/window_length_distribution_" . $this->{threshold} . ".csv";
	open(my $fh, '>', $filename);
	print $fh "Window length,Amount of windows\n";
	foreach my $length (sort {$a <=> $b} keys %$length_distribution) {
		print $fh "$length,$length_distribution->{$length}->[0]$length_distribution->{$length}->[1]\n";
	}
	close $fh;
}


sub export_windows_to_csv {
	my ($this, $windows_per_chr, $output_folder) = @_;

	my $filename = "$output_folder/window_list_" . $this->{threshold} . ".csv";
	open(my $fh, '>', $filename);
	print $fh "Chromosome,Starting position,Ending position (included),Locus\n";
	foreach my $chr (keys %{ $this->{chr_info} }) {
		my $windows = $windows_per_chr->{$chr};
		foreach my $window (@$windows) { # Already sorted
			print $fh $chr,',',$window->{begin},',',$window->{end}-1,",$chr:$window->{begin}-", $window->{end}-1, "\n";
		}
	}
	close $fh;
}


sub export_windows_to_gff {
	my ($this, $windows_per_chr, $output_folder, $gff_annotation_file) = @_;

	my $filename = "$output_folder/window_list_" . $this->{threshold} . ".gff3";
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

	my $filename = "$output_folder/miRna_list_" . $this->{threshold} . ".gff3";
	open(my $fh, '>', $filename);
	foreach my $chr (keys %{ $this->{chr_info} }) {
		my $miRnaPos = $miRnaPosPerChr->{$chr};
		foreach my $pair (@$miRnaPos) { # Already sorted
			my $miRNA_5p = $pair->[0]{strand} eq '+' ? $pair->[0] : $pair->[1];
			print $fh $chr,"\t.\tmiRna-5p\t",$miRNA_5p->{begin},"\t",$miRNA_5p->{end}-1,"\t.\t", $miRNA_5p->{strand} ,"\t.\t\n";
			my $miRNA_3p = $pair->[0]{strand} eq '+' ? $pair->[1] : $pair->[0];
			print $fh $chr,"\t.\tmiRna-3p\t",$miRNA_3p->{begin},"\t",$miRNA_3p->{end}-1,"\t.\t", $miRNA_3p->{strand} ,"\t.\t\n";
		}
	}
	close $fh;

	if ($gff_annotation_file) {
		system("intersectBed -a $gff_annotation_file -b $filename -v > $output_folder/non_detected_miRna_$this->{threshold}.csv");
		system("intersectBed -a $gff_annotation_file -b $filename -wa -u > $output_folder/detected_miRna_$this->{threshold}.csv");
		system("intersectBed -b $gff_annotation_file -a $filename -v > $output_folder/false_positive_miRna_$this->{threshold}.csv");
	}
}


sub export_windows_to_bed {
	my ($this, $windows_per_chr, $output_folder, $gff_annotation_file) = @_;

	my $filename = "$output_folder/window_list_" . $this->{threshold} . ".bed";
	open(my $fh, '>', $filename);
	foreach my $chr (keys %{ $this->{chr_info} }) {
		my $windows = $windows_per_chr->{$chr};
		foreach my $window (@$windows) { # Already sorted
			print $fh $chr,"\t",$window->{begin},"\t",$window->{end},"\n";
		}
	}
	close $fh;

	if ($gff_annotation_file) {
		system("intersectBed -a $gff_annotation_file -b $filename -v > $output_folder/non_detected_windows_$this->{threshold}.csv");
		system("intersectBed -b $gff_annotation_file -a $filename -v > $output_folder/false_positive_windows_$this->{threshold}.csv");
	}
}


1;
