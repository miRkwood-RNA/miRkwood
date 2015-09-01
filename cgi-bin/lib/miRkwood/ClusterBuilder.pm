package miRkwood::ClusterBuilder;

# ABSTRACT: ClusterBuilder object

use strict;
use warnings;

use parent 'miRkwood::LociBuilder';

use miRkwood::Utils;
use miRkwood::HairpinBuilder;
use miRkwood::ClusterJob;
use miRkwood::KMean;
use List::Util qw(max min);
use Log::Message::Simple qw[msg error debug];

sub new {
    my ( $class,
    $genome_db,
    $bed_file ) = @_;
    my $self = bless {
        'genome_db' => $genome_db,
        'bed_file' => $bed_file,
        'accepting_time' => 350,
        'train_detection_threshold' => 3, # How much overlapping reads (with depth) there should be to capture a train
        'loci_read_coverage_threshold' => 10, # The read coverage threshold below which the locus is discarded
        'peak_padding' => 100 # How much nt we add on both sides of miRNAs to create a locus
    }, $class;
    $self->{'chr_info'} = $self->get_chromosomes_info_from_genome_file();
    return $self;
}

sub get_parsed_bed {
	my $this = shift;
	return $this->{'parsed_bed'};
}


sub build_loci_per_chr {
    my ($this, @args) = @_;
    my $chromosome = shift @args;
    my $average_coverage = shift @args;

    debug( '    - Get read distribution...' . ' [' . localtime() . ']', miRkwood->DEBUG() );
	my ($reads_per_chr, $parsed_bed) = $this->get_read_distribution_per_chr_from_bed( $this->{'bed_file'}, $chromosome );
	$this->{'parsed_bed'} = $parsed_bed;

    debug( '    - Get trains...' . ' [' . localtime() . ']', miRkwood->DEBUG() );
	my $trains_hash_per_chr = $this->__get_trains_for_chr( $reads_per_chr );

	my $cluster_job = miRkwood::ClusterJob->new($this->{'genome_db'});
	$cluster_job->init_from_clustering($this);

    debug( '    - Get spikes...' . ' [' . localtime() . ']', miRkwood->DEBUG() );
	my $spikes_per_chr = $cluster_job->extract_spike_train_per_chr($trains_hash_per_chr);

	undef $trains_hash_per_chr;

    debug( '    - Process spikes...' . ' [' . localtime() . ']', miRkwood->DEBUG() );
	my $putative_miRna = $cluster_job->process_spikes_for_chr($chromosome, $spikes_per_chr);

	undef $spikes_per_chr;

    debug( '    - Get loci...' . ' [' . localtime() . ']', miRkwood->DEBUG() );
	my $loci_per_chr = $cluster_job->compute_candidate_precursors_from_miRnaPos_for_chr($chromosome,
                                                                                        $putative_miRna,
                                                                                        $this->{'chr_info'}{$chromosome},
                                                                                        $average_coverage,
                                                                                        $this->{'loci_read_coverage_threshold'},
                                                                                        $this->{'peak_padding'},
                                                                                        $parsed_bed);
    debug( '    - End of get loci...' . ' [' . localtime() . ']', miRkwood->DEBUG() );
	undef $putative_miRna;

	miRkwood::Utils::display_var_sizes_in_log_file( '..... ClusterBuilder : build_loci' );

	return $loci_per_chr;
}


sub add_read_info_in_loci {
	my $this = shift;
	my $loci = shift;

	foreach my $chr (keys %{$loci}) {
		foreach my $locus (@{$loci->{$chr}}) {
			$locus->{'reads'} = miRkwood::HairpinBuilder::get_contained_reads($this->get_parsed_bed, $chr, $locus->{'begin'},
			  $locus->{'end'}, $locus->{'strand'});
		}
	}

    return;
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
    die ('Not supported anymore');
    # go chr by chr, using the .fai index file to get the chr names
    #~ my %reads = ();
    #~ my %parsed_reads = ();
    #~ my @chrs = keys %{ $this->{chr_info} };
#~ 
    #~ foreach my $chr (@chrs) {
        #~ my $samtools_cmd = "samtools view $bam_file $chr";
        #~ open( my $DEPTH, '-|', $samtools_cmd ) or die "Can't open '$samtools_cmd': $!";
        #~ $reads{$chr} = __get_read_distribution_from_bam_for_chr($DEPTH);
        #~ close $DEPTH;
    #~ }
    #~ return (\%reads, \%parsed_reads);
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
	die ('Not supported anymore');
    #~ my $HANDLE = shift;
    #~ my %chr_reads = ();
    #~ my %parsed_reads = ('+' => [], '-' => []);
    #~ while (<$HANDLE>) {
        #~ chomp;
        #~ my @fields = split( /\t/ );
        #~ my $pos = $fields[3]-1;
        #~ my $strand = '+';
        #~ $strand = '-' if $fields[1] & 0x10;
        #~ my $end = $pos + length $fields[9];
        #~ if (defined $chr_reads{$pos}) {
			#~ $chr_reads{$pos}->{'read_count'}++;
			#~ $chr_reads{$pos}->{'end'} = $end if $end > $chr_reads{$pos}->{'end'};
			#~ $chr_reads{$pos}->{'forward_read_count'}++ if $strand eq '+';
        #~ }
        #~ else {
			#~ $chr_reads{$pos} = {read_count => 1, end => $end, forward_read_count => ($strand eq '+') ? 1 : 0};
        #~ }
        #~ my $added = 0;
         #~ if (scalar @{$parsed_reads{$strand}}) {
			#~ my $read_ref = $parsed_reads{$strand}[-1];
			#~ if ($read_ref->{'begin'} == $pos) {
				#~ $read_ref->{'depth'}++;
				#~ if (defined $read_ref->{'ends'}{$end}) {
					#~ $read_ref->{'ends'}{$end}++;
				#~ }
				#~ else {
					#~ $read_ref->{'ends'}{$end} = 1;
				#~ }
				#~ $added = 1;
			#~ }
        #~ }
        #~ if ($added == 0) {
			#~ push @{$parsed_reads{$strand}}, {'begin' => $pos, 'depth' => 1, 'ends' => {$end => 1}};
        #~ }
    #~ }
    #~ return (\%chr_reads, \%parsed_reads);
}


=method get_read_distribution_per_chr_from_bed

Retrieve the reads from a miRkwood-normalized bed file. The bed file must have the following columns (tab-separated):
# 0: chromosome name,
# 1: starting position (0-based, starting at 0),
# 2: end positing (excluded)
# 3: sequence name
# 4: depth (integer)
# 5: strand ('+' or '-')

 Usage : my $reads = $self->get_read_distribution_per_chr_from_bed('reads.bed', $chr),
 Input : The bed file path
 Return: A hash reference {
					begin_pos => {
						read_count => read depth,
						end => maximum end coordinates of all reads starting at begin_pos,
						forward_read_count => read depth of reads mapped onto forward strand
						}
					}

=cut
sub get_read_distribution_per_chr_from_bed {
    my ($this, $bed_file, $chromosome) = @_;
    my $reads = [];
    my $parsed_reads = { '+' => [], '-' => [] };

    open( my $HANDLE, '<', $bed_file) or die "Can't open '$bed_file': $!";
    while (<$HANDLE>) {
        chomp;
        my @fields = split( /\t/ );
        if (scalar @fields != 6) {
			next;
        }
        my $chr = $fields[0];
        my $pos = $fields[1];
        my $end = $fields[2];
        my $depth = $fields[4];
        my $strand = $fields[5];
        if ($strand ne '+' && $strand ne '-') {
			die ("Error while parsing bed: Incorrect strand, expected '+' or '-', got '$strand'.");
		}

        if ( $chr eq $chromosome ) {
            # First data structure
            if (scalar @{$reads} && $reads->[-1]{'begin'} == $pos) {
                my $locus = $reads->[-1];
                if ($locus->{'begin'} == $pos) {
                    $locus->{'read_count'} += $depth;
                    $locus->{'end'} = $end if $end > $locus->{'end'};
                    $locus->{'forward_read_count'} += $depth if $strand eq '+';
                }
            }
            else {
                push @{$reads}, { 'begin' => $pos,
                                  'read_count' => $depth,
                                  'end' => $end,
                                  'forward_read_count' => ($strand eq '+') ? $depth : 0};
            }
            # Second data structure
            my $added = 0;
            if (scalar @{$parsed_reads->{$strand}}) {
                my $read_ref = $parsed_reads->{$strand}[-1];
                if ($read_ref->{'begin'} == $pos) {
                    $read_ref->{'depth'} += $depth;
                    if (defined $read_ref->{'ends'}{$end}) {
                        $read_ref->{'ends'}{$end}+=$depth;
                    }
                    else {
                        $read_ref->{'ends'}{$end} = $depth;
                    }
                    $added = 1;
                }
            }
            if ($added == 0) {
                push @{$parsed_reads->{$strand}}, {'begin' => $pos, 'depth' => $depth, 'ends' => {$end => $depth}};
            }
        }
    }
    close $HANDLE;
    return ($reads, $parsed_reads);
}


=method get_chromosomes_info_from_genome_file

Retrieve chromosomes name and length and from FAI file.

 Usage : my %chr_info_output = $self->get_chromosomes_info_from_genome_file(),
 Input : The genome file
 Return: A hash {name => length}

=cut

sub get_chromosomes_info_from_genome_file {
    my ($self, @args) = @_;

    my $chr_lengths;
    foreach my $chr ($self->{'genome_db'}->get_all_primary_ids) {
		$chr_lengths->{$chr} = $self->{'genome_db'}->length($chr);
	}

    return $chr_lengths;
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
#~ sub get_windows {
	#~ my $this = shift;
	#~ my $read_distribution_per_chr = shift;
	#~ my $train_detection_threshold = shift;
	#~ my %windows_per_chr = ();
	#~ foreach my $chr (keys %{ $this->{chr_info} }) {
		#~ print "\t", $chr, "\n";
		#~ $windows_per_chr{$chr} = $this->__get_windows_for_chr($read_distribution_per_chr->{$chr}, $train_detection_threshold);
	#~ }
	#~ return \%windows_per_chr;
#~ }

sub get_trains {
	my $this = shift;
	my $read_distribution_per_chr = shift;
	my %trains_per_chr = ();
	foreach my $chr (keys %{ $this->{chr_info} }) {
		$trains_per_chr{$chr} = $this->__get_trains_for_chr($read_distribution_per_chr->{$chr});
	}
	return \%trains_per_chr;
}


=method get_strand

Static private helper function that returns the strand of a region based on read statistics.

 Usage : my $strand = get_strand($forward_read_count, $read_count),
 Input : The number of reads mapped onto the forward strand, and mapped on both strand
 Return: '+' (if more than 70% of reads are on the forward strand) or '-' (if more than 70% of reads are on reverse strand) or '?' (otherwise)

=cut
sub get_strand {
	my ($forward_read_count, $read_count) = @_;
	return $forward_read_count >= $read_count*.7 ? '+' : $forward_read_count <= $read_count*.3 ? '-' : '?';
}


sub __get_trains_for_chr {
	my $this = shift;
	my $read_distribution = shift;
	my $train_detection_threshold = $this->{'train_detection_threshold'};

	my @trains = ();
	if (scalar @{$read_distribution} == 0) {
		return \@trains;
	}

	my %current_train = (begin => $read_distribution->[0]{'begin'},
		end => $read_distribution->[0]{'end'}, read_count => 0, forward_read_count => 0,
		spikes => [], classifier => KMean->new());
	$current_train{'classifier'}->add_point(0);

	my $last_read_count = 0;
	my $spike_detection = 2;
	my $total_read_count = 0;
	my @end_reads = ();

	foreach my $read_locus (@{$read_distribution}) {
		my $position = $read_locus->{'begin'};
		if ($current_train{'begin'} <= $position && $position < $current_train{'end'}) {
			$current_train{'end'} = max($current_train{'end'}, $read_locus->{'end'});
			$current_train{'read_count'} += $read_locus->{'read_count'};
			$current_train{'forward_read_count'} += $read_locus->{'forward_read_count'};
			static__get_trains__process_train_spikes(\%current_train, $position, $read_locus, \$total_read_count, \@end_reads,
				\$last_read_count, $spike_detection);
		}
		else { # the read train ended
			static__get_trains__maintain_read_count(\%current_train, $position, \$total_read_count, \$last_read_count, \@end_reads,
				$spike_detection);
			static__get_trains__finish_train_spikes(\%current_train);
			static__get_trains__add_candidate_train(\@trains, \%current_train, $train_detection_threshold);
			$current_train{'begin'} = $position;
			$current_train{'end'} = $read_locus->{'end'};
			$current_train{'read_count'} = $read_locus->{'read_count'};
			$current_train{'forward_read_count'} = $read_locus->{'forward_read_count'};
			$current_train{spikes} = [];
			$total_read_count = 0;
			@end_reads = ();
			$last_read_count = 0;
			static__get_trains__process_train_spikes(\%current_train, $position, $read_locus, \$total_read_count, \@end_reads,
				\$last_read_count, $spike_detection);
		}
		$last_read_count = $total_read_count;
	}
	static__get_trains__finish_train_spikes(\%current_train);
	static__get_trains__add_candidate_train(\@trains, \%current_train, $train_detection_threshold);
	return \@trains;
}


=method __add_candidate_window_from_train

Private helper function. You shouldnt use this function.

=cut
sub static__get_trains__add_candidate_train {
	my ($trains, $current_train, $train_detection_threshold) = @_;
	my $reverse_read_count = $current_train->{'read_count'} - $current_train->{'forward_read_count'};
	if ($current_train->{'forward_read_count'} >= $train_detection_threshold || $reverse_read_count >= $train_detection_threshold) {
		push @{$trains}, {%{$current_train}, 'strand' => get_strand($current_train->{'forward_read_count'}, $current_train->{'read_count'})};
	}
    return;
}

=method static__get_trains__maintain_read_count

Static private helper function. You shouldnt use this function.

=cut
sub static__get_trains__maintain_read_count {
	my ($current_train, $position, $total_read_count, $last_read_count, $end_reads, $spike_detection) = @_;
	my $spikes = $current_train->{spikes};

	# Maintaining read count
	if (scalar @{$spikes} && $spikes->[-1]{'end'} == -1) {
		my $last_spike = $spikes->[-1];
		while (scalar @{$end_reads} && $end_reads->[0]{'end'} <= $position) {
			${$last_read_count} = ${$total_read_count};
			${$total_read_count} -= $end_reads->[0]{'read_count'};
			my $trigger = $end_reads->[0]{'end'};
			shift @{$end_reads};
			my $added = $current_train->{'classifier'}->add_point(${$total_read_count});
			if ($added == KMean::ASSIGNED_FIRST_CLASS) {
				$last_spike->{'end'} = $trigger;
				$last_spike->{'strand'} = get_strand $last_spike->{'forward_read_count'}, $last_spike->{'read_count'};
				$current_train->{'classifier'}->clear_points();
				$current_train->{'classifier'}->add_point(${$total_read_count});
				last;
			}
		}
	}
	while (scalar @{$end_reads} && $end_reads->[0]{'end'} <= $position) {
		${$last_read_count} = ${$total_read_count};
		${$total_read_count} -= $end_reads->[0]{'read_count'};
		shift @{$end_reads};
	}
    return;
}


=method static__get_trains__finish_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub static__get_trains__finish_train_spikes {
	my ($current_train) = @_;
	my $spikes = $current_train->{spikes};
	if (scalar @{$spikes} && $spikes->[-1]{'end'} == -1) {
		my $last_spike = $spikes->[-1];
		$last_spike->{'end'} = $current_train->{'end'};
		$last_spike->{'strand'} = get_strand $last_spike->{'forward_read_count'}, $last_spike->{'read_count'};
	}
	elsif (scalar @{$spikes} == 0 && $current_train->{'end'} - $current_train->{'begin'} <= 40) {
		my $strand = get_strand $current_train->{'forward_read_count'}, $current_train->{'read_count'};
		push @{$spikes}, {begin => $current_train->{'begin'}, end => $current_train->{'end'}, trigger => 0,
		read_count => $current_train->{'read_count'}, forward_read_count => $current_train->{'forward_read_count'}, strand => $strand};
	}
	$current_train->{'classifier'}->clear_points();
	$current_train->{'classifier'}->add_point(0);
    return;
}


=method static__get_trains__process_train_spikes

Static private helper function. You shouldnt use this function.

=cut
sub static__get_trains__process_train_spikes {
	my ($current_train, $position, $read_locus, $total_read_count, $end_reads, $last_read_count, $spike_detection) = @_;
	my $spikes = $current_train->{spikes};

	# Maintaining read count
	static__get_trains__maintain_read_count($current_train, $position, $total_read_count, $last_read_count, $end_reads, $spike_detection);
	${$total_read_count} += $read_locus->{'read_count'};
	push @{$end_reads}, {end => $read_locus->{'end'}, read_count => $read_locus->{'read_count'}};
	@{$end_reads} = sort {$a->{'end'} <=> $b->{'end'}} @{$end_reads};
	# Looking for spikes
	my $added = $current_train->{'classifier'}->add_point(${$total_read_count});
	if ($added == KMean::ASSIGNED_SECOND_CLASS) {
		if (scalar @{$spikes}) {
			my $last_spike = $spikes->[-1];
			if ($last_spike->{'end'} == -1) {
				if ($current_train->{'classifier'}->class_of($last_spike->{'trigger'}) == KMean::ASSIGNED_FIRST_CLASS) {
					$last_spike->{'begin'} = $position;
					$last_spike->{'trigger'} = ${$total_read_count};
					$last_spike->{'read_count'} = $read_locus->{'read_count'};
					$last_spike->{'forward_read_count'} = $read_locus->{'forward_read_count'};
				}
				return;
			}
		}
		push @{$spikes}, {begin => $position, end => -1, trigger => ${$total_read_count}, read_count => $read_locus->{'read_count'},
		forward_read_count => $read_locus->{'forward_read_count'}};
	}
	elsif ($added == KMean::ASSIGNED_FIRST_CLASS) {
		my $last_spike = $spikes->[-1];
		if ($last_spike->{'end'} == -1) {
			$last_spike->{'end'} = $position;
			$last_spike->{'strand'} = get_strand $last_spike->{'forward_read_count'}, $last_spike->{'read_count'};
		}
		$current_train->{'classifier'}->clear_points();
		$current_train->{'classifier'}->add_point(${$total_read_count});
	}
	elsif (scalar @{$spikes} && $spikes->[-1]{'end'} == -1) {
		my $last_spike = $spikes->[-1];
		$last_spike->{'read_count'} += $read_locus->{'read_count'};
		$last_spike->{'forward_read_count'} += $read_locus->{'forward_read_count'};
	}
    return;
}

1;
