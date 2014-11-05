use strict;
use warnings;

use MiRnaDuplexDetector;
use Storable;

use Time::HiRes ();

# my $w_len = 12;
# my $mismatch = 2;

# mtr-MIR166d
# my $non_detected_hairpin = 'GUUGA' . 'GGGGAAUGUUGCCUGGCUCAAGGCUUU' . 'UCAUAGUUUGAUCUAUCUACUAUGAUUAAGGCUUCGGGCCAGGCUUCAUCCCCCUCAAA'; # Starting at 5
# my $non_detected_miRNAa = 'UCGGGCCAGGCUUCAUCCCCC';
# my $read = $non_detected_miRNAa, my $seq = $non_detected_hairpin;

# my %final_result = %{retrieve '../data/mirbase/matching_results.dat'};
# my $final_result_10_2 = $final_result{12}{2};
# my %miRNAs = %{retrieve('../data/mirbase/miRNAs.dat')};
# my $total = keys %miRNAs;
# 
# foreach my $miRNA (keys %miRNAs) {
# 	my $results = $final_result_10_2->{$miRNA};
# 	if (scalar(@$results)) {}
# 	else {
# 		print $miRNA, "\n";
# 	}
# }
# 
# my %hairpins = %{retrieve('../data/mirbase/hairpins.dat')};
# my %miRNAs = %{retrieve('../data/mirbase/miRNAs.dat')};
# my $total = keys %miRNAs;
# # 
# # my @w_lens = (15..24);
# my @mismatches = (1..7);
# # 
# # my %final_result = ();
# # 
# # foreach my $w_len (@w_lens) {
# # # 	$final_result{$w_len} = {};
# 	foreach my $mismatch (@mismatches) {
# 		my $found = 0;
# 		foreach my $miRNA (keys %miRNAs) {
# 			my $results = MiRnaDuplexDetector::local_align_forward_strand($miRNAs{$miRNA}[0], $hairpins{$miRNAs{$miRNA}[1]}, $mismatch); #$final_result{$miRNA};
# 			if (scalar(@$results)) {
# 				$found++;
# 			}
# # 			elsif ($mismatch == 6) {
# # 				print $miRNA, "\n";
# # 			}
# 		}
# 		print "For allowed mismatches = $mismatch, ", $found*100/$total, "% of miRNAs detected.\n";
# # 		print "For window size = $w_len, allowed mismatches = $mismatch, ", $found*100/$total, "% of miRNAs detected.\n";
# 	}
# # }


# 		$final_result{$w_len}{$mismatch} = {};
# 		my $found = 0;
# 		foreach my $miRNA (keys %miRNAs) {
# 			my $results = MiRnaDuplexDetector::local_align_forward_strand($miRNAs{$miRNA}[0], $hairpins{$miRNAs{$miRNA}[1]});
# 			if ($results->[0] >= $w_len && $results->[1] <= $mismatch) {
# 				$found++;
# 			}
# 			$final_result{$miRNA} = $results;
# 			if (scalar(@$results)) {
# 				$found++;
# 			}
# 		}
# 		print $found, "\n";
# 		print "For window size = $w_len, allowed mismatches = $mismatch, ", $found*100/$total, "% of miRNAs detected.\n";
# 	}
# }

# store \%final_result, '../data/mirbase/matching_results.dat';

# my %miRNAs = ();
# my $id, my $miRNA_id;
# my @fields;
# open (my $fh, '<', '../data/mirbase/hairpin.fa');
# while (my $line = <$fh>) {
# 	chomp $line;
# 	if ($line =~ /^>/) {
# 		@fields = split(' ', substr($line, 1));
# 		$id = lc $fields[0];
# 		$hairpins{$id} = '';
# 		next;
# 	}
# 	$hairpins{$id} .= $line;
# }
# 
# open (my $fh2, '<', '../data/mirbase/mature.fa');
# while (my $line = <$fh2>) {
# 	chomp $line;
# 	if ($line =~ /^>/) {
# 		@fields = split(' ', substr($line, 1));
# 		$miRNA_id = lc $fields[0];
# 		$id = $miRNA_id;
# 		$id =~ s/-[35]p$//g;
# 		if (!defined $hairpins{$id}) {
# 			$miRNA_id = undef;
# 			next;
# 		}
# 		$miRNAs{$miRNA_id} = ['', $id];
# 		next;
# 	}
# 	if ($miRNA_id) {
# 		$miRNAs{$miRNA_id}->[0] .= $line;
# 	}
# }
# 
# store \%hairpins, '../data/mirbase/hairpins.dat';
# store \%miRNAs, '../data/mirbase/miRNAs.dat';

my $total_gen = 0, my $total_run = 0, my $total_run_optimized = 0, my $total_results = 0, my $total_mismatches = 0;
for (my $i = 0; $i < 10000; $i++) {
	my $start = Time::HiRes::gettimeofday();
	my ($read, $seq) = MiRnaDuplexDetector::generate_candidates(22, 300);

	my $end_reads = Time::HiRes::gettimeofday();
	my $results = MiRnaDuplexDetector::local_align_forward_strand($read, $seq, 5);

	my $end_run = Time::HiRes::gettimeofday();
# 	my $results_opt = MiRnaDuplexDetector::run_optimized_cpp($read, $seq, $w_len, $mismatch);
#
# 	my $end_run_optimized = Time::HiRes::gettimeofday();
	$total_gen += $end_reads - $start;
	$total_run += $end_run - $end_reads;
# 	$total_run_optimized += $end_run_optimized - $end_run;
	if (scalar @$results) {
		$total_results++;
	}
}

print "Total gen time: $total_gen", 's', " Total run time: $total_run", "s Average match success: ", $total_results*100/10000, "\n";#, " Total run optimized: $total_run_optimized", "s.\n";

# my ($read, $seq) = MiRnaDuplexDetector::generate_candidates(24, 300);
# my $results = MiRnaDuplexDetector::run_reverse_strand($read, $seq, $w_len, $mismatch);
# if (scalar(@$results)) {
	# forward strand
# 	foreach my $result (@$results) {
# 		print "Read:", $result->[0]+$w_len, "-",$result->[0]," followed by Sequence:$result->[1]-", $result->[1]+$w_len, ":\n";
# 		my $read_rev = reverse substr($read, $result->[0], $w_len);
# 		print $read_rev, "\n";
# 		for (my $i = 0; $i < $w_len; $i++) {
# 			if (MiRnaDuplexDetector::__equiv(substr($read_rev, $i, 1), substr($seq, $result->[1]+$i, 1))) {
# 				print '|';
# 			}
# 			else {
# 				print ' ';
# 			}
# 		}
# 		print "\n". substr($seq, $result->[1], $w_len), "\n\n";
# 	}
	# reverse strand
# 	print "Read: $read\n";
# 	print "Read [-]: ",  MiRnaDuplexDetector::reverse_complement($read), "\n";
# 	print "Sequence: $seq\n";
# 	print "Sequence [-]: ",  MiRnaDuplexDetector::reverse_complement($seq), "\n";
# 	foreach my $result (@$results) {
# 		print "Read:", $result->[0], "-",$result->[0]+$w_len," followed by Sequence:", $result->[1]+$w_len, "-", $result->[1], ":\n";
# 		my $seq_rev_comp = MiRnaDuplexDetector::reverse_complement(substr($seq, $result->[1], $w_len));
# 		my $read_comp =  MiRnaDuplexDetector::complement(substr($read, $result->[0], $w_len));
# 		print $read_comp, "\n";
# 		for (my $i = 0; $i < $w_len; $i++) {
# 			if (MiRnaDuplexDetector::equiv_forward_strand(substr($read_comp, $i, 1), substr($seq_rev_comp, $i, 1))) {
# 				print '|';
# 			}
# 			else {
# 				print ' ';
# 			}
# 		}
# 		print "\n". $seq_rev_comp, "\n\n";
# 	}
# }
# else {
# 	print "No result.\n"
# }