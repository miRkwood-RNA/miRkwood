# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use Time::HiRes ();

use Clusters;
use Storable;
use File::Basename;

# ../data/ShortStack/ShortStack-2.0.9/SRR

my @runs = ({bam => "../data/shortstack_miRNA.bam", fa => "../data/Athaliana_167.fa", out => "output_miRNA_trains"},
{bam => "../data/shortstack_others.bam", fa => "../data/Athaliana_167.fa", out => "output_others_trains"},
{bam => "../data/ShortStack/ShortStack-2.0.9/SRR051927.bam", fa => "../data/ShortStack/Arabidopsis_thaliana.TAIR10.23.dna.toplevel.fa", out => "SRR051927"},
{bam => "../data/ShortStack/ShortStack-2.0.9/SRR275588.bam", fa => "../data/ShortStack/Arabidopsis_thaliana.TAIR10.23.dna.toplevel.fa", out => "SRR275588"},
{bam => "../data/ShortStack/ShortStack-2.0.9/SRR032102.bam", fa => "../data/ShortStack/Oryza_sativa.IRGSP-1.0.23.dna.toplevel.fa", out => "SRR032102"},
{bam => "../data/ShortStack/ShortStack-2.0.9/SRR488774.bam", fa => "../data/ShortStack/Zea_mays.AGPv3.23.dna.toplevel.fa", out => "SRR488774"});

# foreach my $current_run (@runs) {

	my $current_run = $runs[0];

	my $clustering = Clusters->new($current_run->{bam}, $current_run->{fa}, 2);
	print "Getting the reads...";
	my $start = Time::HiRes::gettimeofday();
	my $read_loci;
	my $read_loci_filename = 'cache/perl_read_distrib__' . basename($clustering->{bam_file}). '.dat';
	if (-e $read_loci_filename) {
		$read_loci = retrieve($read_loci_filename);
	}
	else {
		$read_loci = $clustering->get_read_distribution_per_chr();
		store $read_loci, $read_loci_filename;
	}
	# my $read_distrib = $clustering->compute_read_distribution_per_chr_from_read_loci($reads);
	my $end_reads = Time::HiRes::gettimeofday();

	# begin
	# spaces
	# my $spaces = $clustering->compute_read_spaces_chr_free($reads);
	# my $distrib = $clustering->get_read_spaces_distribution_chr_free($spaces);
	# $clustering->plot_space_distribution_chr_free($distrib);

	#density
	# my $densities = $clustering->compute_read_density_distribution($reads, 10);
	# $clustering->plot_read_density_distribution($densities, 10);
	# die;

# 	trains
# 	my $reads = $clustering->get_reads_per_chr();
# 	my $read_trains = $clustering->compute_read_train_per_chr($reads);
# 	my $read_train_distributions = $clustering->compute_train_distributions_per_chr($read_trains);
# 	$clustering->plot_read_train_distribution($read_train_distributions, $current_run->{out});
# 	die;
# 	end

	print " Done. (in ", $end_reads-$start, "s)\n";
	# Computing parameters...";
	# # my $params;
	# # my $params = $clustering->compute_poisson_parameters($reads);
	# my $end_params = Time::HiRes::gettimeofday();
	# print " Done. (in ", $end_params-$end_reads,"s)\n";

	# print "$chr\t$params->[0]\t$params->[1]\\n";

	my @thresholds = (6);

	foreach my $threshold (@thresholds) {
		$clustering->{threshold} = $threshold;
		print "Computing windows (threshold=$threshold)...";
		my $end_params = Time::HiRes::gettimeofday();
		my $windows = $clustering->get_windows_from_train_analysis_with_read_distribution($read_loci, 2);
		my $miRnaPos = $clustering->process_window_spikes($windows, $threshold);
		my $window_dists = $clustering->compute_window_length_distribution($windows);
		mkdir $current_run->{out};
		$clustering->plot_window_distribution($window_dists, $current_run->{out});
		$clustering->export_windows_to_gff($windows, $current_run->{out});
		$clustering->export_miRnaPos_to_gff($miRnaPos, $current_run->{out}, '../data/TAIR10_miRNA_only_covered_by_reads.gff3');
		my $end_windows = Time::HiRes::gettimeofday();
		print " Done. (in ", $end_windows-$end_params,"s)\n";
	}

# }

