# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use Time::HiRes ();

use Clusters;
use Storable;
use File::Basename;

use miRkwood::Programs;

miRkwood::Programs::init_programs();

# ../data/ShortStack/ShortStack-2.0.9/SRR

my @runs = (
{reads => "../data/shortstack_miRNA.bam", bam => 1, fa => "../data/Athaliana_167.fa", out => "shortstack_miRNA", gff => '../data/TAIR10_miRNA.gff3',
extended_bed => "../data/shortstack_miRNA.bed"},
{reads => "../data/Bed_isabelle/SRR051927_filtered.bed", fa => "../data/Athaliana_167.fa", out => "SRR051927", gff => '../data/TAIR10_miRNA.gff3',
extended_bed => "../data/Bed_isabelle/SRR051927_filtered_extended.bed"},
{reads => "../data/Bed_isabelle/SRR275588_filtered.bed", fa => "../data/Athaliana_167.fa", out => "SRR275588",  gff => '../data/TAIR10_miRNA.gff3',
extended_bed => "../data/Bed_isabelle/SRR275588_filtered_extended.bed"},
{reads => "../data/Bed_isabelle/SRR065156_filtered.bed", fa => "../data/Bed_isabelle/Genomes/Solanum_lycopersicum_2.4.fasta", out => "SRR065156",
extended_bed => "../data/Bed_isabelle/SRR065156_filtered_extended.bed"},
{reads => "../data/Bed_isabelle/SRR032102_filtered.bed", fa => "../data/Bed_isabelle/Genomes/Osativa_204_v7.0.formated_2.fa", out => "SRR032102",
extended_bed => "../data/Bed_isabelle/SRR032102_filtered_extended.bed"}
# {reads => "../data/Bed_isabelle/SRR488774.bed", fa => "../data/ShortStack/Zea_mays.AGPv3.23.dna.toplevel.fa", out => "SRR488774"}
);

foreach my $current_run (@runs) {

	my $clustering = Clusters->new($current_run->{fa});
	print "Running $current_run->{out}...\n";
	mkdir $current_run->{out};
	mkdir $current_run->{out}. '/rnalfold_cache';
	mkdir $current_run->{out}. '/rnastemloop_cache';
	print "Getting the reads...";
	my $start = Time::HiRes::gettimeofday();
	my $read_distrib;
	my $read_distrib_filename = 'cache/perl_read_distrib__' . basename($current_run->{reads}). '.dat';
	if (-e $read_distrib_filename) {
		$read_distrib = retrieve($read_distrib_filename);
	}
	else {
		if (defined $current_run->{bam} && $current_run->{bam} == 1) {
			$read_distrib = $clustering->get_read_distribution_from_bam($current_run->{reads});
		}
		else {
			$read_distrib = $clustering->get_read_distribution_from_bed($current_run->{reads});
		}
		store $read_distrib, $read_distrib_filename;
	}
	
	my $end_reads = Time::HiRes::gettimeofday();

	print " Done. (in ", $end_reads-$start, "s)\n";

	my $annotation_filename = $current_run->{out} . '/annotation.gff3';
	if (-e $annotation_filename) {
		$current_run->{annotation} = $annotation_filename;
	}
	elsif (defined $current_run->{gff}) {
		system("intersectBed -a $current_run->{gff} -b $current_run->{reads} -wa -s -u > $annotation_filename");
	}
	else {
		$annotation_filename = undef;
	}

# 	my $region = $clustering->{genome_db}->seq('Chr1', 28435, 28825);
# 	my $rnafold_result = Clusters::run_rnalfold('Test', $region, 'tmp');
# 	my $rnastemloop_result = Clusters::run_rnastemloop($rnafold_result, 'tmp/stemloop_result.txt', 'tmp/stemloop_opt_result.txt');
# 	print $rnastemloop_result, "\n";
# 	my $results = Clusters::parse_stemloop_output('Test', 28435, 'tmp/stemloop_opt_result.txt');
# 	print scalar @$results, " results found!\n";
# 	foreach my $result (@$results) {
# 		print "Begin = $result->{begin}, End = $result->{end}\n";
# 		print "Sequence : $result->{sequence}\n";
# 		print "Sequence : ", $clustering->{genome_db}->seq('Chr1', $result->{begin}, $result->{end}), "\n";
# 		print "Structure: $result->{secondary_structure}\n";
# 	}
# 	die;

	print "Computing windows...\n";
	my $end_params = Time::HiRes::gettimeofday();
	my $windows = $clustering->get_windows($read_distrib, 2);
	print "Processing windows...\n";
	my $miRnaPos = $clustering->process_window_spikes($windows);
	print "Creating candidate precursors...\n";
	my $regions = $clustering->compute_candidate_precursors_from_miRnaPos($miRnaPos);
	print "Running RNALfold and RNAstemloop...\n";
	my $new_regions = $clustering->apply_structure_criterion($regions, $current_run->{out});
	print "Exporting...";
	$clustering->export_windows_to_gff($windows, $current_run->{out});
	$clustering->export_miRnaPos_to_gff($miRnaPos, $current_run->{out}, $annotation_filename);
	$clustering->export_precursors_to_gff($regions, $current_run->{out}, 'precursors_list', $annotation_filename, $current_run->{extended_bed});
	$clustering->export_precursors_to_gff($new_regions, $current_run->{out}, 'precursors_list_with_RNAstemloop', $annotation_filename, $current_run->{extended_bed});
	my $end_windows = Time::HiRes::gettimeofday();
	
	print "Done. (in ", $end_windows-$end_params,"s)\n";
	
}

