#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Clusters');
}
require_ok('PipelineMiRNA::Clusters');

my $DEFAULT_mindepth = 20;
my $DEFAULT_pad      = 100;

my $bamfile = input_file('Clusters.reads-Athaliana_167-ChrC.bam');
file_exists_ok($bamfile);

my $genome_file = input_file('Clusters.Athaliana_167-ChrC.fa');
file_exists_ok($genome_file);

## get_faidx_file() ##

ok( my $faidx_file = PipelineMiRNA::Clusters->get_faidx_file($genome_file),
    'Can get existing FAIDX file' );
file_exists_ok($faidx_file);

my $dummy_genome_file = input_file('dummy_genome.fa');
link( $genome_file, $dummy_genome_file );
file_exists_ok($dummy_genome_file);

ok(
    my $dummy_faidx_file =
      PipelineMiRNA::Clusters->get_faidx_file($dummy_genome_file),
    'Can get unexisting FAIDX file'
);
file_exists_ok($dummy_faidx_file);
my $contents1 = slurp_file($faidx_file);
my $contents2 = slurp_file($dummy_faidx_file);
is( $contents1, $contents2, 'get_faidx_file created the correct index file' );
unlink ($dummy_faidx_file, $dummy_genome_file);

## get_islands() ##

ok(
    my @get_islands_output =
      PipelineMiRNA::Clusters->get_islands( $bamfile, $DEFAULT_mindepth, $faidx_file ),
    'Can call get_islands()'
);
my @get_islands_expected =
  ( ['ChrC', 6, 77], ['ChrC', 292, 318], ['ChrC', 448, 478], ['ChrC', 487, 515] );
is_deeply( \@get_islands_output, \@get_islands_expected,
    'get_islands returns the correct values' );

## merge_clusters() ##

ok(
    my @merge_clusters_output = PipelineMiRNA::Clusters->merge_clusters(
        \@get_islands_output, $DEFAULT_pad, $genome_file
    ),
    'Can call merge_clusters()'
);
my @merge_clusters_expected = ( ['ChrC', 6, 77], ['ChrC', 292, 515] );
is_deeply( \@merge_clusters_output, \@merge_clusters_expected,
    'merge_clusters returns the correct values' );

## get_clusters() ##

ok(
    my @get_clusters_output =
      PipelineMiRNA::Clusters->get_clusters(
        $bamfile, $genome_file, $DEFAULT_mindepth, $DEFAULT_pad
      ),
    'Can call get_clusters()'
);
is_deeply(\@get_clusters_output, \@merge_clusters_expected,
    'get_clusters returns the correct values' );

## get_sequences_from_clusters() ##

my @get_sequences_from_clusters_expected1 = (
    'ChrC__6-77',
'GAACGACGGGAATTGAACCCGCGATGGTGAATTCACAATCCACTGCCTTAATCCACTTGGCTACATCCGCCCCTACG'
);
my @get_sequences_from_clusters_expected2 = (
    'ChrC__292-515',
'AAAACGAGTAATATAAGCCCTCTTTCTTATTTAAAGAAGGCTTATATTGCTCGTTTTTTACTAAACTAGATCTAGACTAACACTAACGAATTATCCATTTGTAGATGGAGCCTCAACAGCAGCTAGGTCTAGAGGGAAGTTGTGAGCATTACGTTCATGCATAACTTCCATACCAAGGTTAGCACGGTTAATAATATCAGCCCAAGTATTAATAACACGTCCTTGACTATCAACTACTGATTGGTTGAAATTGAAACCATTTAGGTTGAAAGCCATAGTACTAATACCTAAAGCAGTAAACCAAATACCTACTACCGGCCAAGCCGCTAAGAAGAAATGTAAAGAACGAGAATTGTTGAAACTAGCATATTGGAAAATCAATCGGCCAAAATAACCGTGAGCAGCTACAATGTTGTAAGTTTCTTCTTCTTGCCCGAATCTGTAACCTTCATTAGCAGATTCATTTTCTGTGGTTTCCCTGATCAAACTAGAAGTTACCAAGGAACCATGCATAG'
);
my @get_sequences_from_clusters_expected = (
    \@get_sequences_from_clusters_expected1,
    \@get_sequences_from_clusters_expected2
);
ok(
    my @get_sequences_from_clusters_output =
      PipelineMiRNA::Clusters->get_sequences_from_clusters(
        $genome_file, \@merge_clusters_output
      ),
    'Can call get_sequences_from_clusters()'
);
is_deeply(
    \@get_sequences_from_clusters_output,
    \@get_sequences_from_clusters_expected,
    'get_sequences_from_clusters returns the correct values'
);

## get_clustered_sequences_from_bam() ##

ok(
    my @get_clustered_sequences_from_bam_output =
      PipelineMiRNA::Clusters->get_clustered_sequences_from_bam(
        $bamfile, $genome_file, $DEFAULT_mindepth, $DEFAULT_pad
      ),
    'Can call get_clustered_sequences_from_bam()'
);
is_deeply(
    \@get_clustered_sequences_from_bam_output,
    \@get_sequences_from_clusters_expected,
    'get_clustered_sequences_from_bam returns the correct values'
);
