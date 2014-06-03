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

## Constructor ##
my @args = ($bamfile, $genome_file);
my $clustering_obj = new_ok('PipelineMiRNA::Clusters' => \@args );

## get_faidx_file() ##
can_ok( $clustering_obj, 'get_faidx_file' );
ok( my $faidx_file = $clustering_obj->get_faidx_file(),
    'Can get existing FAIDX file' );
file_exists_ok($faidx_file);

my $dummy_genome_file = input_file('dummy_genome.fa');
link( $genome_file, $dummy_genome_file );
file_exists_ok($dummy_genome_file);

my @dummy_args = ($bamfile, $dummy_genome_file);
my $dummy_clustering_obj = new_ok('PipelineMiRNA::Clusters' => \@dummy_args );
ok(
    my $dummy_faidx_file =
      $dummy_clustering_obj->get_faidx_file($dummy_genome_file),
    'Can get unexisting FAIDX file'
);
file_exists_ok($dummy_faidx_file);
my $contents1 = slurp_file($faidx_file);
my $contents2 = slurp_file($dummy_faidx_file);
is( $contents1, $contents2, 'get_faidx_file created the correct index file' );
unlink ($dummy_faidx_file, $dummy_genome_file);

## get_islands() ##
can_ok( $clustering_obj, 'get_islands' );
ok(
    my @get_islands_output =
      $clustering_obj->get_islands(),
    'Can call get_islands()'
);
my @get_islands_expected =
  ( ['ChrC', 6, 77], ['ChrC', 292, 318], ['ChrC', 448, 478], ['ChrC', 487, 515] );
is_deeply( \@get_islands_output, \@get_islands_expected,
    'get_islands returns the correct values' );

## merge_clusters() ##
can_ok( $clustering_obj, 'merge_clusters' );
ok(
    my @merge_clusters_output = $clustering_obj->merge_clusters(
        \@get_islands_output
    ),
    'Can call merge_clusters()'
);
my @merge_clusters_expected = ( ['ChrC', 6, 77], ['ChrC', 292, 515] );
is_deeply( \@merge_clusters_output, \@merge_clusters_expected,
    'merge_clusters returns the correct values' );

## get_clusters() ##
can_ok( $clustering_obj, 'get_clusters' );
ok(
    my @get_clusters_output =
      $clustering_obj->get_clusters(),
    'Can call get_clusters()'
);
is_deeply(\@get_clusters_output, \@merge_clusters_expected,
    'get_clusters returns the correct values' );

## get_sequences_from_clusters() ##
can_ok( $clustering_obj, 'get_sequences_from_clusters' );
my @get_sequences_from_clusters_expected1 = (
    'ChrC__1-192',
'TGGGCGAACGACGGGAATTGAACCCGCGATGGTGAATTCACAATCCACTGCCTTAATCCACTTGGCTACATCCGCCCCTACGCTACTATCTATTCTTTTTTGTATTGTCTAAAAAAAAAAAAAAATACAAATTTCAATAAAAAATAAAAAAAGGTAGCAAATTCCACCTTATTTTTTTTCTAATAAAAAATA'
);
my @get_sequences_from_clusters_expected2 = (
    'ChrC__254-554',
'TATGATACTCTATAAAAATTTGCTCATTTTTATAGAAAAAAACGAGTAATATAAGCCCTCTTTCTTATTTAAAGAAGGCTTATATTGCTCGTTTTTTACTAAACTAGATCTAGACTAACACTAACGAATTATCCATTTGTAGATGGAGCCTCAACAGCAGCTAGGTCTAGAGGGAAGTTGTGAGCATTACGTTCATGCATAACTTCCATACCAAGGTTAGCACGGTTAATAATATCAGCCCAAGTATTAATAACACGTCCTTGACTATCAACTACTGATTGGTTGAAATTGAAACCATTTAGGTTGAAAGCCATAGTACTAATACCTAAAGCAGTAAACCAAATACCTACTACCGGCCAAGCCGCTAAGAAGAAATGTAAAGAACGAGAATTGTTGAAACTAGCATATTGGAAAATCAATCGGCCAAAATAACCGTGAGCAGCTACAATGTTGTAAGTTTCTTCTTCTTGCCCGAATCTGTAACCTTCATTAGCAGATTCATTTTCTGTGGTTTCCCTGATCAAACTAGAAGTTACCAAGGAACCATGCATAGC'
);
my @get_sequences_from_clusters_expected = (
    \@get_sequences_from_clusters_expected1,
    \@get_sequences_from_clusters_expected2
);
ok(
    my @get_sequences_from_clusters_output =
      $clustering_obj->get_sequences_from_clusters( \@merge_clusters_output ),
    'Can call get_sequences_from_clusters()'
);
is_deeply(
    \@get_sequences_from_clusters_output,
    \@get_sequences_from_clusters_expected,
    'get_sequences_from_clusters returns the correct values'
);

## get_clustered_sequences_from_bam() ##
can_ok( $clustering_obj, 'get_clustered_sequences_from_bam' );
ok(
    my @get_clustered_sequences_from_bam_output =
      $clustering_obj->get_clustered_sequences_from_bam(
        $bamfile, $genome_file, $DEFAULT_mindepth, $DEFAULT_pad
      ),
    'Can call get_clustered_sequences_from_bam()'
);
is_deeply(
    \@get_clustered_sequences_from_bam_output,
    \@get_sequences_from_clusters_expected,
    'get_clustered_sequences_from_bam returns the correct values'
);

## extend_cluster() ##
can_ok( $clustering_obj, 'extend_cluster' );
my @extend_cluster_input1 = ( 'Chr1', 200, 250 );
ok(
    my @extend_cluster_output1 =
      $clustering_obj->extend_cluster( \@extend_cluster_input1 ),
    'Can call extend_cluster()'
);
my @extend_cluster_expected1 = ( 'Chr1', 75, 375 );
is_deeply( \@extend_cluster_output1, \@extend_cluster_expected1,
    'extend_cluster returns the correct values' );

## get_chromosomes_info_from_genome_file() ##
can_ok( $clustering_obj, 'get_chromosomes_info_from_genome_file' );
ok(
    my %chr_info_output =
      $clustering_obj->get_chromosomes_info_from_genome_file(
        $genome_file),
    'Can call get_chromosomes_info_from_genome_file()'
);
my %chr_info_expected = ( 'ChrC' => 153091 );
is_deeply( \%chr_info_output, \%chr_info_expected,
    'get_chromosomes_info_from_genome_file returns the correct values' );
