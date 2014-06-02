#!/usr/bin/perl -w

# PODNAME: mirkwood.pl
# ABSTRACT: miRkwood - A micro-RNA analysis pipeline for sRNAseq analysis

use warnings;
use strict;

use Pod::Usage;
use Getopt::Long;
use File::Copy;
use File::Spec;

my $man  = 0;
my $help = 0;

# Pipeline options
my $both_strands = 0;
my $shuffles     = 0;
my $mfei         = 0;
my $align        = 0;
my $species_mask = '';
my $genome_file  = '';
my $no_varna     = 0;
my $no_process   = 0;
my $mask          = 0;
my $output_folder = 'results_directory';

## Parse options
GetOptions(
    shuffles         => \$shuffles,
    mfei             => \$mfei,
    align            => \$align,
    'both-strands'   => \$both_strands,
    'no-varna'       => \$no_varna,
    'no-process'     => \$no_process,
    'species-mask=s' => \$species_mask,
    'output=s'       => \$output_folder,
    'genome=s'       => \$genome_file,
    'help|?'         => \$help,
    man              => \$man
) || pod2usage( -verbose => 0 );
pod2usage( -verbose => 1 ) if ($help);
pod2usage( -verbose => 2 ) if ($man);

pod2usage("$0: No BAM file given.") if ( @ARGV == 0 );
pod2usage("$0: No genome file given.") if ( ! $genome_file );

if ($species_mask) {
    $mask = 1;
}

my $varna = 1;
if ($no_varna) {
    $varna = 0;
}

my $bam_file = $ARGV[0];
( -e $bam_file ) or die("$bam_file is not a file");

( -e $genome_file ) or die("Genome file $genome_file is not a file");

my $abs_output_folder = File::Spec->rel2abs($output_folder);
if ( !-e $abs_output_folder ) {
    mkdir $output_folder or die("Error when creating $abs_output_folder");
}

# Importing modules after directory creation
use PipelineMiRNA;
use PipelineMiRNA::CLI;
use PipelineMiRNA::MainPipeline;
use PipelineMiRNA::Paths;
use PipelineMiRNA::Clusters;

my $run_options_file =
  PipelineMiRNA::Paths->get_job_config_path($abs_output_folder);
PipelineMiRNA->CONFIG_FILE($run_options_file);
PipelineMiRNA::write_config( $run_options_file, $both_strands, $mask, $mfei, $shuffles,
    $align, "", $species_mask, $varna, 'bam' );

PipelineMiRNA::MainPipeline::bam_pipeline($abs_output_folder, $bam_file, $genome_file);

unless ($no_process) {
	my $tmp_pieces_folder = File::Spec->catdir( $abs_output_folder, 'pieces' );
	if ( !-e $tmp_pieces_folder ) {
        mkdir $tmp_pieces_folder or die("Error when creating $tmp_pieces_folder");
    }
	PipelineMiRNA::CLI::process_results_dir_for_offline($abs_output_folder) unless $no_process;
}

__END__

=head1 SYNOPSIS

mirkwood [options] --genome GENOME [BAM file]

=head1 OPTIONS

=over 8

=item B<--genome>

The genome file to use

=item B<--both-strands>

Process both strands

=item B<--species-mask>

Mask coding regions against the given organism

=item B<--shuffles>

Compute thermodynamic stability

=item B<--mfei>

Compute MFE/MFEI/AMFE (minimal folding energy)

=item B<--align>

Align against mature microRNAs miRBase

=item B<--no-varna>

Disable the structure genration using Varna

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<miRkwood> is a micro-RNA analysis pipeline.

=cut
