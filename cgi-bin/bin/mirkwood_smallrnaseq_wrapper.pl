#!/usr/bin/perl
# PODNAME: mirkwood_smallrnaseq_wrapper.pl
# ABSTRACT: A wrapper for miRkwood - A micro-RNA analysis pipeline

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Archive::Zip;
use File::Spec;
use File::Basename;
use Cwd 'abs_path';

my $man  = 0;
my $help = 0;

### Wrapper options
my $log_file = '';
my $novel_results_html = '';
my $known_results_html = '';
my $zip_folder = '';

### miRkwood options
my $genome_file   = '';
my $bed_file      = '';
my $output_folder = '';
my $mirbase_file  = '';
my @annotation_gff;
my $shuffles = 0;
my $align    = 0;
my $no_mfei  = 0;
my $filter_bad_hairpins    = 1;
my $no_filter_bad_hairpins = 0;
my $min_read_positions_nb  = 0;
my $max_read_positions_nb  = 5;
my $varna = 0;
my $force = 0;


### Parse options
GetOptions(
    # wrapper options
    'log=s'          => \$log_file,
    'novel-html=s'   => \$novel_results_html,
    'known-html=s'   => \$known_results_html,
    'zip=s'          => \$zip_folder,
    # miRkwood options
    'genome=s'                => \$genome_file,
    'input=s'                 => \$bed_file,
    'output=s'                => \$output_folder,
    'shuffles'                => \$shuffles,
    'align'                   => \$align,
    'no-filter-mfei'          => \$no_mfei,
    'mirbase=s'               => \$mirbase_file,
    'gff=s'                   => \@annotation_gff,
    'no-filter-bad-hairpins'  => \$no_filter_bad_hairpins,
    'min-read-positions-nb=s' => \$min_read_positions_nb,
    'max-read-positions-nb=s' => \$max_read_positions_nb,
    'varna'                   => \$varna,
    'help|?'                  => \$help,
    'force'                   => \$force,
    'man'                     => \$man
) || pod2usage( -verbose => 0 );
pod2usage( -verbose => 1 ) if ($help);
pod2usage( -verbose => 2 ) if ($man);


### Check mandatory options
pod2usage("$0: No zip path given.") if ( $zip_folder eq '' );
pod2usage("$0: No log path given.") if ( $log_file eq '' );
pod2usage("$0: No html path given for novel miRNAs.") if ( $novel_results_html eq '' );
pod2usage("$0: No html path given for known miRNAs.") if ( $known_results_html eq '' );


### Run miRkwood
my $flag_options = '';
if ( $shuffles ){
    $flag_options .= ' --shuffles';
}
if ( $align ){
    $flag_options .= ' --align';
}
if ( $no_mfei ){
    $flag_options .= ' --no-filter-mfei';
}
if ( $no_filter_bad_hairpins ){
    $flag_options .= ' --no-filter-bad-hairpins';
}
if ( $varna ){
    $flag_options .= ' --varna';
}
if ( $help ){
    $flag_options .= ' --help';
}
if ( $force ){
    $flag_options .= ' --force';
}
if ( $man ){
    $flag_options .= ' --man';
}


my $mirkwood_path = File::Spec->catfile( abs_path(dirname(__FILE__)), "mirkwood-bed.pl" );
my $cmd = $mirkwood_path;
$cmd .= " --input $bed_file";
$cmd .= " --output $output_folder";
$cmd .= " --genome $genome_file";
if ( $mirbase_file ne '' and -r $mirbase_file ){
    $cmd .= " --mirbase $mirbase_file";
}
$cmd .= " --min-read-positions-nb $min_read_positions_nb";
$cmd .= " --max-read-positions-nb $max_read_positions_nb";
foreach my $gff_file (@annotation_gff){
    $cmd .= " -gff $gff_file";
}


$cmd .= $flag_options;
#~ print STDERR "$cmd\n";
system( $cmd );


### Copy log and html results files
system( "cp $output_folder/log.log $log_file" );
system( "cp $output_folder/results/novel_miRNA/results_novel_miRNA.html $novel_results_html" );
system( "cp $output_folder/results/known_miRNA/results_known_miRNA.html $known_results_html" );


### Create a zip with miRkwood output directory
my $zip = Archive::Zip->new();
$zip->addTree( $output_folder );
$zip->writeToFileNamed($zip_folder)



__END__

=head1 SYNOPSIS

perl mirkwood-bed.pl [options]

=head1 OPTIONS

=head2 Wrapper options

=over 8

=item B<--log>
Path to log file

=item B<--novel_html>
Path to HTML results file corresponding to novel miRNAs

=item B<--known_html>
Path to HTML results file corresponding to known miRNAs

=item B<--zip>
Path to zipped output directory

=back

=head2 miRkwood mandatory options

=over 8

=item B<--input>

Path to the BED file (created with our script mirkwood-bam2bed.pl).

=item B<--genome>

Path to the genome (fasta format).

=item B<--output>

Output directory. If non existing it will be created. The directory must be empty.

=back

=head2 miRkwood additional options

=over 8

=item B<--shuffles>

Compute thermodynamic stability (shuffled sequences).

=item B<--align>

Flag conserved mature miRNAs (alignment with miRBase + miRdup).

=item B<--no-filter-mfei>

Don't filter out sequences with MFEI >= -0.6.
Default : only keep sequences with MFEI < -0.6.

=item B<--mirbase>

If you have a gff file containing known miRNAs for this assembly,
use this option to give the path to this file.

=item B<--gff>

List of annotation files (gff or gff3 format).
Reads matching with an element of these files will be filtered out.
For instance you can filter out CDS by providing a suitable GFF file.

=item B<--no-filter-bad-hairpins>

By default the candidates with a quality score of 0 and no
conservation are discarded from results and are stored in a BED file.
Use this option to keep all results.

=item B<--min-read-positions-nb>

Minimum number of positions for each read to be kept.
Default : 0.

=item B<--max-read-positions-nb>

Maximum number of positions for each read to be kept.
Default : 5 (reads that map at more than 5 positions are filtered out).

=item B<--varna>

Allow the structure generation using Varna.

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<miRkwood> is a micro-RNA analysis pipeline.

=cut
