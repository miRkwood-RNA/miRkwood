#!/usr/bin/perl

# PODNAME: mirkwood-bed.pl
# ABSTRACT: miRkwood - A micro-RNA analysis pipeline for sRNAseq analysis

my $script_directory;
my $lib_directory;

BEGIN {
    use File::Spec;
    use FindBin;
    $script_directory = $FindBin::Bin;
    $lib_directory = File::Spec->catdir( $script_directory, File::Spec->updir(), 'lib');
}

use lib '../lib/';
use lib $lib_directory;

use warnings;
use strict;

use Pod::Usage;
use Getopt::Long;
use File::Spec;

use miRkwood;
use miRkwood::Paths;
use miRkwood::BEDPipeline;
use miRkwood::CLI;


##### Variables
my $man  = 0;
my $help = 0;
my $job_title = 0;
my $output_folder = '';
my $bed_file = '';
my $genome_file;
my $varna = 0;
my $mirbase_file = '';
my @annotation_gff;
my $annotation_gff;
my $filter_bad_hairpins = 1;
my $no_filter_bad_hairpins = 0;
my $min_read_positions_nb = 0;
my $max_read_positions_nb = 5;
my $read_positions_nb_interval = '';
my $randfold;
my $mfei = 1;
my $no_mfei = 0;
my $align = 0;
my $species = '';
my $force = 0;

##### Get options
GetOptions(
    'genome=s'                => \$genome_file,
    'input=s'                 => \$bed_file,
    'output=s'                => \$output_folder,
    'shuffles'                => \$randfold,
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

pod2usage("$0: No BED file given.") if ( $bed_file eq '' );
pod2usage("$0: No genome file given.") if ( ! $genome_file );

$annotation_gff = join( '&', @annotation_gff);

if ( $no_mfei ){
    $mfei = 0;
}
if ( $no_filter_bad_hairpins ){
    $filter_bad_hairpins = 0;
}

if ( $min_read_positions_nb !~ /^-?(\d+)$/ ){
    my $msg = "ERROR: --min-read-positions-nb must be digit.\n";
    $msg .= "Enter 0 if you don't want to specify a minimum number of positions for each read.";
    die $msg;
}
if ( $max_read_positions_nb !~ /^-?(\d+)$/ ){
    my $msg = "ERROR: --max-read-positions-nb must be digit.\n";
    $msg .= "Enter 0 if you don't want to specify a maximum number of positions for each read.";
    die $msg;
}
if ( $min_read_positions_nb > $max_read_positions_nb ){
    die "ERROR: --min-read-positions-nb must be lower than --max-read-positions-nb.\n";
}
if ( $min_read_positions_nb < 0 || $max_read_positions_nb < 0 ){
    $min_read_positions_nb = 0;
    $max_read_positions_nb = 0;
}
$read_positions_nb_interval = "[$min_read_positions_nb;$max_read_positions_nb]";


# Check output folder
if ($output_folder eq ''){
    die("You must indicate an empty directory with the --output option.");
}

$output_folder = miRkwood::Paths::create_folder( $output_folder );

if( my @files = glob("$output_folder/*") ) {
    if ( $force ){
        print "Directory $output_folder is not empty. It will be cleared out.\n";
        system("rm -Rf $output_folder");
        $output_folder = miRkwood::Paths::create_folder( $output_folder );
    }
    else{
        die("Directory $output_folder is not empty. Please clear it out or choose another directory.");
    }
}

my $abs_output_folder = miRkwood::Paths::create_folder( File::Spec->rel2abs($output_folder) );


# Check input files
( -e $bed_file ) or die("$bed_file is not a file");

open (my $BED, '<', $bed_file) or die "Error when opening BED file: $!";
my $previous_position = '';
my $previous_chromosome = '';
while ( <$BED> ){
    my @fields = split( /\t/ );
    if ( $previous_chromosome ne '' && $fields[0] eq $previous_chromosome 
         && $previous_position ne '' && $fields[1] < $previous_position ){
        die 'Your BED file is not valid: reads are not sorted.';
    }
    $previous_chromosome = $fields[0];
    $previous_position = $fields[1];
    if ( ! miRkwood::Utils::is_correct_BED_line($_) ){
        die "Your BED file is not valid. Make sure you correctly used our provided script to convert BAM into BED file (file : $bed_file).\n";
    }
}
close $BED;

( -e $genome_file ) or die("Genome file $genome_file is not a file");
if ( $genome_file =~ /([^\/]+)[.](fa|fasta|dat)/ ){
    $species = $1;
}

my $basename_bed = '';
if ( $bed_file =~ /.*\/([^\/]+)[.](bed|dat)/ ){
    $basename_bed = $1;
}


##### Create config file
my $run_options_file = miRkwood::Paths->get_job_config_path($abs_output_folder);
miRkwood->CONFIG_FILE($run_options_file);
miRkwood::write_config_for_bam_pipeline( $run_options_file,
                                         $job_title,
                                         $species,
                                         'smallRNAseq',
                                         $basename_bed,
                                         $align,
                                         $mirbase_file,
                                         $annotation_gff,
                                         $filter_bad_hairpins,
                                         $read_positions_nb_interval,
                                         $mfei,
                                         $randfold,
                                         $varna,
                                         $abs_output_folder);


##### Check external softwares
miRkwood::Programs::init_programs();
my @unavailable = miRkwood::Programs::list_unavailable_programs();
if (@unavailable){
    my $error = "Cannot find required third-party software: @unavailable.";
    die($error);
}


##### Launch pipeline
my $pipeline = miRkwood::BEDPipeline->new($output_folder, $bed_file, $genome_file);
$pipeline->run_pipeline();


##### Write results
miRkwood::CLI::process_results_dir_for_offline($abs_output_folder, 'smallRNAseq', 'novel_miRNA');
miRkwood::CLI::process_results_dir_for_offline($abs_output_folder, 'smallRNAseq', 'known_miRNA');

my $summary_page = File::Spec->catfile( $abs_output_folder, 'summary.txt' );
miRkwood::Results::create_summary_page( $abs_output_folder, $summary_page );
#~ miRkwood::Results::clean_job_dir_for_cli_pipeline( $abs_output_folder );

__END__

=head1 SYNOPSIS

perl mirkwood-bed.pl [options]

=head1 OPTIONS

=head2 Mandatory options

=over 8

=item B<--input>

Path to the BED file (created with our script mirkwood-bam2bed.pl).

=item B<--genome>

Path to the genome (fasta format).

=item B<--output>

Output directory. If non existing it will be created. The directory must be empty.

=back

=head2 Additional options

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
