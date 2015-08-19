#!/usr/bin/perl

# PODNAME: mirkwood-bed.pl
# ABSTRACT: miRkwood - A micro-RNA analysis pipeline for sRNAseq analysis

use lib '/vagrant/cgi-bin/lib/';

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
my $genome_file;
my $varna = 1;
my $no_varna = 0;
my $filter_tRNA_rRNA = 1;
my $no_filter_tRNA_rRNA = 0;
my $filter_CDS = 1;
my $no_filter_CDS = 0;
my $filter_multimapped = 1;
my $no_filter_multimapped = 0;
my $randfold;
my $mfei = 1;
my $no_mfei = 0;
my $align = 0;
my $species = '';
my $species_db = '';
my $force = 0;

##### Get options
GetOptions(
    'genome=s'              => \$genome_file,
    'output=s'              => \$output_folder,
    'shuffles'              => \$randfold,
    'align'                 => \$align,
    'no-filter-mfei'        => \$no_mfei,
    'no-filter-CDS'         => \$no_filter_CDS,
    'no-filter-t-r-snoRNA'  => \$no_filter_tRNA_rRNA,
    'no-filter-multimapped' => \$no_filter_multimapped,
    'no-varna'              => \$no_varna,
    'help|?'                => \$help,
    'force'                 => \$force,
    man                     => \$man
) || pod2usage( -verbose => 0 );
pod2usage( -verbose => 1 ) if ($help);
pod2usage( -verbose => 2 ) if ($man);

pod2usage("$0: No BED file given.") if ( @ARGV == 0 );
pod2usage("$0: No genome file given.") if ( ! $genome_file );

if ( $no_mfei ){
    $mfei = 0;
}
if ( $no_filter_tRNA_rRNA ){
    $filter_tRNA_rRNA = 0;
}
if ( $no_filter_CDS ){
    $filter_CDS = 0;
}
if ( $no_filter_multimapped ){
    $filter_multimapped = 0;
}

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


# Image
if ( $no_varna ){
    $varna = 0;
}

# Check input files
my $bed_file = $ARGV[0];
( -e $bed_file ) or die("$bed_file is not a file");

( -e $genome_file ) or die("Genome file $genome_file is not a file");
if ( $genome_file =~ /([^.\/]+)[.](fa|fasta)/ ){
    $species = $1;
}

my $basename_bed = '';
if ( $bed_file =~ /.*\/([^\/.]+)[.]bed/ ){
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
                                         $species_db,
                                         $filter_CDS,
                                         $filter_tRNA_rRNA,
                                         $filter_multimapped,
                                         $mfei,
                                         $randfold,
                                         $varna);


##### Launch pipeline
my $pipeline = miRkwood::BEDPipeline->new($output_folder, $bed_file, $genome_file);
$pipeline->run_pipeline();


##### Write results
miRkwood::CLI::process_results_dir_for_offline($abs_output_folder, 'smallRNAseq', 'novel_miRNA');
miRkwood::CLI::process_results_dir_for_offline($abs_output_folder, 'smallRNAseq', 'known_miRNA');


__END__

=head1 SYNOPSIS

./mirkwood-bed.pl [options] [BED file]

=head1 OPTIONS

=over 8

=head2 Mandatory options

=item B<--genome>

Path to the genome (fasta format).

=item B<--output>

Output directory. If non existing it will be created. The directory must be empty.

=head2 Additional options

=item B<--shuffles>

Compute thermodynamic stability (shuffled sequences).

=item B<--align>

Flag conserved mature miRNAs (alignment with miRBase + miRdup).

=item B<--no-filter-mfei>

Don't filter out sequences with MFEI >= -0.6.
Default : only keep sequences with MFEI < -0.6.

=item B<--no-filter-CDS>

Don't filter out CDS.
Default: if an annotation GFF file is available CDS are filtered out.

=item B<--no-filter-t-r-snoRNA>

Don't filter out rRNA, tRNA, snoRNA.
Default: if an annotation GFF file is available rRNA, tRNA, snoRNA are filtered out.

=item B<--no-filter-multimapped>

Don't filter out multimapped reads.
Default: reads that map at more than 5 positions are filtered out.

=item B<--no-varna>

Disable the structure generation using Varna.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<miRkwood> is a micro-RNA analysis pipeline.

=cut
