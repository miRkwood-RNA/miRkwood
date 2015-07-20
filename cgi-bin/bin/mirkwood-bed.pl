#!/usr/bin/perl

# PODNAME: mirkwood-bed.pl
# ABSTRACT: miRkwood - A micro-RNA analysis pipeline for sRNAseq analysis

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
    'shuffles'           => \$randfold,
    'no-mfei'            => \$no_mfei,
    'align'              => \$align,
    'no-filter_otherRNA' => \$no_filter_tRNA_rRNA,
    'no-filter_CDS'      => \$no_filter_CDS,
    'no-filter_multimapped' => \$no_filter_multimapped,
    'no-varna'           => \$no_varna,
    'output=s'           => \$output_folder,
    'genome=s'           => \$genome_file,
    'help|?'             => \$help,
    'force'              => \$force,
    man                  => \$man
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

if (! -d $output_folder){
	mkdir $output_folder, 0777;
}

if( my @files = glob("$output_folder/*") ) {
    if ( $force ){
        print "Directory $output_folder is not empty. It will be cleared out.\n";
        system("rm -Rf $output_folder");
        mkdir $output_folder, 0777;
    }
    else{
        die("Directory $output_folder is not empty. Please clear it out or choose another directory.");
    }
} 

my $abs_output_folder = File::Spec->rel2abs($output_folder);
if ( !-e $abs_output_folder ) {
    mkdir $output_folder or die("Error when creating $abs_output_folder");
}

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
                                         'WebBAM',
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
my $results_folder       = miRkwood::Paths::create_folder( File::Spec->catdir( $abs_output_folder, 'Results' ) );
my $new_results_folder   = miRkwood::Paths::create_folder( File::Spec->catdir( $results_folder, 'new_miRNAs' ) );
my $known_results_folder = miRkwood::Paths::create_folder( File::Spec->catdir( $results_folder, 'known_miRNAs' ) );

my $new_candidates_dir   = miRkwood::Paths::get_new_candidates_dir_from_job_dir( $abs_output_folder );
my $known_candidates_dir = miRkwood::Paths::get_known_candidates_dir_from_job_dir( $abs_output_folder );

miRkwood::CLI::process_results_dir_for_offline($new_candidates_dir, $new_results_folder, 'smallRNAseq', 'New');
miRkwood::CLI::process_results_dir_for_offline($known_candidates_dir, $known_results_folder, 'smallRNAseq', 'Known');
