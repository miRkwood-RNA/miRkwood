#!/usr/bin/perl -w

# PODNAME: mirkwood.pl
# ABSTRACT: miRkwood - A micro-RNA analysis pipeline

use lib '/vagrant/cgi-bin/lib/';

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
my $varna        = 0;
my $no_process   = 0;
my $mask         = 0;
my $trna         = 0;
my $rrna         = 0;
my $output_folder = '';
my $fasta_file = '';

## Parse options
GetOptions(
    'shuffles'       => \$shuffles,
    'input=s'        => \$fasta_file,
    'filter-mfei'    => \$mfei,
    'align'          => \$align,
    'both-strands'   => \$both_strands,
    'varna'          => \$varna,
    'no-process'     => \$no_process,
    'species-mask=s' => \$species_mask,
    'filter-trna'    => \$trna,
    'filter-rrna'    => \$rrna,
    'output=s'       => \$output_folder,
    'help|?'         => \$help,
    'man'            => \$man
) || pod2usage( -verbose => 0 );
pod2usage( -verbose => 1 ) if ($help);
pod2usage( -verbose => 2 ) if ($man);

pod2usage("$0: No FASTA file given.") if ( $fasta_file eq '' );

if ($output_folder eq ''){
    die("You must indicate an empty directory with the --output option.");
}

$output_folder = miRkwood::Paths::create_folder( $output_folder );

if( my @files = glob("$output_folder/*") ) {
     die("Directory $output_folder is not empty. Please clear it out or choose another directory.");
}   

if ($species_mask) {
    $mask = 1;
}

( -e $fasta_file ) or die("$fasta_file is not a file");

my $abs_output_folder = miRkwood::Paths::create_folder( File::Spec->rel2abs($output_folder) );


# Importing modules after directory creation
use miRkwood;
use miRkwood::FastaPipeline;
use miRkwood::CLI;
use miRkwood::Paths;

my $seq_name = 'input_sequences.fas';
my $seq_path = File::Spec->catfile( $abs_output_folder, $seq_name );

File::Copy::copy( $fasta_file, $seq_path );

my $run_options_file =
  miRkwood::Paths->get_job_config_path($abs_output_folder);
miRkwood->CONFIG_FILE($run_options_file);
miRkwood::write_config( $run_options_file, $both_strands, $mask, $trna, $rrna, $mfei, $shuffles,
    $align, '', $species_mask, $varna, 'abinitio', $abs_output_folder );

my $pipeline = miRkwood::FastaPipeline->new($abs_output_folder);
$pipeline->run_pipeline();


unless ($no_process) {
    my $results_folder = miRkwood::Paths::create_folder( File::Spec->catdir( $abs_output_folder, 'results' ) );
    my $candidates_dir = miRkwood::Paths::get_dir_candidates_path_from_job_dir( $abs_output_folder );
    miRkwood::CLI::process_results_dir_for_offline($abs_output_folder, 'abinitio') unless $no_process;
}

__END__

=head1 SYNOPSIS

./mirkwood.pl [options]

=head1 OPTIONS

=over 8

=item B<--input>

Path to the fasta file.

=item B<--output>

Output directory. If non existing it will be created. The directory must be empty.

=item B<--both-strands>

Scan both strands.

=item B<--species-mask>

Mask coding regions against the given organism.

=item B<--shuffles>

Compute thermodynamic stability (shuffled sequences).

=item B<--filter-mfei>

Select only sequences with MFEI < -0.6.

=item B<--filter-rrna>

Filter out ribosomal RNAs (using RNAmmer).

=item B<--filter-trna>

Filter out tRNAs (using tRNAscan-SE).

=item B<--align>

Flag conserved mature miRNAs (alignment with miRBase + miRdup).

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
