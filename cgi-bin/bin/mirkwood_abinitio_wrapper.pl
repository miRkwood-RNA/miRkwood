#!/usr/bin/perl
# PODNAME: mirkwood_abinitio_wrapper.pl
# ABSTRACT: A wrapper for miRkwood - A micro-RNA analysis pipeline

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Archive::Zip;

my $man  = 0;
my $help = 0;

### Wrapper options
my $log_file = '';
my $html_results_file = '';
my $zip_folder = '';

### miRkwood options
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
my $force = 0;


### Parse options
GetOptions(
    # wrapper options
    'log=s'          => \$log_file,
    'html=s'         => \$html_results_file,
    'zip=s'          => \$zip_folder,
    # miRkwood options
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
    'force'          => \$force,
    'man'            => \$man
) || pod2usage( -verbose => 0 );
pod2usage( -verbose => 1 ) if ($help);
pod2usage( -verbose => 2 ) if ($man);


### Check mandatory options
pod2usage("$0: No zip path given.") if ( $zip_folder eq '' );
pod2usage("$0: No log path given.") if ( $log_file eq '' );
pod2usage("$0: No html path given.") if ( $html_results_file eq '' );


### Run miRkwood
my $flag_options = '';
if ( $shuffles ){
    $flag_options .= ' --shuffles';
}
if ( $mfei ){
    $flag_options .= ' --filter-mfei';
}
if ( $align ){
    $flag_options .= ' --align';
}
if ( $both_strands ){
    $flag_options .= ' --both-strands';
}
if ( $varna ){
    $flag_options .= ' --varna';
}
if ( $no_process ){
    $flag_options .= ' --no-process';
}
if ( $trna ){
    $flag_options .= ' --filter-trna';
}
if ( $rrna ){
    $flag_options .= ' --filter-rrna';
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


my $mirkwood_path = './';
my $cmd = $mirkwood_path . "mirkwood.pl --input $fasta_file --output $output_folder $flag_options";
#~ print STDERR "$cmd\n";
system( $cmd );


### Copy log and html results files
system( "cp $output_folder/log.log $log_file" );
system( "cp $output_folder/results/results.html $html_results_file" );


### Create a zip with miRkwood output directory
my $zip = Archive::Zip->new();
$zip->addTree( $output_folder );
$zip->writeToFileNamed($zip_folder)



__END__

=head1 SYNOPSIS

perl mirkwood_abinitio_wrapper.pl [options]

=head1 OPTIONS

=head2 Wrapper options

=over 8

=item B<--log>
Path to log file

=item B<--html>
Path to HTML results file

=item B<--zip>
Path to zipped output directory

=back

=head2 miRkwood mandatory options

=over 8

=item B<--input>

Path to the fasta file.

=item B<--output>

Output directory. If non existing it will be created. The directory must be empty.

=back

=head2 miRkwood additional options

=over 8

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
