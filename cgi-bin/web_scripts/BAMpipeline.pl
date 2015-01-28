#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Spec;
use Log::Message::Simple qw[msg error debug];

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood;
use miRkwood::Paths;
use miRkwood::WebPaths;
use miRkwood::WebTemplate;
use miRkwood::Results;
use miRkwood::Utils;
use miRkwood::BEDHandler;
use miRkwood::BamPipelineSebastien;

my $error_url = miRkwood::WebTemplate::get_cgi_url('error.pl');


##### Create job id and job directory
my $jobId = miRkwood::Results->make_job_id( 'BAM' );
my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($jobId);
mkdir $absolute_job_dir;


##### Create log file
my $log_file = File::Spec->catfile( $absolute_job_dir, 'log.log' );
local $Log::Message::Simple::DEBUG_FH = miRkwood->LOGFH($log_file);


##### Check results directory
my $root = miRkwood::Paths->get_results_filesystem_path();

if (! -e $root) {
    my $error = "Designated directory ($root) for results does not exist. Please contact the system administrator";
    miRkwood::WebTemplate::web_die($error);
}

if (! -W $root) {
    my $error = "Cannot write results in designated directory $root. Please contact the system administrator";
    miRkwood::WebTemplate::web_die($error);
}


##### Check external softwares
miRkwood::Programs::init_programs();
my @unavailable = miRkwood::Programs::list_unavailable_programs();
if (@unavailable){
    my $error = "Cannot find required third-party software: @unavailable. Please contact the system administrator";
    miRkwood::WebTemplate::web_die($error);
}


##### Parameters
my $cgi        = CGI->new();
my $job_title  = $cgi->param('job');
my $mail       = $cgi->param('mail');
my $species    = $cgi->param('species');
my $species_db = $cgi->param('db');
my $filter_CDS = $cgi->param('CDS');
my $filter_tRNA_rRNA   = $cgi->param('filter-tRNA-rRNA');
my $filter_multimapped = $cgi->param('filter_multimapped');
my $mfei       = $cgi->param('mfei');
my $randfold   = $cgi->param('randfold');
my $align      = $cgi->param('align');
my $varna = 1;

if ( $filter_tRNA_rRNA   ) { $filter_tRNA_rRNA   = 1 } else { $filter_tRNA_rRNA   = 0 }
if ( $filter_multimapped ) { $filter_multimapped = 1 } else { $filter_multimapped = 0 }
if ( $filter_CDS ) { $filter_CDS = 1 } else { $filter_CDS = 0 }
if ( $mfei       ) { $mfei       = 1 } else { $mfei       = 0 }
if ( $randfold   ) { $randfold   = 1 } else { $randfold   = 0 }
if ( $align      ) { $align      = 1 } else { $align      = 0 }
if ( !$job_title ) { $job_title  = 0 }


##### Download BED file, check it and write it in the results directory
my $bedFile = '';
$bedFile   = $cgi->upload('bedFile') or miRkwood::WebTemplate::web_die("Error when getting BED file: $!");
my $localBED = $absolute_job_dir . '/' . $cgi->param('bedFile');
open (my $BED, '>', $localBED) or miRkwood::WebTemplate::web_die("Error when creating BED file: $!");
while ( <$bedFile> ){
    print $BED $_;
    if ( ! miRkwood::Utils::is_correct_BED_line($_) ){
        print $cgi->redirect($error_url . '?type=noBED');
        exit;        
    }
}
close $BED;


##### Get species or reference sequence
my $seqArea = $cgi->param('seqArea');
my $genome = '';
my $max_length = 100000;
if ( $seqArea eq q{} )    # case model organism
{
    debug('Reference species is a model organism.', 1);
    $genome = File::Spec->catfile( miRkwood::Paths->get_data_path(), 'genomes/', $species . '.fa');
}
else{
    debug('Reference sequence is provided by the user.', 1);
    
    # Check if genome is a valid fasta   
    $seqArea = miRkwood::Utils::cleanup_fasta_sequence($seqArea);
    
    if ( ! miRkwood::Utils::is_fasta($seqArea) )
    {
        print $cgi->redirect($error_url . '?type=noFasta');
        exit;
    }    
    if ( ! miRkwood::Utils::check_nb_sequences($seqArea) )
    {
        print $cgi->redirect($error_url . '?type=severalSequences');
        exit;
    }    
    if ( ! miRkwood::Utils::check_sequence_length($seqArea, $max_length) )
    {
        print $cgi->redirect($error_url . '?type=tooLongSequence');
        exit;
    }

    $genome = $absolute_job_dir . '/genome.fa';
    open (my $GENOME, '>', $genome) or miRkwood::WebTemplate::web_die("Error when creating genome file: $!");
    print $GENOME $seqArea;
    close $GENOME;    
       
}


##### Redirect to the wait.pl page until the job is done
my $arguments = '?jobId=' . $jobId . '&nameJob=' . $job_title . '&mail=' . $mail . '&mode=BAM';
my $waiting_url = miRkwood::WebTemplate::get_cgi_url('wait.pl') . $arguments;

print $cgi->redirect( $waiting_url );


##### Create config file
my $run_options_file = miRkwood::Paths->get_job_config_path($absolute_job_dir);
miRkwood->CONFIG_FILE($run_options_file);
miRkwood::write_config_for_bam_pipeline( $run_options_file, $job_title, $species, 
                                         'WebBAM', $align, $species_db, $filter_CDS, 
                                         $filter_tRNA_rRNA, $filter_multimapped, $mfei, 
                                         $randfold, $varna);


##### Launch pipeline
my $pipeline = miRkwood::BamPipelineSebastien->new($absolute_job_dir, $localBED, $genome);
$pipeline->run_pipeline();


close $log_file;
