#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use Cwd qw( abs_path );
use File::Basename qw(dirname);
use File::Spec;
use Log::Message::Simple qw[msg error debug];
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA;
use PipelineMiRNA::Paths;
use PipelineMiRNA::Programs;
use PipelineMiRNA::Results;
use PipelineMiRNA::Utils;
use PipelineMiRNA::WebTemplate;

my $cgi = new CGI;
my $mail = $cgi->param('mail');
my $filter   = $cgi->param('CDS');
my $mfei     = $cgi->param('mfei');
my $randfold = $cgi->param('randfold');
my $align    = $cgi->param('align');
my $plant    = $cgi->param('db');
my $strand   = $cgi->param('strand');
my $job_title=  $cgi->param('job');
if ( !$plant ) {
    $plant = "NonePlant";
}

my $local_dir = dirname( abs_path($0) );

my $dirScript = PipelineMiRNA::Paths->get_scripts_path();
my $dirLib    = PipelineMiRNA::Paths->get_lib_path();

my $error_url = PipelineMiRNA::WebTemplate::get_cgi_url('error.pl');
my $root = PipelineMiRNA::Paths->get_results_filesystem_path();

if (! -e $root) {
    my $error = "Designated directory ($root) for results does not exist. Please contact the system administrator";
    print PipelineMiRNA::WebTemplate::get_error_page($error);
    die($error);
}

if (! -W $root) {
    my $error = "Cannot write results in designated directory $root. Please contact the system administrator";
    print PipelineMiRNA::WebTemplate::get_error_page($error);
    die($error);
}

my @unavailable = PipelineMiRNA::Programs::list_unavailable_programs();
if (@unavailable){
    my $error = "Cannot find required third-party software: @unavailable. Please contact the system administrator";
    print PipelineMiRNA::WebTemplate::get_error_page($error);
    die($error);
}

my $jobId = PipelineMiRNA::Results->make_job_id();
my $absolute_job_dir = PipelineMiRNA::Results->jobId_to_jobPath($jobId);
mkdir $absolute_job_dir;

my $log_file = File::Spec->catfile( $absolute_job_dir, 'log.log' );
local $Log::Message::Simple::DEBUG_FH = PipelineMiRNA->LOGFH($log_file);

my $sequence_origin_file = File::Spec->catfile( $absolute_job_dir, 'sequence.fas' );
my $sequence_upload = File::Spec->catfile( $absolute_job_dir, 'input_sequences.fas' );

my $seqArea = $cgi->param('seqArea');
my $seq = q{};
if ( $seqArea eq q{} )    # cas upload fichier
{
    debug("Sequences are a FASTA file uploaded", 1);
    my $upload = $cgi->upload('seqFile')
      || die "Error when getting seqFile: $!";
    while ( my $ligne = <$upload> ) {
        $seq .= $ligne;
    }
}else{
    debug('Sequences are provided through the text area', 1);
    $seq = $seqArea;
}

$seq = PipelineMiRNA::Utils::cleanup_fasta_sequence($seq);

if ( ! PipelineMiRNA::Utils::is_fasta($seq) )
{
    print $cgi->redirect($error_url);
    exit;
}

open( my $OUTPUT, '>>', $sequence_upload )
  or die "Error when opening -$sequence_upload-: $!";
my $name;
my @lines = split /\n/, $seq;
foreach my $line (@lines) {
    if ( $line =~ /^>/smx ) {
        $name = $line;
    }
    else {
        my $cleaned_line = PipelineMiRNA::Utils::mask_sequence_nucleotide($line);
        print $OUTPUT $name . "\n" . $cleaned_line . "\n";
        $name = "";
    }
}
close $OUTPUT or die("Error when closing $sequence_upload: $!");

# redirection vers la page wait en attendant le calcul
my $arguments = '?jobId=' . $jobId . '&nameJob=' . $job_title . '&mail=' . $mail;
my $waiting_url = PipelineMiRNA::WebTemplate::get_cgi_url('wait.pl') . $arguments;


print $cgi->redirect( -uri => $waiting_url  );
print "Location: $waiting_url \n\n";

if ( $strand   ) { $strand   = 1 } else { $strand   = 0 }
if ( $mfei     ) { $mfei     = 1 } else { $mfei     = 0 }
if ( $randfold ) { $randfold = 1 } else { $randfold = 0 }
if ( $align    ) { $align    = 1 } else { $align    = 0 }
if ( $filter   ) { $filter   = 1 } else { $filter  = 0 }
if ( !$job_title ) {
    $job_title = 0;
}
my $varna = 1;
my $run_options_file = PipelineMiRNA::Paths->get_job_config_path($absolute_job_dir);
PipelineMiRNA->CONFIG_FILE($run_options_file);
PipelineMiRNA::write_config( $run_options_file, $strand, $filter, $mfei, $randfold, $align, $job_title, $plant, $varna, 'fasta' );

# execution de tous les scripts de traitements
my $perl_script = File::Spec->catfile( $dirScript, 'execute_scripts.pl' );
my $cmd =
"perl -I$dirLib $perl_script $absolute_job_dir";
debug("Running perl script $cmd", 1);
system($cmd);
debug("Getting back from Perl script", 1);

close $log_file;