#!/usr/bin/perl -w
use Class::Struct;
use CGI;
use File::Spec;
use FindBin qw($Bin);
use Cwd qw( abs_path );
use File::Basename qw(dirname);
use CGI::Carp qw(fatalsToBrowser);
use Log::Message::Simple qw[msg error debug];

use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory

use PipelineMiRNA;
use PipelineMiRNA::Paths;
use PipelineMiRNA::Programs;
use PipelineMiRNA::WebFunctions;
use PipelineMiRNA::WebTemplate;

my $cgi = new CGI;
my $mail = $cgi->param('mail');
my $nameJob = $cgi->param('job');

my $local_dir = dirname( abs_path($0) );
my $rootdir = PipelineMiRNA::Paths->get_server_root_dir();

$dirScript = PipelineMiRNA::Paths->get_absolute_path('scripts' );    # chemin script
$dirLib    = PipelineMiRNA::Paths->get_absolute_path( 'lib' );

my $formatFasta_bin = File::Spec->catfile( $dirScript, 'formatFasta.sh' );

my $jobId = PipelineMiRNA::WebFunctions->make_job_id();
my $dirJob_name = PipelineMiRNA::WebFunctions->jobId_to_jobPath($jobId);
my $dirJob_path = PipelineMiRNA::Paths->get_absolute_path($dirJob_name);
my $root = PipelineMiRNA::Paths->get_absolute_path(PipelineMiRNA::Paths->get_results_dir_name());

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

mkdir $dirJob_path;

my $log_file = File::Spec->catfile( $dirJob_path, 'log.log' );
local $Log::Message::Simple::DEBUG_FH = PipelineMiRNA->LOGFH($log_file);

my $sequence_origin = File::Spec->catfile( $dirJob_path, 'sequence.fas' );
my $sequence_load   = File::Spec->catfile( $dirJob_path, 'sequenceLoad.fas' );
my $sequence_upload = File::Spec->catfile( $dirJob_path, 'sequenceUpload.fas' );

my $seqArea = $cgi->param('seqArea');

if ( $seqArea eq "" )    # cas upload fichier
{
    debug("Sequences are a FASTA file uploaded", 1);
    $seq = "";
    my $upload = $cgi->upload('seqFile')
      || die "Error when getting seqFile: $!";
    while ( my $ligne = <$upload> ) {
        $seq .= $ligne;
    }
    open( INPUT, '>>', $sequence_origin )
      || die "Error when opening -$sequence_origin-: $!";
    print INPUT lc($seq) . "\n";    ##mise en minuscules
    close INPUT;
    chmod 777, $sequence_origin;
    $seq = lc($seq) . "\n";

    my $fasta_cmd = "sh $formatFasta_bin $sequence_origin $sequence_load";
    debug("Run shell script $fasta_cmd", 1);
    system($fasta_cmd);             # script qui elimine les retour a la ligne

    if ( $seq !~ /^( *>.+[\r\n]+([-\. atcgunwkmsydr0-9]+[\r\n]+)+){1,}$/ )
    {                               # erreur de syntaxe
        print $cgi->redirect(
            'http://' . $ENV{SERVER_NAME} . '/cgi-bin/error.pl' );
        exit;
    }
}
else                                #cas textArea
{
    debug('Sequences are provided through the text area', 1);
    open( INPUT, '>>', $sequence_origin )
      or die "Error when opening -$sequence_origin-: $!";
    print INPUT lc($seqArea);
    system( "sed -i 'N;s/\r//g' " . $sequence_origin )
      ;    #commande sed permettant d'eliminer les caracteres du (copier-coller)
    system("sh $formatFasta_bin $sequence_origin $sequence_load")
      ;    # script qui elimine les retour a la ligne
    while ( my $ligne = <INPUT> ) {
        $seqArea .= $ligne;
    }
    $seqArea = lc($seqArea) . "\n";    #mise en miniscule

    if ( $seqArea !~ /^( *>.+[\r\n]+([-\. atcgunwkmsydr0-9]+[\r\n]+)+){1,}$/ )
    {                                  # erreur de syntaxe
        print $cgi->redirect(
            'http://' . $ENV{SERVER_NAME} . '/cgi-bin/error.pl' );
        exit;
    }
    close INPUT;
}

open( LOAD, '<', $sequence_load )
  or die "Error when opening -$sequence_load-: $!";
open( OUTPUT, '>>', $sequence_upload )
  or die "Error when opening -$sequence_upload-: $!";
while ( my $line = <LOAD> ) {
    if ( $line =~ /^>/ ) {
        $name = $line;
    }
    else {
        $line =~ s/w/n/g;
        $line =~ s/r/n/g;
        $line =~ s/k/n/g;
        $line =~ s/m/n/g;
        $line =~ s/s/n/g;
        $line =~ s/y/n/g;
        $line =~ s/d/n/g;
        print OUTPUT $name . $line;
        $name = "";
    }
}
unlink($sequence_load);

# redirection vers la page wait en attendant le calcul
my $waiting_url ='http://' . $ENV{SERVER_NAME} . "/cgi-bin/wait.pl?jobId=$jobId&nameJob=$nameJob&mail=$mail";
print $cgi->redirect( -uri => $waiting_url  );
print "Location: $waiting_url \n\n";

my $check    = $cgi->param('check');
my $mfei     = $cgi->param('mfei');
my $randfold = $cgi->param('randfold');

my $align    = $cgi->param('align');
my $plant    = $cgi->param('db');
if ( !$plant ) {
    $plant = "NonePlant";
}
if ( $mfei     eq "" ) { $mfei     = 0 }
if ( $randfold eq "" ) { $randfold = 0 }
if ( $align    eq "" ) { $align    = 0 }
if ( $check    eq "" ) { $check    = 0 }

#execution de tous les scripts de traitements
my $perl_script = File::Spec->catfile( $dirScript, 'execute_scripts.pl' );
my $cmd =
"perl -I$dirLib $perl_script $check $mfei $randfold $align $dirJob_path $plant";
debug("Running perl script $cmd", 1);
system($cmd);
my $finish_file = File::Spec->catfile( $dirJob_path, 'finished' );
open( my $finish, '>', $finish_file ) || die "$!";
close $finish_file;
debug("Writing finish file $finish_file", 1);

close $log_file;
