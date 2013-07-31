#!/usr/bin/perl -w
use Class::Struct;
use CGI;
use Time::gmtime;
use File::Spec;
use FindBin qw($Bin);
use Cwd qw( abs_path );
use File::Basename qw(dirname);
use CGI::Carp qw(fatalsToBrowser);
use Log::Message::Simple qw[msg error debug];

use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory

use PipelineMiRNA::Paths;

my $cgi = new CGI;
$now = gmctime();

$now =~ s/[: ]//g;
$now = substr( $now, 3 );

my $local_dir = dirname( abs_path($0) );
my $rootdir = PipelineMiRNA::Paths->get_server_root_dir();

$dirScript = PipelineMiRNA::Paths->get_absolute_path('scripts' );    # chemin script
$dirLib    = PipelineMiRNA::Paths->get_absolute_path( 'lib' );

my $formatFasta_bin = File::Spec->catfile( $dirScript, 'formatFasta.sh' );

$dirJob_name = 'job' . $now;
$dirJob = PipelineMiRNA::Paths->get_absolute_path('results', $dirJob_name );

#$dirJob = $dirData."job".$now."/"; # chemin séquence de base
mkdir $dirJob;

my $log_file = File::Spec->catfile( $dirJob, 'log.log' );
open( my $LOG, '>>', $log_file ) or die "Error when opening $log_file: $!";

my $sequence_origin = File::Spec->catfile( $dirJob, 'sequence.fas' );
my $sequence_load   = File::Spec->catfile( $dirJob, 'sequenceLoad.fas' );
my $sequence_upload = File::Spec->catfile( $dirJob, 'sequenceUpload.fas' );

my $seqArea = $cgi->param('seqArea');
if ( $seqArea eq "" )    # cas upload fichier
{
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
    print $LOG "Captured $seq";
    my $fasta_cmd = "sh $formatFasta_bin $sequence_origin $sequence_load";
    print $LOG "Run $fasta_cmd";
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
    open( INPUT, '>>', $sequence_origin )
      || die "Error when opening -$sequence_origin-: $!";
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
my $mail = $cgi->param('mail');

my $nameJob = $cgi->param('job');

print $cgi->redirect( -uri => 'http://'
      . $ENV{SERVER_NAME}
      . '/cgi-bin/wait.pl?jobId='
      . $now
      . '&mail='
      . $mail
      . '&nameJob='
      . $nameJob );

print $cgi->redirect( -uri => 'http://'
      . $ENV{SERVER_NAME}
      . '/cgi-bin/wait.pl?jobId='
      . $now
      . '&mail='
      . $mail );
print "Location: http://"
  . $ENV{SERVER_NAME}
  . "wait.pl?jobId="
  . $now
  . "&mail="
  . $mail . "\n\n";

my $check    = $cgi->param('check');
my $mfei     = $cgi->param('mfei');
my $randfold = $cgi->param('randfold');
my $SC       = $cgi->param('selfContain');
my $align    = $cgi->param('align');
my $plant    = $cgi->param('db');
if ( !$plant ) {
    $plant = "NonePlant";
}
if ( $mfei     eq "" ) { $mfei     = "notChecked" }
if ( $randfold eq "" ) { $randfold = "notChecked" }
if ( $SC       eq "" ) { $SC       = "notChecked" }
if ( $align    eq "" ) { $align    = "notChecked" }
if ( $check    eq "" ) { $check    = "notChecked" }

#execution de tous les scripts de traitements
my $perl_script = File::Spec->catfile( $dirScript, 'execute_scripts.pl' );
my $cmd =
"perl -I$dirLib $perl_script $check $mfei $randfold $SC $align $dirJob $plant";

print $LOG "The command is: $cmd";

system($cmd);

my $finish_file = File::Spec->catfile( $dirJob, 'finished' );
open( my $finish, '>', $finish_file ) || die "$!";
close $finish_file;
close $log_file;
