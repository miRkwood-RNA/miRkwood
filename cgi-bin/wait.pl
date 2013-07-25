#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use Time::gmtime;
use File::stat;
use Socket;
use File::Spec;
use Cwd qw( abs_path );
use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory
use File::Basename qw(dirname);
use CGI::Carp qw(fatalsToBrowser);

use POSIX ':sys_wait_h';
use IO::Handle;

use PipelineMiRNA::WebTemplate;

my $local_dir = dirname( abs_path($0) );
my $rootdir = abs_path( File::Spec->catdir( $local_dir, '..' ) );

my $dirScript  = File::Spec->catdir( $rootdir, 'scripts' );    # chemin script
my $dirData    = File::Spec->catdir( $rootdir, 'data' );
my $dirResults = File::Spec->catdir( $rootdir, 'results' );
my $dirLib     = File::Spec->catdir( $rootdir, 'lib' );

my $html    = CGI->new();
my $jobId   = $html->param('jobId');
my $email   = $html->param('mail');
my $nameJob = $html->param('nameJob');
my $name    = $nameJob;

my $check;
if ( $nameJob eq q{} ) { $check = 'noTitle' }

my $results_link = 'resultsWithID.pl?run_id=' . $jobId . '&nameJob=' . $name;
my $results_full_link =
  'http://' . $ENV{SERVER_NAME} . '/cgi-bin/' . $results_link;
my $dirJob_name = 'job' . $jobId;
my $dirJob      = File::Spec->catdir( $rootdir, 'results', $dirJob_name );
my $is_finished = File::Spec->catfile( $dirJob, 'finished' );

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();

#le calcul est fini
if ( -e $is_finished ) {
    my $email_script = File::Spec->catfile( $dirScript, 'email.pl' );
    my $email_cmd = "perl $email_script $nameJob $jobId";
    system($email_cmd);
    print $html->redirect( -uri => $results_full_link )
      or die("Error when redirecting: $!");
    exit;
}

my $url =
    'http://'
  . $ENV{SERVER_NAME}
  . "/cgi-bin/wait.pl?jobId=$jobId&nameJob=$name&mail=$email";

my $email_HTML;
if ( $email ne q{} ) {
    $email_HTML =
"<p>An E-mail notification will be sent to <strong>$email</strong> as soon as the job is completed.</p>";
}

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
        <meta http-equiv='Refresh' content='6;URL=$url'>
        <meta name="keywords" content="RNA, ARN, mfold, fold, structure, prediction, secondary structure" />
        <link title="test" type="text/css" rel="stylesheet" href="/arn/style/script.css" />
		<script src="/arn/js/miARN.js" type="text/javascript" LANGUAGE="JavaScript"></script>
		<title>miREST :: identification of miRNA/miRNA hairpins in plants</title>
    </head>
    <body>
        $bioinfo_menu
        <div class="bloc_droit">
            $header_menu
            <div class="main">
                <div class="dialog">
                    <br/><br/>
                    <div class="waitMessage">
                        <p>Your request has been successfully submitted.<p>
                        <p>Your ID is <B>$jobId</B>.</p>
                        $email_HTML

                        <p>This page is updated every five seconds.</p>
                        <p>You will be redirect to the <a href="./$results_link">results page</a> once the job is completed.</p>
                    </div>
                </div><!-- dialog -->
            </div><!-- main -->
            $footer
        </div><!-- bloc droit-->
    </body>
</html>
DATA
###End###
