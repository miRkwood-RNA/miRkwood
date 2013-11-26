#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Cwd qw( abs_path );
use File::Basename qw(dirname);
use File::Spec;
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA;
use PipelineMiRNA::WebTemplate;

my $dirScript = PipelineMiRNA::Paths->get_scripts_path();
my $dirLib    = PipelineMiRNA::Paths->get_lib_path();
my $web_root   = PipelineMiRNA::Paths->get_web_root();
my $html    = CGI->new();
my $jobId   = $html->param('jobId');
my $mail   = $html->param('mail');
my $nameJob = $html->param('nameJob');
my $name    = $nameJob;

my $check;
if ( $nameJob eq q{} ) { $check = 'noTitle' }

my $res_arguments = '?run_id=' . $jobId . '&nameJob=' . $name;
my $results_page  = 'resultsWithID.pl';
my $results_link  = $results_page . $res_arguments;
my $results_baseurl = PipelineMiRNA::WebTemplate::make_url($results_page);
my $results_url   = $results_baseurl . $res_arguments;

my $wait_arguments = '?jobId=' . $jobId . '&nameJob=' . $name . '&mail=' . $mail;
my $waiting_url = PipelineMiRNA::WebTemplate::make_url('wait.pl') . $wait_arguments;

my $dirJob_name = 'job' . $jobId;

my $results_dir = PipelineMiRNA::Paths->get_results_filesystem_path();
my $job_dir     = File::Spec->catdir( $results_dir, $dirJob_name );
my $is_finished = File::Spec->catfile( $job_dir, 'finished' );

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();

#le calcul est fini
if ( -e $is_finished ) {
    my $email_script = File::Spec->catfile( $dirScript, 'email.pl' );
    my $email_cmd = "perl $email_script $results_baseurl $jobId $mail $nameJob ";
    system($email_cmd);
    print $html->redirect( -uri => $results_url )
      or die("Error when redirecting: $!");
    exit;
}

my $email_HTML = '';
if ( $mail ne q{} ) {
    $email_HTML =
"<p>An E-mail notification will be sent to <strong>$mail</strong> as soon as the job is completed.</p>";
}

my $bioinfo_css = PipelineMiRNA::WebTemplate->get_server_css_file();
my $project_css = PipelineMiRNA::WebTemplate->get_css_file();
my $js  = PipelineMiRNA::WebTemplate->get_js_file();

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
        <meta http-equiv='Refresh' content='6;URL=$waiting_url'>
        <meta name="keywords" content="RNA, ARN, mfold, fold, structure, prediction, secondary structure" />
        <link title="test" type="text/css" rel="stylesheet" href="$project_css" />
        <link title="test" type="text/css" rel="stylesheet" href="$bioinfo_css" />
        <script src="$js" type="text/javascript" LANGUAGE="JavaScript"></script>
        <title>miREST :: identification of miRNA/miRNA hairpins in plants</title>
    </head>
    <body>
       <div class="logo"></div> 
		
       <div class="theme-border"></div>
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
        </div><!-- bloc droit-->
      $footer
    </body>
</html>
DATA
###End###