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
use miRkwood;
use miRkwood::WebTemplate;
use miRkwood::Results;

my $html    = CGI->new();
my $jobId   = $html->param('jobId');
my $mail   = $html->param('mail');
my $nameJob = $html->param('nameJob');
$nameJob =~ s/ /_/g;#replace spaces with '_' in name job
my $name    = $nameJob;

my $res_arguments = '?run_id=' . $jobId;
my $results_page  = 'resultsWithID.pl';
my $results_link  = $results_page . $res_arguments;
my $results_baseurl = miRkwood::WebTemplate::get_cgi_url($results_page);
my $results_url   = $results_baseurl . $res_arguments;

my $wait_arguments = '?jobId=' . $jobId . '&nameJob=' . $name . '&mail=' . $mail;
my $waiting_url = miRkwood::WebTemplate::get_cgi_url('wait.pl') . $wait_arguments;

if ( miRkwood::Results->is_job_finished($jobId) ) {
    if ( $mail ne q{} ) {
        miRkwood::WebTemplate::send_email($mail, $jobId, $name);
    }
    print $html->redirect( -uri => $results_url )
      or die("Error when redirecting: $!");
    exit;
}

my $email_HTML = '';
if ( $mail ne q{} ) {
    $email_HTML =
"<p>An E-mail notification will be sent to <strong>$mail</strong><br/>as soon as the job is completed.</p>";
}

my @css = (miRkwood::WebTemplate->get_server_css_file(), miRkwood::WebTemplate->get_css_file());
my @js  = (miRkwood::WebTemplate->get_js_file());

my $page = <<"END_TXT";
<div class="main">
  <div class="dialog">
    <div class="waitMessage">
      <p>Your request has been successfully submitted.<p>
      <p>Your ID is <B>$jobId</B>.</p>
      $email_HTML

      <p>You will be redirect to the <a href="./$results_link">results page</a> once the job is completed.</p>
      <p>Please wait...</p>
    </div><!-- waitMessage -->
	<div id="waitGif"></div>
  </div><!-- dialog -->
</div><!-- main -->
END_TXT

my $html_text = miRkwood::WebTemplate::get_HTML_page_for_content($page, \@css, \@js);

$html_text =~ s/<meta/<meta http-equiv='Refresh' content='10;URL=$waiting_url'><meta/;

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html_text
DATA
###End###
