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


##### Page settings
my @css = (miRkwood::WebTemplate->get_server_css_file(), miRkwood::WebTemplate->get_css_file());
my @js  = (miRkwood::WebTemplate->get_js_file());


##### Parameters
my $html    = CGI->new();
my $jobId   = $html->param('jobId');
my $mail    = $html->param('mail');
my $nameJob = $html->param('nameJob');
$nameJob =~ s/ /_/g; # replace spaces with '_' in name job

my $page = '';
my $html_text = '';

my $valid = miRkwood::Results->is_valid_jobID($jobId);


##### Create page
if ( $valid ){

    my $results_page  = 'nonImplemented.pl';
    my $res_arguments = '?run_id=' . $jobId;

    if ( $jobId =~ /BAM/ ){
        $results_page  = 'BAMresults.pl'
    }
    else {
        $results_page  = 'resultsWithID.pl';
    }

    my $results_baseurl = miRkwood::WebTemplate::get_cgi_url($results_page);
    my $results_url     = $results_baseurl . $res_arguments;
    my $wait_arguments  = '?jobId=' . $jobId . '&nameJob=' . $nameJob . '&mail=' . $mail;
    my $waiting_url = miRkwood::WebTemplate::get_cgi_url('wait.pl') . $wait_arguments;

    if ( miRkwood::Results->is_job_finished($jobId) ) {
        if ( $mail ne q{} ) {
            miRkwood::WebTemplate::send_email($mail, $jobId, $nameJob);
        }
        print $html->redirect( -uri => $results_url )
          or die("Error when redirecting: $!");
        exit;
    }

    my $email_HTML = '';
    if ( $mail ne q{} ) {
        $email_HTML = "<p>An E-mail notification will be sent to <strong>$mail</strong><br/>as soon as the job is completed.</p>";
    }


    $page = <<"END_TXT";
<div class="main">
  <div class="dialog">
    <div class="waitMessage">
      <p>Your request has been successfully submitted.<p>
      <p>Your ID is <B>$jobId</B>.</p>
      $email_HTML

      <p>You will be redirected to the <a href="$results_url">results page</a> once the job is completed.</p>
      <p>Please wait...</p>
    </div><!-- waitMessage -->
	<div id="waitGif"></div>
  </div><!-- dialog -->
</div><!-- main -->
END_TXT

    $html_text = miRkwood::WebTemplate::get_HTML_page_for_content($page, \@css, \@js);

    $html_text =~ s/<meta/<meta http-equiv='Refresh' content='10;URL=$waiting_url'><meta/;

}
else{
    my $HTML_additional = "<p class='header-results' id='job_id'><b>Job ID:</b>$jobId</p>";
	$page = <<"END_TXT";
<div class="main">
    $HTML_additional
    <p>No results available for the given job identifier $jobId.</p>
</div><!-- main -->
END_TXT

	$html_text = miRkwood::WebTemplate::get_HTML_page_for_content( $page, \@css, \@js, 1 );

}

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html_text
DATA

###End###
