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
use miRkwood::WebPaths;
use miRkwood::Results;


##### Page settings
my @css = (miRkwood::WebTemplate->get_server_css_file(), miRkwood::WebTemplate->get_css_file());
my @js  = (miRkwood::WebTemplate->get_js_file());


##### Parameters
my $html    = CGI->new();
my $jobId   = '';
my $mail    = '';
my $nameJob = '';
$jobId   = $html->param('jobId');
$mail    = $html->param('mail');
$nameJob = $html->param('nameJob');
$nameJob =~ s/ /_/g; # replace spaces with '_' in name job
$jobId   =~ s/\s//g; # delete possible spaces before or after the job ID

my $page = '';
my $html_text = '';

my $valid = miRkwood::Results->is_valid_jobID($jobId);
my $retrieve_results_page = File::Spec->catfile( miRkwood::WebPaths->get_html_path(), 'id.php');

##### Create page
if ( $valid ){

    my $results_page  = 'nonImplemented.pl';
    my $res_arguments = '?run_id=' . $jobId;
    my $mode = '';

    if ( $jobId =~ /BAM/ ){
        $results_page  = 'BAMresults.pl';
        $mode = 'WebBAM';
    }
    else {
        $results_page  = 'resultsWithID.pl';
        $mode = 'fasta';
    }

    my $results_baseurl = miRkwood::WebTemplate::get_cgi_url($results_page);
    my $results_url     = $results_baseurl . $res_arguments;
    my $wait_arguments  = '?jobId=' . $jobId . '&nameJob=' . $nameJob . '&mail=' . $mail;
    my $waiting_url = miRkwood::WebTemplate::get_cgi_url('wait.pl') . $wait_arguments;

    if ( miRkwood::Results->is_job_finished($jobId) ) {
        if ( $mail ne q{} ) {
            miRkwood::WebTemplate::send_email($mode, $mail, $jobId, $nameJob);
        }
        print $html->redirect( -uri => $results_url )
          or die("Error when redirecting: $!");
        exit;
    }
    else{
        my $HTML_additional = "<p class='header-results' id='job_id'><b>Job ID:</b> $jobId</p>";
        $page = <<"END_TXT";
<div class="main">
    $HTML_additional
    <p>Sorry, something went wrong with miRkwood. Your job is not finished or it has crashed. <br />
    No results available for the given job identifier $jobId.</p>
</div><!-- main -->
END_TXT

        $html_text = miRkwood::WebTemplate::get_HTML_page_for_content( 'static/', $page, \@css, \@js, 1 );        
    }


}
else{
    my $HTML_additional = "<p class='header-results' id='job_id'><b>Job ID:</b> $jobId</p>";
	$page = <<"END_TXT";
<div class="main">
    $HTML_additional
    <p>No results available for the given job identifier $jobId.</p>

    <br />

    <input type="button" name="nom" value="Submit another ID" onclick="window.location.href='$retrieve_results_page'" />
</div><!-- main -->
END_TXT

	$html_text = miRkwood::WebTemplate::get_HTML_page_for_content( 'static/', $page, \@css, \@js, 1 );

}

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html_text

DATA

###End###
