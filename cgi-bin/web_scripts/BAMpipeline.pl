#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Spec;
use Log::Message::Simple qw[msg error debug];

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::WebTemplate;
use miRkwood::Results;
use miRkwood::Utils;

#~ my @css = (miRkwood::WebTemplate->get_server_css_file(), miRkwood::WebTemplate->get_css_file());
#~ my @js  = (miRkwood::WebTemplate->get_js_file());


my $error_url = miRkwood::WebTemplate::get_cgi_url('error.pl');

my $cgi = new CGI;
my $job_title = $cgi->param('job');
my $mail      = $cgi->param('mail');
my $species   = $cgi->param('species');

my $jobId = miRkwood::Results->make_job_id();
my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($jobId);
mkdir $absolute_job_dir;

my $log_file = File::Spec->catfile( $absolute_job_dir, 'log.log' );
local $Log::Message::Simple::DEBUG_FH = miRkwood->LOGFH($log_file);

my $seqArea = $cgi->param('seqArea');
my $seq = q{};
if ( $seqArea eq q{} )    # case model organism
{
    debug("Reference species is a model organism", 1);

}else{
    debug("Reference sequence is provided by the user", 1);
    $seq = $seqArea;
}

$seq = miRkwood::Utils::cleanup_fasta_sequence($seq);

if ( ! miRkwood::Utils::is_fasta($seq) )
{
    print $cgi->redirect($error_url);
    exit;
}

##### Redirect to the wait.pl page until the job is done
my $arguments = '?jobId=' . $jobId . '&nameJob=' . $job_title . '&mail=' . $mail . '&mode="BAM"';
my $waiting_url = miRkwood::WebTemplate::get_cgi_url('wait.pl') . $arguments;

print $cgi->redirect( -uri => $waiting_url  );
print "Location: $waiting_url \n\n";

sleep 10;

my $is_finished_file = File::Spec->catfile( $absolute_job_dir, 'finished' );
open( my $finish, '>', $is_finished_file )
    or die "Error when opening $is_finished_file: $!";
close $finish;


#~ my $page = <<"END_TXT";
#~ <div class="main">
    #~ <p>Job title : $job_title</p>
    #~ <p>mail : $mail</p>
    #~ <p>species : $species</p>
#~ 
#~ </div><!-- main -->
#~ END_TXT
#~ 
#~ my $html = miRkwood::WebTemplate::get_HTML_page_for_content($page, \@css, \@js);
#~ print <<"DATA" or die("Error when displaying HTML: $!");
#~ Content-type: text/html
#~ 
#~ $html
#~ DATA

close $log_file;
