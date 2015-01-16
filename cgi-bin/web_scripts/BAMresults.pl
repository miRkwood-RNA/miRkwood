#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Spec;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood;
use miRkwood::Paths;
use miRkwood::Pipeline;
use miRkwood::WebTemplate;
use miRkwood::Results;
use miRkwood::ResultsExporterMaker;


##### Page settings
my $header_menu = miRkwood::WebTemplate::get_header_menu();
my $footer      = miRkwood::WebTemplate::get_footer();
my $web_scripts = miRkwood::WebPaths->get_web_scripts();

my @css = (
	miRkwood::WebTemplate->get_server_css_file(),
	miRkwood::WebTemplate->get_css_file()
);
my @js = (
	File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'results.js' ),
	File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'graphics.js' ),
	File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'miARN.js' ),
	File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'jquery.min.js' ),
	File::Spec->catfile( miRkwood::WebPaths->get_js_path(),
		'imgpreview.full.jquery.js' )

);

##### Parameters
my $cgi = CGI->new();
my $id_job = $cgi->param('run_id');    # get id job



##### Create page
my $valid = miRkwood::Results->is_valid_jobID($id_job);
my $html = '';
my $HTML_additional = "<div class='forms'>";
my $page = '';


$HTML_additional .= '<p class="header-results" id="job_id"><b>Job ID:</b> ' . $id_job . '</p>';

if ( $valid ){

    my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($id_job);

    my $run_options_file = miRkwood::Paths->get_job_config_path($absolute_job_dir);
    miRkwood->CONFIG_FILE($run_options_file);
    my $cfg = miRkwood->CONFIG();

    my $nb_results = 0;
    my $nb_known_results = 0;

    if ( $cfg->param('job.title') ) {
        $HTML_additional .= "<p class='header-results' id='job_title'><b>Job title:</b> " . $cfg->param('job.title') . '</p>';
    }

	unless ( miRkwood::Results->is_job_finished($id_job) ) {
		$HTML_additional .= "<p class='warning'>Still processing...</p>";
	} else {
        $nb_results = miRkwood::Results->number_of_results_bis( $id_job, 'new' );
        $nb_known_results = miRkwood::Results->number_of_results_bis( $id_job, 'known' );

        my $arguments = '?jobID=' . $id_job;
        my $known_url = miRkwood::WebTemplate::get_cgi_url('BAMresults_for_mirnas.pl') . $arguments . '&type=known';
        my $new_url = miRkwood::WebTemplate::get_cgi_url('BAMresults_for_mirnas.pl') . $arguments . '&type=new';

        $HTML_additional .= "<ul>";
        $HTML_additional .= "<li>Total number of reads: XXX (XXX unique reads)</li>";
        $HTML_additional .= "<li>rRNA/tRNA: XXX reads (download)</li>";
        $HTML_additional .= "<li>CoDing Sequences: XXX reads (download)</li>";
        $HTML_additional .= "<li>Frequent reads: XXX reads (download)</li>";
        $HTML_additional .= "<li>Known miRNAs: $nb_known_results sequence(s) (<a href=$known_url>see results</a>)</li>";
        $HTML_additional .= "<li>Novel miRNAs: $nb_results sequence(s) (<a href=$new_url>see results</a>)</li>";
        $HTML_additional .= "</ul>";

    }

    $HTML_additional .= "</div>";


    $page = <<"END_TXT";
<body onload="main('all');">
    <div class="theme-border"></div>
    <a href="/">
        <div class="logo"></div>
    </a>
    <div class="bloc_droit">
        $header_menu
        <div class="main main-full">
            $HTML_additional
        </div><!-- main -->
    </div><!-- bloc droit-->
    $footer
</body>
END_TXT
    
    $html = miRkwood::WebTemplate::get_HTML_page_for_body($page, \@css, \@js);

}   # end if valid
else{   # job id is not a valid ID
    $HTML_additional .= '</div>';
	$page = <<"END_TXT";
<div class="main">
    $HTML_additional
    <p>No results available for the given job identifier $id_job.</p>
</div><!-- main -->
END_TXT

	$html = miRkwood::WebTemplate::get_HTML_page_for_content( $page, \@css, \@js, 1 );
}


print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
