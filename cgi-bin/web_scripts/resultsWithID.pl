#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Spec;
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::Paths;
use miRkwood::Results;
use miRkwood::WebPaths;
use miRkwood::WebTemplate;
use miRkwood::ResultsExporterMaker;

my $cgi = CGI->new;

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

my $id_job = $cgi->param('run_id');    # récupération id job
my $job_path = miRkwood::Results->jobId_to_jobPath($id_job);

my $run_options_file = miRkwood::Paths->get_job_config_path($job_path);
miRkwood->CONFIG_FILE($run_options_file);
my $cfg = miRkwood->CONFIG();

my $html = '';

my $HTML_additional = "";
$HTML_additional .=
  "<p class='header-results' id='job_id'><b>Job ID:</b> " . $id_job . '</p>';
if ( $cfg->param('job.title') ) {
	$HTML_additional .=
	  "<p class='header-results' id='job_title'><b>Job title:</b> " . $cfg->param('job.title') . '</p>';
}

my $valid = miRkwood::Results->is_valid_jobID($id_job);

if ($valid) {
	my $HTML_results = '';
	my $nb_results = 0;
	unless ( miRkwood::Results->is_job_finished($id_job) ) {
		$HTML_additional .= "<p class='warning'>Still processing...</p>";
	} else {
		$HTML_results = miRkwood::Results->get_basic_pseudoXML_for_jobID($id_job);
		$nb_results = miRkwood::Results->number_of_results_bis( $id_job );
			$HTML_additional .=
	    "<p class='header-results' id='precursors_count'><b>"
	  . $nb_results
	  . "  miRNA precursor(s) found</b></p>";
	}
	my $body = "";
	if ( $nb_results != 0 ) {
		$body = <<"END_TXT";
<body onload="main('all');">
    <div class="theme-border"></div>
    <div class="logo"></div>
    <div class="bloc_droit">
    $header_menu
<div class="main main-full">
    $HTML_additional
    <div id="select" >
    	<div style="width: 510px"  class="forms">
            <p class='text-results'>Export selected entries \(<a id='select-all' onclick='selectAll()' >select all<\/a>/<a id='deselect-all' onclick='deSelectAll()'  >deselect all</a>\) in one of the following formats:</p>
    		<form id= 'exportForm'>
                <input type="radio" name="export" id="export-csv" checked='checked' value="csv" />&#160;<label for='export-csv'>tabular format (CSV)</label><br/>
                <input type="radio" name="export" id="export-fas" value="fas" />&#160;<label for='export-fas'>FASTA format</label><br/>
                <input type="radio" name="export" id="export-dot" value="dot" />&#160;<label for='export-dot'>dot-bracket format (plain sequence + secondary structure)</label><br/>
                <input type="radio" name="export" id="export-odf" value="odf" />&#160;<label for='export-odf'>full report in document format (ODF)</label><br/>
                <input type="radio" name="export" id="export-gff" value="gff" />&#160;<label for='export-gff'>GFF format</label>
                <input style="margin-left:360px" class="myButton" type="button" name="export-button" id='export-button' value="Export" onclick='exportTo("$id_job", "$web_scripts")'/>
    		</form>
    	</div>
    		<p class='helper-results'>Click on a line to see the HTML report of a pre-miRNA prediction. Click on a checkbox to select an entry.<br/>
  			<a id="hrefposition" onclick='sortBy("quality")' >Sort by position <\/a> /  <a id="hrefquality" onclick='sortBy("position")'  >sort by quality</a>
    		</p>
    </div>
    
    <div id="table" ></div>
    <div id="singleCell"> </div>
    $HTML_results
    
    <div id="id_job" >$id_job</div>
</div><!-- main -->
</div><!-- bloc droit-->
$footer
</body>
END_TXT
	}
	else {
		$body = <<"END_TXT";
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
	}
	$html =
	  miRkwood::WebTemplate::get_HTML_page_for_body( $body, \@css, \@js );

}
else {
	my $page = <<"END_TXT";
<div class="main">
    $HTML_additional
    <p>No results available for the given job identifier $id_job.</p>
</div><!-- main -->
END_TXT

	$html =
	  miRkwood::WebTemplate::get_HTML_page_for_content( $page, \@css, \@js,
		1 );
}
print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
