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
use miRkwood::BEDHandler;
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
my $HTML_additional = '';
my $HTML_results = '';
my $page = '';


$HTML_additional .= '<p class="header-results" id="job_id"><b>Job ID:</b> ' . $id_job . '</p>';

if ( $valid ){

    my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($id_job);

    my $run_options_file = miRkwood::Paths->get_job_config_path($absolute_job_dir);
    miRkwood->CONFIG_FILE($run_options_file);
    my $cfg = miRkwood->CONFIG();

    my $initial_bed     = miRkwood::Paths::get_bed_file ( $id_job, '' );
    my $mirna_bed       = miRkwood::Paths::get_bed_file ( $id_job, '_miRNAs' );
    my $final_bed       = miRkwood::Paths::get_bed_file ( $id_job, '_filtered' );
    my $other_bed       = miRkwood::Paths::get_bed_file ( $id_job, '_otherRNA' );
    my $cds_bed         = miRkwood::Paths::get_bed_file ( $id_job, '_CDS' );
    my $multimapped_bed = miRkwood::Paths::get_bed_file ( $id_job, '_multimapped' );

    my $nb_new_results   = 0;
    my $nb_known_results = 0;
    my $nb_total_reads   = 0;
    my $nb_CDS_reads     = 0;
    my $nb_other_reads   = 0;
    my $nb_multi_reads   = 0;
    my $nb_total_reads_unq = 0;
    my $nb_CDS_reads_unq   = 0;
    my $nb_other_reads_unq = 0;
    my $nb_multi_reads_unq = 0;

    if ( $cfg->param('job.title') ) {
        $HTML_additional .= "<p class='header-results' id='job_title'><b>Job title:</b> " . $cfg->param('job.title') . '</p>';
    }

	unless ( miRkwood::Results->is_job_finished($id_job) ) {
		$HTML_additional .= "<p class='warning'>Still processing...</p>";
	} else {
        ##### Summary of options
        my $basename_bed = '';
        if ( $initial_bed =~ /.*(\/|\\)([^\/]+)/ ){
            $basename_bed = $2;
        }
        $HTML_additional .= "<div class='results_summary'><ul>";
        $HTML_additional .= '<h2>Options summary:</h2>';
        $HTML_additional .= '<br />';
        $HTML_additional .= "<li><b>BED file:</b> $basename_bed</li>";

        # Reference species
        if ( $cfg->param('job.plant') ){
            $HTML_additional .= '<li><b>Reference species:</b> ' . $cfg->param('job.plant') . '</li>';
        }

        # Align
        if ( $cfg->param('options.align') ){
            $HTML_additional .= '<li><b>Flag conserved mature miRNAs:</b> yes</li>';
        }
        else{
            $HTML_additional .= '<li><b>Flag conserved mature miRNAs:</b> no</li>';
        }

        # MFEI
        if ( $cfg->param('options.mfei') ){
            $HTML_additional .= '<li><b>Select only sequences with MFEI < -0.6:</b> yes</li>';
        }
        else{
            $HTML_additional .= '<li><b>Select only sequences with MFEI < -0.6:</b> no</li>';
        }

        # Ranfold
        if ( $cfg->param('options.ranfold') ){
            $HTML_additional .= '<li><b>Compute thermodynamic stability:</b> yes</li>';
        }
        else{
            $HTML_additional .= '<li><b>Compute thermodynamic stability:</b> no</li>';
        }                

        # CDS
        if ( $cfg->param('options.filter_CDS') ){
            $HTML_additional .= '<li><b>Filter CoDing Sequences:</b> yes</li>';
        }
        else{
            $HTML_additional .= '<li><b>Filter CoDing Sequences:</b> no</li>';
        }

        # tRNA and rRNA
        if ( $cfg->param('options.filter_tRNA_rRNA') ){
            $HTML_additional .= '<li><b>Filter tRNA and rRNA:</b> yes</li>';
        }
        else{
            $HTML_additional .= '<li><b>Filter tRNA and rRNA:</b> no</li>';
        }

        # Multimapped reads
        if ( $cfg->param('options.filter_multimapped') ){
            $HTML_additional .= '<li><b>Filter multimapped reads:</b> yes</li>';
        }
        else{
            $HTML_additional .= '<li><b>Filter multimapped reads:</b> no</li>';
        }

        $HTML_additional .= '</ul></div>';

        ##### Summary of results
        $nb_new_results   = miRkwood::Results->number_of_results_bis( $id_job, 'New' );
        $nb_known_results = miRkwood::Results->number_of_results_bis( $id_job, 'Known' );
        ($nb_total_reads, $nb_total_reads_unq) = miRkwood::BEDHandler::count_reads_in_bed_file( $initial_bed );

        my $arguments = '?jobID=' . $id_job;
        my $known_url = miRkwood::WebTemplate::get_cgi_url('BAMresults_for_mirnas.pl') . $arguments . '&type=Known';
        my $new_url = miRkwood::WebTemplate::get_cgi_url('BAMresults_for_mirnas.pl') . $arguments . '&type=New';

        $HTML_results .= "<div class='results_summary'><ul>";
        $HTML_results .= '<h2>Results summary:</h2>';
        $HTML_results .= '<br />';
        $HTML_results .= "<li>Total number of reads: $nb_total_reads ($nb_total_reads_unq unique reads)</li>";

        if ( $cfg->param('options.filter_CDS') ){
            ($nb_CDS_reads, $nb_CDS_reads_unq) = miRkwood::BEDHandler::count_reads_in_bed_file( $cds_bed );
            if ( $nb_CDS_reads > 0 ){
                $HTML_results .= "<li>CoDing Sequences: $nb_CDS_reads reads (download)</li>";
            }
            else {
                $HTML_results .= '<li>CoDing Sequences: 0 reads</li>';
            }
        }
        if ( $cfg->param('options.filter_tRNA_rRNA') ){
            ($nb_other_reads, $nb_other_reads_unq) = miRkwood::BEDHandler::count_reads_in_bed_file( $other_bed );
            if ( $nb_other_reads > 0 ){
                $HTML_results .= "<li>rRNA/tRNA: $nb_other_reads reads (download)</li>";
            }
            else {
                $HTML_results .= '<li>rRNA/tRNA: 0 reads</li>';
            }
        }
        if ( $cfg->param('options.filter_multimapped') ){
            ($nb_multi_reads, $nb_multi_reads_unq) = miRkwood::BEDHandler::count_reads_in_bed_file( $multimapped_bed );
            if ( $nb_multi_reads > 0 ){
                $HTML_results .= "<li>Frequent reads: $nb_multi_reads reads (download)</li>";
            }
            else {
                $HTML_results .= '<li>Frequent reads: 0 reads</li>';
            }
        }

        $HTML_results .= "<li>Known miRNAs: $nb_known_results sequence(s) (<a href=$known_url>see results</a>)</li>";
        $HTML_results .= "<li>Novel miRNAs: $nb_new_results sequence(s) (<a href=$new_url>see results</a>)</li>";
        $HTML_results .= "</ul></div>";

    }


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
            
            $HTML_results
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
