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
my $page = '';

my $mirna_bed;
my $final_bed;
my $other_bed;
my $cds_bed;
my $multimapped_bed;
my $initial_bed;

$HTML_additional .= '<p class="header-results" id="job_id"><b>Job ID:</b> ' . $id_job . '</p>';

if ( $valid ){

    my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($id_job);

    my $run_options_file = miRkwood::Paths->get_job_config_path($absolute_job_dir);
    miRkwood->CONFIG_FILE($run_options_file);
    my $cfg = miRkwood->CONFIG();
    
    opendir (my $dh, $absolute_job_dir) or die "Cannot open $absolute_job_dir : $!";
    while (readdir $dh) {
        if ( /_miRNAs.bed/ ){
            $mirna_bed = File::Spec->catfile($absolute_job_dir, $_);
        }
        elsif ( /_filtered.bed/ ){
            $final_bed = File::Spec->catfile($absolute_job_dir, $_);
        }        
        elsif ( /_other.bed/ ){
            $other_bed = File::Spec->catfile($absolute_job_dir, $_);
        }
        elsif ( /_CDS.bed/ ){
            $cds_bed = File::Spec->catfile($absolute_job_dir, $_);
        }
        elsif ( /_multimapped.bed/ ){
            $multimapped_bed = File::Spec->catfile($absolute_job_dir, $_);
        }
        elsif ( /.bed/ ){
            $initial_bed = File::Spec->catfile($absolute_job_dir, $_);
        }        
    }
    closedir $dh;    

    my $nb_new_results   = 0;
    my $nb_known_results = 0;
    my $nb_total_reads   = 0;
    my $nb_CDS_reads     = 0;
    my $nb_other_reads   = 0;
    my $nb_multi_reads   = 0;

    if ( $cfg->param('job.title') ) {
        $HTML_additional .= "<p class='header-results' id='job_title'><b>Job title:</b> " . $cfg->param('job.title') . '</p>';
    }

	unless ( miRkwood::Results->is_job_finished($id_job) ) {
		$HTML_additional .= "<p class='warning'>Still processing...</p>";
	} else {
        $nb_new_results   = miRkwood::Results->number_of_results_bis( $id_job, 'new' );
        $nb_known_results = miRkwood::Results->number_of_results_bis( $id_job, 'known' );
        $nb_total_reads   = miRkwood::BEDHandler::count_reads_in_bed_file( $initial_bed );

        my $arguments = '?jobID=' . $id_job;
        my $known_url = miRkwood::WebTemplate::get_cgi_url('BAMresults_for_mirnas.pl') . $arguments . '&type=known';
        my $new_url = miRkwood::WebTemplate::get_cgi_url('BAMresults_for_mirnas.pl') . $arguments . '&type=new';

        $HTML_additional .= "<ul>";
        $HTML_additional .= "<li>Total number of reads: $nb_total_reads (XXX unique reads)</li>";

        if ( -r $cds_bed ){
            $nb_CDS_reads     = miRkwood::BEDHandler::count_reads_in_bed_file( $cds_bed );
            $HTML_additional .= "<li>CoDing Sequences: $nb_CDS_reads reads (download)</li>";            
        }
        if ( -r $nb_other_reads ){
            $nb_other_reads   = miRkwood::BEDHandler::count_reads_in_bed_file( $other_bed );
            $HTML_additional .= "<li>rRNA/tRNA: $nb_other_reads reads (download)</li>";
        }
        if ( -r $nb_multi_reads ){
            $nb_multi_reads   = miRkwood::BEDHandler::count_reads_in_bed_file( $multimapped_bed ); 
            $HTML_additional .= "<li>Frequent reads: $nb_multi_reads reads (download)</li>";
        }

        $HTML_additional .= "<li>Known miRNAs: $nb_known_results sequence(s) (<a href=$known_url>see results</a>)</li>";
        $HTML_additional .= "<li>Novel miRNAs: $nb_new_results sequence(s) (<a href=$new_url>see results</a>)</li>";
        $HTML_additional .= "</ul>";

    }

    #~ $HTML_additional .= "</div>";


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
