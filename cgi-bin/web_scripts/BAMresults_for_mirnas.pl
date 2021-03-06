#!/usr/bin/perl
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Which;
use File::Spec;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood;
use miRkwood::Paths;
use miRkwood::Utils;
use miRkwood::Pipeline;
use miRkwood::WebTemplate;
use miRkwood::Results;
use miRkwood::ResultsExporterMaker;


##### Page settings
my $header_menu = miRkwood::WebTemplate::get_header_menu( 'smallrnaseq' );
my $footer      = miRkwood::WebTemplate::get_footer();
my $web_scripts = miRkwood::WebPaths->get_web_scripts();

my @css = (
	miRkwood::WebTemplate->get_server_css_file(),
	miRkwood::WebTemplate->get_css_file(),
    miRkwood::WebTemplate->get_mirkwood_css_file()
);
my @js = ( File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'results.js' ),
	       File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'graphics.js' ),
	       File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'miARN.js' ),
	       File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'jquery.min.js' )
);

my $help_page = File::Spec->catfile( File::Spec->catdir( miRkwood::WebPaths->get_html_path(), 'smallRNAseq'), 'help.php');

##### Parameters
my $cgi = CGI->new();
my $id_job = $cgi->param('jobID');    # get id job
my $mirnas_type = $cgi->param('type');      # 'known_miRNA' or 'novel_miRNA'
my $returnlink = miRkwood::WebTemplate::get_link_back_to_BAM_results($id_job);
my $return_html = "<p><a class='returnlink' href='$returnlink'>Back to main results page</a></p>";


##### Create page
my $valid = miRkwood::Results->is_valid_jobID($id_job);
my $html = '';
my $HTML_additional = '';
my $main_contents;
my $page = '';
my $arborescence = "<a href='/'>home</a> &gt; ";
$arborescence .= "<a href='/topic/rna.html'>RNA</a> &gt; ";
$arborescence .= "<a href='/mirkwood/mirkwood.php'>mirkwood</a>";

my $mirnas_results = '';
my $mirna = '';
if ( $mirnas_type eq 'novel_miRNA' ){
    $mirna = 'Novel';
}
elsif ( $mirnas_type eq 'known_miRNA' ){
    $mirna = 'Known';
}


$HTML_additional .= '<p class="header-results" id="job_id"><b>Job ID:</b> ' . $id_job . '</p>';

if ( $valid ){

    my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($id_job);

    my $run_options_file = miRkwood::Paths->get_job_config_path($absolute_job_dir);
    miRkwood->CONFIG_FILE($run_options_file);
    my $cfg = miRkwood->CONFIG();

    my $nb_results = 0;

    my $pdf_export = '';
    #~ if ( which( 'pandoc' ) ){
        #~ $pdf_export = "<input type='radio' name='export' id='export-pdf' value='pdf' />&#160;<label for='export-pdf'>full report in PDF format</label><br/>";
    #~ }

    if ( $cfg->param('job.title') ) {
        $HTML_additional .= "<p class='header-results' id='job_title'><b>Job title:</b> " . $cfg->param('job.title') . '</p>';
    }

	unless ( miRkwood::Results->is_job_finished($id_job) ) {
		$HTML_additional .= "<p class='warning'>Still processing...</p>";
	} else {
        $nb_results = miRkwood::Results->number_of_results_bis( $id_job, $mirnas_type );
        $nb_results = miRkwood::Utils::make_numbers_more_readable( $nb_results );

        $mirnas_results .= miRkwood::Results->get_basic_pseudoXML_for_jobID($id_job, 'smallRNAseq', $mirnas_type);
    }

    if ( $nb_results != 0 ) {

        $main_contents = <<"END_TXT";
            $return_html
            <div id='new_mirnas'>
            <p class='header-results' id='precursors_count' style='font-size: 150%;'>
                <b>$mirna miRNAs : $nb_results miRNA precursor(s) found</b>
            </p> 

            <div id="select" >
                <div style="width: 510px"  class="forms">
                    <p class='text-results'>Export selected entries \(<a id='select-all' onclick='selectAll()' >select all<\/a>/<a id='deselect-all' onclick='deSelectAll()'  >deselect all</a>\) in one of the following formats:</p>
                    <form id= 'exportForm'>
                        <input type="radio" name="export" id="export-csv" checked='checked' value="csv" />&#160;<label for='export-csv'>tabular format (CSV)</label><br/>
                        <input type="radio" name="export" id="export-fas" value="fas" />&#160;<label for='export-fas'>FASTA format</label><br/>
                        <input type="radio" name="export" id="export-dot" value="dot" />&#160;<label for='export-dot'>dot-bracket format (plain sequence + secondary structure)</label><br/>
                        <input type="radio" name="export" id="export-org" value="org" />&#160;<label for='export-org'>full report in ORG-mode format</label><br/>
                        $pdf_export
                        <input type="radio" name="export" id="export-gff" value="gff" />&#160;<label for='export-gff'>GFF format</label><br/>
                        <input type="radio" name="export" id="export-reads" value="reads" />&#160;<label for='export-reads'>read cloud format</label>
                        <input style="margin-left:360px" class="myButton" type="button" name="export-button" id='export-button' value="Export" onclick='exportTo("$id_job", "$web_scripts", "smallRNAseq", "$mirnas_type")'/>
                    </form>
                </div>

                <p class='helper-results'>Click on a line to see the HTML report of a pre-miRNA prediction. Click on a checkbox to select an entry.<br/>

                <a id="hrefposition" class='on' style="color:blue;" onclick='sortBy("position","smallRNAseq")' >Sort by position <\/a> /  <a id="hrefquality" class='off' style="color:black;" onclick='sortBy("quality","smallRNAseq")'  >sort by quality</a> <a href="$help_page#quality">[?]</a>
                </p>
            </div> 

            <div id="table" > </div>
            <div id="singleCell"> </div>
            $mirnas_results

            <div id="id_job" >$id_job</div>
        </div>

        <br />
END_TXT

    } # end if nb results != 0
    else {  # 0 results
        $main_contents .= "<p class='header-results' id='precursors_count'><b>0 miRNA precursor(s) found</b></p>";
    }

}   # end if valid
else{   # job id is not a valid ID
    $main_contents .= "<p>No results available for the given job identifier $id_job.</p>";
}

$page = <<"END_TXT";
<body onload="main('all');">
    <div class="frametitle">
        <h1 id="title">miRkwood small RNA-seq</h1>
    </div>

    <div class="tabs" id="menu_central" style="display:inline-block"> 
        $header_menu
    </div>
    <div id="arborescence">$arborescence</div>

    <div id="main">
        $HTML_additional
        $main_contents
    </div>
    $footer
</body>
END_TXT

my $title = 'miRkwood - '. lc($mirna).' candidates';
$html = miRkwood::WebTemplate::get_HTML_page_for_body($page, \@css, \@js, $title);


print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
