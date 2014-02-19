#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Spec;
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::Results;
use PipelineMiRNA::WebTemplate;

my $cgi = CGI->new;

my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();
my $web_root     = PipelineMiRNA::Paths->get_web_root();


my @css = (PipelineMiRNA::WebTemplate->get_server_css_file(), PipelineMiRNA::WebTemplate->get_css_file());
my @js  = ( File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'results.js'),
            File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'graphics.js'),
            File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'miARN.js'),
            File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'jquery.min.js'),
        	File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'imgpreview.full.jquery.js')
            
          );

my $id_job = $cgi->param('run_id'); # récupération id job
my $name_job = $cgi->param('nameJob'); # récupération id job
$name_job =~ s/_/ /g;#replace _ with space in name job
my $html = '';

my $HTML_additional = "";
$HTML_additional .= "<p><b>Job ID  : </b>".$id_job.'</p>';
unless (!$name_job)
{
    $HTML_additional .= "<p><b>Job title  : </b>".$name_job.'</p>';
}

my $valid = PipelineMiRNA::Results->is_valid_jobID($id_job);

if($valid){
    my %myResults = PipelineMiRNA::Results->get_structure_for_jobID($id_job);
    my $nb_results = PipelineMiRNA::Results->number_of_results( \%myResults);
    my $HTML_results = PipelineMiRNA::Results->resultstruct2pseudoXML( \%myResults);
	$HTML_additional .= "<p><b>".$nb_results."  miRNA precursors found</b></p>";
   	my $body ="";
   	if ($nb_results != 0) {
    $body = <<"END_TXT";
<body onload="main('all');">
    <div class="theme-border"></div>
    <div class="logo"></div>
    <div class="bloc_droit">
    $header_menu
<div class="main main-full">
    $HTML_additional
    <div  id="select" > 
    	<div style="width: 500px"  class="forms">
    		<p  ><b>Export selected entries \( <a onclick='selectAll()' >Select all<\/a> /  <a  onclick='deSelectAll()'  >Deselect all</a> \) :</b></p> 
    		<form id= 'exportForm'>
    		<input type="radio" name="export" checked='checked' value="csv"  />tab-delimited format (csv)<br/>
    		<input type="radio" name="export" value="fas"/>FASTA format<br/>
    		<input type="radio" name="export" value="dot"/>dot-bracket format (plain sequence + secondary structure)<br/>
    		<input type="radio" name="export" value="odf"/>full report in document format (odf)<br/>
    		<input type="radio" name="export" value="gff"/>GFF format<br/><br/>
    		<input style="margin-left:360px" class="myButton" type="button" name="bout" value="Export" onclick='exportTo("$id_job", "$web_root")'/>
    		</form>
    	</div>
    		<p style='font-size:14px' ><br/>	Click on the line to see the HTML report of pre-miRNA. Click on the checkbox to select an entry.<br/><br/>
    		
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
   	} else {
   		$body = <<"END_TXT";
<body onload="main('all');">
    <div class="theme-border"></div>
    <div class="logo"></div>
    <div class="bloc_droit">
    $header_menu
<div class="main main-full">
    $HTML_additional
    <br>
    <h2>No result found</h2>
   
</div><!-- main -->
</div><!-- bloc droit-->
$footer
</body>
END_TXT
   	}
    $html = PipelineMiRNA::WebTemplate::get_HTML_page_for_body($body, \@css, \@js);

}else{
    my $page = <<"END_TXT";
<div class="main">
    $HTML_additional
    <p>No results available for the given job identifier $id_job: $valid </p>
</div><!-- main -->
END_TXT

    $html = PipelineMiRNA::WebTemplate::get_HTML_page_for_content($page, \@css, \@js, 1);
}
print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
