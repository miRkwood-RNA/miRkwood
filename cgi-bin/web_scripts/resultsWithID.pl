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

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();

my $css = PipelineMiRNA::WebTemplate->get_css_file();
my $js1 = File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'results.js');
my $js2 = File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'graphics.js');
my $js3 = File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'miARN.js');
my $js4 = File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'imgpreview.full.jquery.js');

my $id_job = $cgi->param('run_id'); # récupération id job
my $name_job = $cgi->param('nameJob'); # récupération id job

my $HTML_header = <<"END_TXT";
Content-type: application/xhtml+xml

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

    <head>
        <link title="test" type="text/css" rel="stylesheet" href="/Style/bioinfo.css" />
        <link rel="stylesheet" type="text/css" href="$css" />
        <script type="text/javascript" language="Javascript" src="$js1"> </script>
        <script type="text/javascript" src="$js2"></script>
        <script type="text/javascript" src="$js3"></script>
        <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.3.1/jquery.min.js" type="text/javascript"></script>
        <script src="$js4" type="text/javascript"></script>
    </head>
END_TXT

my $HTML_additional = "";
if ($name_job ne "")
{
    $HTML_additional .= "<div class='titleJob' ><li>Title Job : ".$name_job."</li></div>";
}

my $valid = PipelineMiRNA::Results->is_valid_jobID($id_job);

if($valid){
    my %myResults = PipelineMiRNA::Results->get_structure_for_jobID($id_job);
    my $HTML_results = PipelineMiRNA::Results->resultstruct2pseudoXML( \%myResults);

    print <<"HTML";
$HTML_header    <body onload="main('all');">
    <div class="theme-border"></div>
    <div class="logo"></div>
    $bioinfo_menu
<div class="bloc_droit">

$header_menu

<div class="main main-full">
$HTML_additional
<div  id="select" > 
	<div style="width: 500px"  class="forms">
		<p style='font-size:17px;font-family: "Times New Roman", Serif' >Export selected entries \( <a onclick='selectAll()' >Select all<\/a> /  <a  onclick='deSelectAll()'  >Deselect all</a> \) :</p> 
		<form id= 'exportForm'>
		<input type="radio" name="export" checked='checked' value="csv"  />tab-delimited format (csv)<br/>
		<input type="radio" name="export" value="fas"/>fasta format (plain sequence)<br/>
		<input type="radio" name="export" value="dot"/>dot-bracket format (plain sequence + secondary structure)<br/>
		<input type="radio" name="export" value="odf"/>full report in document format (odf)<br/>
		<input type="radio" name="export" value="gff"/>gff format<br/><br/>
		<input style="margin-left:360px" class="myButton" type="button" name="bout" value="Export" onclick='exportTo("$id_job")'/>
		</form>
	</div>
		<p style='font-size:16px;font-family: "Times New Roman", Serif' ><input class="myButton" type="button" id="sort" value="Sort by quality" onclick='changeValue();'/><br/>	<br/>Click on a name to see the full HTML report. Click on the checkbox to select an entry.
		
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
</html>
HTML

}else{
    print <<"HTML";
$HTML_header
<body>
    <div class="theme-border"></div>
    <div class="logo"></div>
    $bioinfo_menu
    <div class="bloc_droit">
        $header_menu
        <div class="main">
        $HTML_additional
            <p>No results available for the given job identifier $id_job: $valid </p>
        </div><!-- main -->
    </div><!-- bloc droit-->
    $footer
    </body>
</html>
HTML
}