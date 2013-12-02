#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::Paths;
use PipelineMiRNA::Utils;
use PipelineMiRNA::Results;
use PipelineMiRNA::Candidate;
use PipelineMiRNA::WebTemplate;

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();

my $bioinfo_css = PipelineMiRNA::WebTemplate->get_server_css_file();
my $project_css = File::Spec->catfile(PipelineMiRNA::Paths->get_css_path(), 'results.css');
my $js  = PipelineMiRNA::WebTemplate->get_js_file();

my $cgi            = CGI->new();
my $jobId          = $cgi->param('jobID');
my $name           = $cgi->param('nameSeq');
my $position       = $cgi->param('position');

my $candidate_name = $name.'__'.$position;
my $job = PipelineMiRNA::Results->jobId_to_jobPath($jobId);
my $returnlink = PipelineMiRNA::WebTemplate::get_link_back_to_results($jobId);
my $return_html = "<a class='returnlink' href='$returnlink'>Back to main results page</a>";

my %candidate;
my $html_contents;

if (! eval {%candidate = PipelineMiRNA::Candidate->retrieve_candidate_information($job, $name, $candidate_name);}) {
    # Catching exception
    $html_contents = "No results for the given identifiers";
}else{

    my $image_url = PipelineMiRNA::Candidate->get_relative_image(\%candidate);

    my $size = length $candidate{'DNASequence'};

    my $linkFasta = "./getCandidateFasta.pl?jobId=$jobId&name=$name&position=$position";
    my $linkVienna = "./exportVienna.pl?jobId=$jobId&name=$name&position=$position";
    my $linkAlternatives = "./exportAlternativesVienna.pl?jobId=$jobId&name=$name&position=$position";
    my $linkViennaOptimal = $linkVienna . '&optimal=1';

    my $Vienna_HTML = "<ul><li><b>Stem-loop structure (dot-bracket format):</b> <a href='$linkVienna'>download</a>";
    if($candidate{'Vienna'} ne $candidate{'Vienna_optimal'}){
        $Vienna_HTML .= "</li><li><b>Optimal MFE secondary structure (dot-bracket format):</b> <a href='$linkViennaOptimal'>download</a></li></ul>"
    } else {
        $Vienna_HTML .= "<br/><i>(This stem-loop structure is the MFE structure)</i></li></ul>"
    }

    my $alternatives_HTML = '<b>Alternative candidates:</b> ';
    if($candidate{'alternatives'}){
        $alternatives_HTML .= "<a href='$linkAlternatives'>download</a>"
    } else {
        $alternatives_HTML .= "<i>None</i>"
    }

    my $alignmentHTML;
    if($candidate{'alignment'}){
         $alignmentHTML = PipelineMiRNA::Candidate->make_alignments_HTML(\%candidate, $job, $name, $candidate_name);
    } else {
        $alignmentHTML = "No alignment has been found.";
    }

    $html_contents ="
            <div id = 'showInfo'>
        <ul>
        <li>
          <b>Name: </b>$candidate{'name'}
        </li>
        <li>
          <b>Position:</b> $position ($size nt)
        </li>
        <li>
          <b>Strand:</b>
        </li>
        <li>
          <b>Sequence (FASTA format):</b> <a href='$linkFasta'>download</a>
        </li>
        <li>
          $alternatives_HTML
        </li>
        </ul>
        <h2>Secondary structure</h2>
        <img id='structure' src='$image_url' height='400px' alt='$candidate_name secondary structure'>
        $Vienna_HTML
        <h2>Thermodynamics stability</h2>
        <ul>
        <li>
          <b>MFE:</b> $candidate{'mfe'} kcal/mol
        </li>
        <li>
          <b>AMFE:</b> $candidate{'amfe'}
        </li>
        <li>
          <b>MFEI:</b> $candidate{'mfei'}
        </li>
        </ul>
        <h2>Mirbase alignments</h2>
        $alignmentHTML

    </div><!-- showInfo -->
"
}

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
	<head>
        <link title="test" type="text/css" rel="stylesheet" href="$project_css" />
		<script src="$js" type="text/javascript" LANGUAGE="JavaScript"></script>
		<title>MicroRNA identification</title>
	</head>
	<body>
	    $return_html
        $html_contents
        $return_html
	</body>
</html>

DATA
###End###
