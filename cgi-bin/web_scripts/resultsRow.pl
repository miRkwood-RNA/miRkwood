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

my $cgi            = CGI->new();
my $jobId          = $cgi->param('jobID');
my $candidate_id   = $cgi->param('id');

my @css = (File::Spec->catfile(PipelineMiRNA::Paths->get_css_path(), 'results.css'));
my @js  = (PipelineMiRNA::WebTemplate->get_js_file());

my $job = PipelineMiRNA::Results->jobId_to_jobPath($jobId);
my $returnlink = PipelineMiRNA::WebTemplate::get_link_back_to_results($jobId);
my $return_html = "<a class='returnlink' href='$returnlink'>Back to main results page</a>";

my %candidate;
my $html_contents;

if (! eval {%candidate = PipelineMiRNA::Candidate->retrieve_candidate_information($job, $candidate_id);}) {
    # Catching exception
    $html_contents = "No results for the given identifiers";
}else{

    my $image_url = PipelineMiRNA::Candidate->get_relative_image(\%candidate);

    my $size = length $candidate{'DNASequence'};

    my $linkFasta = "./getCandidateFasta.pl?jobId=$jobId&id=$candidate_id";
    my $linkVienna = "./exportVienna.pl?jobId=$jobId&id=$candidate_id";
    my $linkAlternatives = "./exportAlternativesVienna.pl?jobId=$jobId&id=$candidate_id";
    my $linkViennaOptimal = $linkVienna . '&optimal=1';

    my $Vienna_HTML = "<ul><li><b>Stem-loop structure (dot-bracket format):</b> <a href='$linkVienna'>download</a>";
    if($candidate{'Vienna'} ne $candidate{'Vienna_optimal'}){
        $Vienna_HTML .= "</li><li><b>Optimal MFE secondary structure (dot-bracket format):</b> <a href='$linkViennaOptimal'>download</a></li></ul>"
    } else {
        $Vienna_HTML .= "<br/>This stem-loop structure is the MFE structure.</li></ul>"
    }

    my $alternatives_HTML = '<b>Alternative candidates (dot-bracket format):</b> ';
    if($candidate{'alternatives'}){
        $alternatives_HTML .= "<a href='$linkAlternatives'>download</a>"
    } else {
        $alternatives_HTML .= "<i>None</i>"
    }

    my $alignmentHTML;
    if($candidate{'alignment'}){
         $alignmentHTML = PipelineMiRNA::Candidate->make_alignments_HTML(\%candidate);
    } else {
        $alignmentHTML = "No alignment has been found.";
    }

    $html_contents = <<"END_TXT";
            <div id = 'showInfo'>
        <ul>
        <li>
          <b>Name: </b>$candidate{'name'}
        </li>
        <li>
          <b>Position:</b> $candidate{'position'} ($size nt)
        </li>
        <li>
          <b>Strand:</b> $candidate{'strand'}
        </li>
        <li>
          <b>Sequence (FASTA format):</b> <a href='$linkFasta'>download</a>
        </li>
        <li>
          $alternatives_HTML
        </li>
        </ul>
        <h2>Secondary structure</h2>
        <img id='structure' src='$image_url' height='400px' alt='$candidate{'name'} secondary structure'>
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
        <h2>miRBase alignments</h2>
        $alignmentHTML

    </div><!-- showInfo -->
END_TXT
}

my $body  = <<"END_TXT";
    <body>
       <h1>Results for $candidate{'name'}, $candidate{'position'}</h1>
        $return_html
        $html_contents
        $return_html
    </body>
END_TXT

my $html = PipelineMiRNA::WebTemplate::get_HTML_page_for_body($body, \@css, \@js);

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
###End###
