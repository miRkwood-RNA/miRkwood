#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
my $cgi = CGI->new();
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use FindBin;                     # locate this script
use lib "$FindBin::Bin/../lib";  # use the parent directory
use PipelineMiRNA::Components;
use PipelineMiRNA::WebTemplate;
use PipelineMiRNA::WebFunctions;
use PipelineMiRNA::Utils;

my $name     = $cgi->param('nameSeq');
my $factor   = $cgi->param('factor');
my $value    = $cgi->param('value');
my $position = $cgi->param('position');
my $typePage = $cgi->param('typePage');
my $url      = $cgi->param('url');

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();


sub get_first_element_of_split {
    my @args = @_;
    my $value = shift @args;
    my @split = split(/-/, $value);
    return $split[0];
}

=method make_HTML

Returns the HTML page with the given <body>

=cut

sub make_HTML {
    my @args = @_;
    my $body = shift @args;
    my $html = <<"DATA";
Content-type: text/html

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html>
    <head>
        <LINK rel="stylesheet" type="text/css" href="/arn/style/script.css" />
        <link href="/arn/css/basic.css" type="text/css" rel="stylesheet" />
        <script type="text/javascript" language="Javascript" src="/arn/js/results.js"> </script>
        <title>MicroRNA identification</title>
    </head>
    <body>
        <div class="theme-border"></div>
        <div class="logo"></div>

        <div class="bloc_droit">

            $header_menu
            <div class="main main-full">
                $body
            </div><!-- main -->
            $footer
        </div><!-- bloc droit-->
    </body>
</html>
DATA
    return $html;
}

# Script code

my $body;
if ( $typePage eq 'simpleCell' ) {
    $body = <<"DATA";
        <div class="titreDiv"> MicroRNA identification results:</div>
        <div id = 'showInfo'>
            <h2 ><u>Sequence Informations </u></h2><br/>
            <li><b>Name sequence:</b> ${name}$position
            </li>
            <li><b>$factor:</b> $value
            </li>
        </div>
DATA
    my $html = make_HTML($body);
}
elsif ( $typePage eq 'alignement' ) {

    my %results;
    if (! eval {%results = PipelineMiRNA::Components::parse_custom_exonerate_output($url);}) {
        # Catching exception
        my $error = "No alignment available";
        print PipelineMiRNA::WebTemplate::get_error_page($error);
        die($error);
    }

    my @url = split( /\//xms, $url );
    my $length = scalar(@url);
    my $jobId  = $url[$length - 4];
    my $len    = length($jobId);
    my $jobId2 = substr($jobId, 3, $len - 3);
    my $dir    = $url[$length - 3];
    my $subDir = $url[$length - 2];

    my $job = PipelineMiRNA::WebFunctions->jobId_to_jobPath($jobId2);
    my %candidate = PipelineMiRNA::WebFunctions->retrieve_candidate_information($job, $dir, $subDir);
    my $sequence = $candidate{'DNASequence'};
    my $structure = $candidate{'Vienna'};
    my $hairpin = PipelineMiRNA::Utils::make_ASCII_viz($sequence, $structure);

    my $contents = "";

    my $predictionCounter = 0;
    my @keys = sort { get_first_element_of_split($a)  <=> get_first_element_of_split($b) } keys %results;
    foreach my $position (@keys) {
        $predictionCounter += 1;
        # Sorting the hit list by descending value of the 'score' element
        my @hits = sort { $b->{'score'} <=> $a->{'score'} } @{$results{$position}};
        $contents .= "<h3>Predition #$predictionCounter: $position</h3>
        <pre style='height: 100px;'>$hairpin</pre>
        <ul>
            <li>Evaluation score of MIRdup: TODO</li>
        </ul>
        <h4>Alignments</h4>
        ";
        foreach my $hit (@hits){
            my $alignment = $hit->{'alignment'};
            my $name = $hit->{'name'};
            my ($top, $middle, $bottom) = split(/\n/, $alignment);
            $top    = sprintf "%-15s %3s %s %s", 'query', $hit->{'begin_query'},  $top,    $hit->{'end_query'};
            $middle = sprintf "%-15s %3s %s %s", '',      '',                     $middle, '';
            $bottom = sprintf "%-15s %3s %s %s", $name,   $hit->{'begin_target'}, $bottom, $hit->{'end_target'};
            $contents .= <<"INNER";
<pre>
$top
$middle
$bottom
</pre>
INNER
        }

    }

    $body = <<"DATA";
<h2> Alignments :</h2>
$contents
DATA

}
elsif ( $typePage eq 'image' ) {
    my @url = split( /\//xms, $url );
    my $image_name = @url[5];
    $body = <<"DATA";
		<div class="titreDiv"> MicroRNA identification results:</div>
		<h2> Structure :</h2>
		<div class="figure" >
		  <img src='$url' border=0 alt='image'>
		  <p>Fig : $image_name sequence
    	</div>	
DATA
}

print make_HTML($body) or die("Error when displaying HTML: $!");

###End###
