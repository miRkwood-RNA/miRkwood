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

=method make_HTML

Returns the HTML page with the given <body>

=cut

sub make_HTML {
    my @args = @_;
    my $body = shift @args;
    my $html = <<"DATA";
Content-type: text/html

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

$bioinfo_menu

<div class="bloc_droit">

$header_menu

<div class="main">
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
    my %results = PipelineMiRNA::Components::parse_custom_exonerate_output($url);
    
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


    keys %results;
    while(my($position, $values) = each %results) {
        $contents .= "<h3>$position</h3><pre style='height: 100px;'>$hairpin</pre>";
        my @values2 = @{$values};
        for my $hit (@values2){
#            my $t = $hit->{'name'};
            $contents .= <<"INNER";
<h4>$hit->{'name'}</h4>
<pre>$hit->{'alignment'}</pre>
INNER
        }

    }

    $body = <<"DATA";
	<body onload="displayFile('$url')">
		<div class="titreDiv"> MicroRNA identification results:</div>
		<h2> Alignments :</h2>
		$contents
            <div id = 'alignement'>
    		<pre id = 'preAlign'>
	</div>	</pre>
	</body>
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
