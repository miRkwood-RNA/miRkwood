#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use FindBin;                     # locate this script
use lib "$FindBin::Bin/../lib";  # use the parent directory
use PipelineMiRNA::WebTemplate;
use PipelineMiRNA::WebFunctions;
use PipelineMiRNA::Paths;

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();

my $cgi            = CGI->new();
my $jobId          = $cgi->param('jobID');
my $name           = $cgi->param('nameSeq');
my $position       = $cgi->param('position');

my $candidate_name = $name.'__'.$position;
my $job = PipelineMiRNA::WebFunctions->jobId_to_jobPath($jobId);

my %candidate;
my $html_contents;

if (! eval {%candidate = PipelineMiRNA::WebFunctions::retrieve_candidate_information($job, $name, $candidate_name);}) {
    # Catching exception
    $html_contents = "No results for the given identifiers";
}else{

    my $image_url = PipelineMiRNA::Paths->get_server_path($candidate{"image"});
    $html_contents ="
            <div id = 'showInfo'>
        <h2><u>Sequence Informations </u></h2><br/>
        <li>
          <b>Name sequence : </b>$candidate{'name'}
        </li>
        <li>
          <b>MFEI :</b> $candidate{'mfei'}
        </li>
        <li>
          <b>MFE :</b> $candidate{'mfe'}
        </li>
        <li>
          <b>AMFE :</b> $candidate{'amfe'}
        </li>
        <li>
          <b>P_Value :</b> $candidate{'p_value'}
        </li>
        <li>
          <b>Position :</b> $position
        </li>
        <li>
          <b>Self-contain :</b> $candidate{'self_contain'}
        </li>
        <div class='figure' >
          <img src='$image_url' border=0 alt='image'>
          <p>Fig : ${name}__$position sequence</p>
        </div>
    </div><!-- showInfo -->
"
}

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

<html>
	<head>
		<LINK rel="stylesheet" type="text/css" href="/arn/style/script.css" />
		<script src="/arn/js/miARN.js" type="text/javascript" LANGUAGE="JavaScript"></script>
		<title>MicroRNA identification</title>
	</head>
	<body>
        <div class="theme-border"></div>
        <div class="logo"></div>
        $bioinfo_menu
        <div class="bloc_droit">
        $header_menu
        <div class="main">
        $html_contents

	   </div><!-- main -->

$footer

</div><!-- bloc droit-->
	</body>
</html>

DATA
###End###
