#!/usr/bin/perl -w
use strict;
use warnings;

use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory

use PipelineMiRNA::WebTemplate;

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html


<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<meta name="keywords" content="RNA, ARN, mfold, fold, structure, prediction, secondary structure" />
<link title="test" type="text/css" rel="stylesheet" href="/arn/style/script.css" />

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

<br>

<h2> Retrieve result with an ID</h2>

<br>

<div class="forms">
<form method="post" action="./resultsWithID.pl">
<p>The ID remains valid 24 hours after sequence submission.</p>
<p>

<b>  Enter the ID :</b>  
<input type="text" name="run_id" size="20">
<input type="hidden" name="command" value="result">
<input type="submit" value="Go"></p>
</form>
</div> <!-- form -->
<br>
</div>

$footer

</body>
</html>

DATA
###End###
