#!/usr/bin/perl -w
use strict;
use warnings;

use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::WebTemplate;

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();

my $css = PipelineMiRNA::WebTemplate->get_css_file();
my $js  = PipelineMiRNA::WebTemplate->get_js_file();

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<meta name="keywords" content="RNA, ARN, mfold, fold, structure, prediction, secondary structure" />
<link title="test" type="text/css" rel="stylesheet" href="$css" />

		<script src="$js" type="text/javascript" LANGUAGE="JavaScript"></script>
	
		
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
