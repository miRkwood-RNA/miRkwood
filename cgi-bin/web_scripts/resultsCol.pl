#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }

my $cgi        = CGI->new();
my $factor     = $cgi->param('factor');
my @CGI_values = $cgi->param();

my @names  = split( /,/xms, $cgi->param( $CGI_values[3] ) );
my @values = split( /,/xms, $cgi->param( $CGI_values[1] ) );

my $css = PipelineMiRNA::WebTemplate->get_css_file();
my $js  = PipelineMiRNA::WebTemplate->get_js_file();

my $HTML_list;
for ( my $i = 0 .. scalar(@values) ) {
    $HTML_list .= '<li><b>' . $names[$i] . '</b> : ' . $values[$i] . '</li>';
}

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

<html>
	<head>
		<LINK rel="stylesheet" type="text/css" href="$css" />
		
		<script src="$js" type="text/javascript" LANGUAGE="JavaScript"></script>
		
		<title>MicroRNA identification</title>
	</head>
	<body>
		<div class="titreDiv"> MicroRNA identification results:</div>
		<div id = 'showInfo'>
		  <h2 class='titre'>List of ${factor}s</h2>
            $HTML_list
        </div>
	</body>
</html>

DATA
###End###
