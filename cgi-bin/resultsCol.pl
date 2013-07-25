#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;

my $cgi        = CGI->new();
my $factor     = $cgi->param('factor');
my @CGI_values = $cgi->param();

my @names  = split( /,/xms, $cgi->param( $CGI_values[3] ) );
my @values = split( /,/xms, $cgi->param( $CGI_values[1] ) );

my $HTML_list;
for ( my $i = 0 .. scalar(@values) ) {
    $HTML_list .= '<li><b>' . $names[$i] . '</b> : ' . $values[$i] . '</li>';
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
		<div class="titreDiv"> MicroRNA identification results:</div>
		<div id = 'showInfo'>
		  <h2 class='titre'>List of ${factor}s</h2>
            $HTML_list
        </div>
	</body>
</html>

DATA
###End###
