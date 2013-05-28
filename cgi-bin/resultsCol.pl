#!/usr/bin/perl -w
use CGI;

my $cgi = new CGI; 
$factor = $cgi->param('factor');
@values = $cgi->param();


print <<DATA;
Content-type: text/html

<html>
	<head>
		<LINK rel="stylesheet" type="text/css" href="/arn/css/script.css" />
		
		<script src="/arn/js/miARN.js" type="text/javascript" LANGUAGE="JavaScript"></script>
		
		<title>MicroRNA identification</title>
	</head>
	<body>
		<div class="titreDiv"> MicroRNA identification results:</div>
	
		<div id = 'showInfo'>
	
		 <h2 class='titre'><u>List of 
DATA

print $factor;
print <<DATA;
s</u></h2>

DATA

	@names = split(/,/,$cgi->param(@values[3]));
	@values = split(/,/,$cgi->param(@values[1]));
	for ($i=0;$i<scalar(@values);$i++)
	{
		print "<li><b>".@names[$i].'</b> : '.@values[$i].'</li>';
	}
	
print <<DATA;


</div>
	</body>
</html>

DATA
###End###
