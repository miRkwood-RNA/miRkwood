#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
my $cgi = CGI->new();

my $name     = $cgi->param('nameSeq');
my $factor   = $cgi->param('factor');
my $value    = $cgi->param('value');
my $position = $cgi->param('position');
my $typePage = $cgi->param('typePage');
my $url      = $cgi->param('url');

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
$body
</html>

DATA
    return $html;
}
my $body;
if ( $typePage eq 'simpleCell' ) {
    $body = <<"DATA";
    <body>
        <div class="titreDiv"> MicroRNA identification results:</div>
        <div id = 'showInfo'>
            <h2 ><u>Sequence Informations </u></h2><br/>
            <li><b>Name sequence:</b> ${name}$position
            </li>
            <li><b>$factor:</b> $value
            </li>
        </div>
    </body>
DATA
    my $html = make_HTML($body);
}
elsif ( $typePage eq 'alignement' ) {
    $body = <<"DATA";
	<body onload="displayFile('$url')">
		<div class="titreDiv"> MicroRNA identification results:</div>
		<h2> Alignments :</h2>
    		<pre id = 'preAlign'>
    		<div id = 'alignement'>
	</div>	</pre>
	</body>
DATA
}
elsif ( $typePage eq 'image' ) {
    my @url = split( /\//xms, $url );
    my $image_name = @url[5];
    $body = <<"DATA";
	<body>
		<div class="titreDiv"> MicroRNA identification results:</div>
		<h2> Structure :</h2>
		<div class="figure" >
		  <img src='$url' border=0 alt='image'>
		  <p>Fig : $image_name sequence
    	</div>	
	</body>
DATA
}

print make_HTML($body) or die("Error when displaying HTML: $!");

###End###
