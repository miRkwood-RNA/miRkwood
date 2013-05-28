#!/usr/bin/perl -w
use CGI;
my $cgi = new CGI; 

$name = $cgi->param('nameSeq');
$factor = $cgi->param('factor');
$value = $cgi->param('value');
$position = $cgi->param('position');
$typePage = $cgi->param('typePage');
$url = $cgi->param('url');
if ($typePage eq 'simpleCell')
{

print <<DATA;
Content-type: text/html

<html>
	<head>
		<LINK rel="stylesheet" type="text/css" href="/arn/css/script.css" />
		<link href="/arn/css/basic.css" type="text/css" rel="stylesheet" />
		<script type="text/javascript" language="Javascript" src="/arn/js/results.js"> </script>
		<title>MicroRNA identification</title>
	</head>
	<body>
		<div class="titreDiv"> MicroRNA identification results:</div>

		<div id = 'showInfo'>
		<h2 ><u>Sequence Informations </u></h2><br/>
		<li><b>Name sequence :</b>
DATA
print $name."_".$position;

print <<DATA;
</li>
<li> <b>
DATA
print $factor;
print <<DATA;
</b>:
DATA
print $value;

print <<DATA;
</li>

	</div>	
	</body>
</html>

DATA
}
if ($typePage eq 'alignement')
{

print <<DATA;
Content-type: text/html

<html>
	<head>
		<LINK rel="stylesheet" type="text/css" href="/arn/css/script.css" />
		
		<script type="text/javascript" src="/arn/js/alignement.js"></script>
			
		<script>
		 
DATA

print "var url ='".$url."';";
print <<DATA;

		</script>
		<title>MicroRNA identification</title>
	</head>
	<body onload="displayFile(url)">
		<div class="titreDiv"> MicroRNA identification results:</div>
		<h2> Alignments :</h2>
		<pre id = 'preAlign'>
		<div id = 'alignement'>
		

	</div>	</pre>
	</body>
</html>

DATA
}

if ($typePage eq 'image')
{

print <<DATA;
Content-type: text/html

<html>
	<head>
		<LINK rel="stylesheet" type="text/css" href="/arn/css/script.css" />
		<script type="text/javascript" language="Javascript" src="/arn/js/results.js"> </script>
		<title>MicroRNA identification</title>
	</head>
	<body >
		<div class="titreDiv"> MicroRNA identification results:</div>
		<h2> Structure :</h2>
	
		<div class="figure" >
DATA
		print "<p><img src='".$url."' border=0 alt='image'>";
		@url = split(/\//,$url);
		$name = $url[5];
		print "<p>Fig : ".$name." sequence";
		
print <<DATA;
	</div>	
	</body>
</html>

DATA
}
###End###
