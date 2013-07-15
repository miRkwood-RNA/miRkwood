#!/usr/bin/perl -w
use CGI;
my $cgi = new CGI; 
$name = $cgi->param('nameSeq');
$pvalue = $cgi->param('pvalue');
$position = $cgi->param('position');
$mfei = $cgi->param('mfei');
$mfe = $cgi->param('mfe');
$amfe = $cgi->param('amfe');
$self_contain = $cgi->param('self_contain');
$Vienna = $cgi->param('Vienna');
$DNASequence = $cgi->param('DNASequence');
$image = $cgi->param('image');
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
		<h2><u>Sequence Informations </u></h2><br/>
		<li><b>Name sequence :</b>
DATA
print $name;
print <<DATA;
</li>
<li><b>MFEI :</b>
DATA
print $mfei;
print <<DATA;
</li>
<li><b>MFE :</b>
DATA
print $mfe;
print <<DATA;
</li>
<li><b>AMFE :</b>
DATA
print $amfe;
print <<DATA;
</li>
<li><b>P_Value :</b>
DATA
print $pvalue;
print <<DATA;
</li>
<li><b>Position :</b>
DATA
print $position;
print <<DATA;
</li>
<li><b>Self-contain :</b>
DATA
print $self_contain;
print <<DATA;
</li>
<li><b>Vienna :</b></li>
DATA
print "<pre>>".$name.'__'.$position."\n".$Vienna."</pre>";
print <<DATA;

<li><b>Sequence :</b></li>
DATA
print "<pre>>".$name.'__'.$position."\n".$DNASequence."</pre>";
print <<DATA;


<li><b>Structure :</b></li>
<a href="../arn/programs/varna.jnlp">Launch</a>
		<div class="figure" >
DATA

		print "<p><img src='".$image."' border=0 alt='image'>";
		print "<p>Fig : ".$name."__".$position." sequence";

print <<DATA;
	</div>	
	</div>	
	</body>
</html>

DATA
###End###
