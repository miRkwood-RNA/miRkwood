#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
my $cgi          = CGI->new();
my $name         = $cgi->param('nameSeq');
my $pvalue       = $cgi->param('pvalue');
my $position     = $cgi->param('position');
my $mfei         = $cgi->param('mfei');
my $mfe          = $cgi->param('mfe');
my $amfe         = $cgi->param('amfe');
my $self_contain = $cgi->param('self_contain');
my $Vienna       = $cgi->param('Vienna');
my $DNASequence  = $cgi->param('DNASequence');
my $image        = $cgi->param('image');
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
		<h2><u>Sequence Informations </u></h2><br/>
		<li>
		  <b>Name sequence :</b> $name
        </li>
        <li>
          <b>MFEI :</b>$mfei
        </li>
        <li>
          <b>MFE :</b>$mfe
        </li>
        <li>
          <b>AMFE :</b>$amfe
        </li>
        <li>
          <b>P_Value :</b>$pvalue
        </li>
        <li>
          <b>Position :</b>$position;
        </li>
        <li>
          <b>Self-contain :</b>$self_contain;
        </li>
        <li>
          <b>Vienna :</b>
        </li>
        <pre>
>${name}__$position
$Vienna
$DNASequence
</pre>
        <li>
          <b>Structure :</b>
        </li>
        <a href="../arn/programs/varna.jnlp">Launch</a>
		<div class="figure" >
		  <img src='$image' border=0 alt='image'>";
		  <p>Fig : ${name}__$position sequence</p>
    	</div>	
	</div>	
	</body>
</html>

DATA
###End###
