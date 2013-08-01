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

<div id="page">


<h2>Error ! </h2>

<div class="section">
MiREST could not recognize the format of one or more sequences.
</div>


<ul>
<li> Your data set should include <b>at least two</b> distinct sequences.</li>
</ul>
<ul>
<li> Sequences should be in <b>FASTA</b> format.
A sequence in FASTA format consists of a
        single-line description, followed by lines of sequence
        data. The first character of the description line is a
        greater-than (&quot;&gt;&quot;) symbol in the first
        column. All lines should be shorter than 80 characters. An
        example  in FASTA format is: 

<pre>
> Name of the sequence 1
ctgcgagcgcgcgatgatagcgcggcgagcatgtagcatgctagctgtcgcgagcact
cggccgagatcaggcgatgcatgcgcagggagcagcgagcgacgagcacagcatgctagctagatgcatgctgtaggcagc
cgccgagagacgatggagctgc
> Name of the sequence 2
gacagatacgataagaggacgggatagaacgtagacatcgccgagagacgatggagctgc
cggccgagatcaggcgatgcatgcgcagggagcaggcgagcatgtagcatgctagctgtcgcgagcact
</pre></li>
</ul>
<ul>
<li>  <b>Lower-case and upper-case</b> letters are both
        accepted.</li>
</ul>
<ul>
<li>  The full  standard <b>IUPAC</b> nucleic acid code is not supported : only <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt> and <tt>U</tt> symbols are recognized.</li>
</ul>
<ul>
<li>  <b>Numerical digits <tt>0</tt>, ..., <tt>9</tt>, <tt>-</tt> and dot <tt>.</tt></b> symbols are accepted. They are simply ignored by MiREST.</li>

</ul>

<br /><br />
<div class="home-made-center">
  <input type="button" name="nom" value="Submit a new job"
onclick="window.location.href='../../carnac/carnac.html'" />

</div>
</div>
</div>
$footer

</div><!-- bloc droit-->

</body>
</html>

DATA
###End###
