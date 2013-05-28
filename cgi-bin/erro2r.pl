#!/usr/bin/perl -w

print <<DATA
Content-type: text/html

<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<link title="test" type="text/css" rel="stylesheet" href="style_bioinfo2.css" />
<title>CARNAC : error message</title>
</head>
<body>
<div id="surpage"> 
<div id="title">
CARNAC
</div>

<div class="menu">
       <a class="menu" href="index.html">home</a>
       <a class="selectionne" href="carnac.html">web server</a>
       <a class="menu" href="help.html">help</a>
       <a class= "menu" href="examples.html">examples</a>
       <a class= "menu" href="id.html">retrieve result with an ID</a>
</div>

<div id="page">


<h2>Error ! </h2>

<div class="section">
Carnac could not recognize the format of one or more sequences.
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
<li>  <b>Numerical digits <tt>0</tt>, ..., <tt>9</tt>, <tt>-</tt> and dot <tt>.</tt></b> symbols are accepted. They are simply ignored by Carnac.</li>
  
</ul>

<br /><br />
<div class="home-made-center">
  <input type="button" name="nom" value="Submit a new job"
onclick="window.location.href='../../carnac/carnac.html'" />

</div>
</div>
<div class="copyright">
For questions about CARNAC or for bug reports, please contact    <a href="mailto: carnac@lifl.fr">carnac@lifl.fr</a><br /> 
Copyright &copy; 2003 Laboratoire d'Informatique Fondamentale de Lille, Universit&eacute; des Sciences et Technologies de Lille</div>
</div>
</body>
</html>
DATA
