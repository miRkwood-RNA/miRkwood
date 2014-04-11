<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
  <head>
        <link type='text/css' rel='stylesheet' href='./Style/bioinfo.css' />
		<link type='text/css' rel='stylesheet' href='./style/script.css' />
		<link type='text/css' rel='stylesheet' href='./style/rna.css' />
        <script type='text/javascript' src='/js/miARN.js'></script>
        <title>miRkwood - MicroRNA identification</title>
    </head>
 <body>
<div class="theme-border"></div>
        <div class="logo"></div>
         <? include("./static/bioinfo_menu.txt") ?>
        <div class="bloc_droit">
        <? include("./static/header_menu.txt") ?>
<div class="main">

    <p>miRkwood is a computational pipeline for the identification of miRNAs and their hairpin precursors. The method is specifically calibrated for plant miRNAs, that have a large structural diversity. It can be used with short (&#60;100 kb) genomic regions or assembled expressed sequences (Sanger or NGS).</p>

    <p>miRkwood first identifies potential pre-miRNAs on the basis of the secondary structure, and then refines this result through a variety of additional complementary features that bring new evidence to the prediction:  Thermodynamical stability, conservation of the mature miRNA, nature of complementarity of the miRNA-miRNA* duplex.</p>


<ul>
 <li><a href='help.php'>read user manual</a></li>
 <li><a href='method.php'>learn about the method</a></li>
 <li><a href='/cgi-bin/interface.pl'>use web server</a></li>
</ul>

<h2>Authors</h2>
<p>Sylvain Legrand, <a href='http://pdv.univ-lille1.fr/labo/sitesadv/'>SADV</a> (University of Lille and INRA)<br/>
Mohcen Benmounah, Jean-Fr&eacute;d&eacute;ric Berthelot, H&eacute;l&egrave;ne Touzet, <a href='http://www.lifl.fr/bonsai'>Bonsai</a> (LIFL and Inria Lille)</p>

<h2>Contact</h2>
 </div>
        </div><!-- bloc droit-->
       <? include("./static/footer.txt") ?>
    </body>
    
</html>
