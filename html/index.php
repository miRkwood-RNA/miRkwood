<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
  <head>
        <link type='text/css' rel='stylesheet' href='./Style/bioinfo.css' />
		<link type='text/css' rel='stylesheet' href='./style/script.css' />
		<link type='text/css' rel='stylesheet' href='./style/rna.css' />
        <script type='text/javascript' src='/js/miARN.js'></script>
        <title>MicroRNA identification</title>
    </head>
 <body>
        <div class="theme-border"></div>
        <div class="logo"></div>
         <? include("./static/bioinfo_menu.txt") ?>
        <div class="bloc_droit">
        <? include("./static/header_menu.txt") ?>
<div class="main">
         <p>  miRkwood is a computational pipeline for the identification of plant miRNAs and their hairpin precursors. It can be used with short  genomicregions or assembled expressed sequences (Sanger or NGS).  </p>
<ul>
 <li>see user manual</li>
 <li>see an example</li>
 <li>use web server</li>
</ul>
<p>The method is is specifically calibrated for plant miRNAs, that have a large structural diversity. It first identifies potential pre-miRNAs on the basis of the secondary structure, and then refines this result through a variety of additional complementary features that bring new evidence to the prediction:  Thermodynamical stability, conservation of the mature miRNA, nature of complementarity of the miRNA-miRNA* duplex. </p>
 </div>
        </div><!-- bloc droit-->
       <? include("./static/footer.txt") ?>
    </body>
    
</html>
