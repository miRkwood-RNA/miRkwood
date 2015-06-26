<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../../Style/bioinfo.css' />
        <link type='text/css' rel='stylesheet' href='../style/script.css' />
        <link type='text/css' rel='stylesheet' href='../style/rna.css' />
        <script type='text/javascript' src='/js/miARN.js'></script>
        <title>miRkwood - MicroRNA identification</title>
    </head>
    <body>
        <div class="theme-border"></div>
        <a href="/">
            <div class="logo"></div>
        </a>
        <? include("../static/bioinfo_menu.txt") ?>
        <div class="bloc_droit">
            <? include("header_menu.txt") ?>
            <div class="main">

                <br /><br />

                <p>miRkwood small RNA-seq is a web application that allows to analyse small RNA deep sequencing data in plants. It takes as input a set of short expressed reads (from 15 to 35 nt) that have been previously mapped on a reference genome, and searches for all expressed microRNAs.</p>

<h2> Web server</h2>
	      
                <ul>
		  <li> <a href='/cgi-bin/mirkwood/web_scripts/BAMinterface.pl'>click here</a></li>
</ul>

<h2> Documentation </h2>

<ul>
		  <li><a href='example.php'>quick start with an example</a></li>
		  <li><a href='BED_file.php'>how to prepare the input file</a></li>
		  <li><a href='help.php'>user manual</a></li>  
		</ul>


		<br /><br /><br />

<img style='width:660px; display: block; margin: 0 auto;' src='../style/tree_smallrnaseq.jpg' alt='' />
		
        <br /><br />
            </div><!-- main-->
        </div><!-- bloc droit-->
        <? include("../static/footer.txt") ?>
    </body>

</html>
