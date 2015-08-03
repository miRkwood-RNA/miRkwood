<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../../Style/bioinfo.css' />
        <link type='text/css' rel='stylesheet' href='../style/script.css' />
        <link type='text/css' rel='stylesheet' href='../style/rna.css' />
        <script type='text/javascript' src='/js/miARN.js'></script>
        <title>miRkwood ab initio</title>
    </head>
    <body>
        <div class="theme-border"></div>
            <a href="/">
                <div class="logo"></div>
            </a>
            <? include("../static/bioinfo_menu.txt") ?>
            <div class="bloc_droit">
                <? include("./header_menu.txt") ?>
                <div class="main">

                    <br /> <br />

                    <p>miRkwood <em>ab initio</em> allows to scan a genomic sequence (up to 100,000 nt) and find all potential microRNAs.</p>

                    <p>miRkwood is able to face the diversity of plant miRNAs and allows the prediction of precursors of both conserved and non-conserved microRNAs. It first identifies potential precursors of  microRNAs on the basis of the secondary structure, and then refines this result through a variety of additional complementary features that bring new evidence to the prediction: Thermodynamical stability, nature of complementarity of the miRNA:miRNA* duplex...</p>


		  
                    <ul>
                        <li><a href='help.php'>read user manual</a></li>
                        <!-- <li><a href='method.php'>learn about the method</a></li> -->
                        <li><a href='/cgi-bin/mirkwood/web_scripts/interface.pl'>use web server</a></li>
                    </ul>

	<br /> <br /><br /><br />

<img style='width:660px; display: block; margin: 0 auto;' src='../style/pipeline_abinitio.jpg' alt='' />
		    
                </div><!-- main-->
            </div><!-- bloc droit-->
    <? include("../static/footer.txt") ?>
    </body>

</html>
