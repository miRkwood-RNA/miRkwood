<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../../Style/css/bioinfo.css' />
        <link type='text/css' rel='stylesheet' href='../style/mirkwood.css' />
        <link href="/Style/css/page_theme.css" rel="stylesheet" style="text/css"/>

        <script type="text/javascript" src="../../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="../js/header.js"></script>
        <title>miRkwood ab initio</title>
    </head>
    <body>

        <div class="frametitle">
            <h1 id="title">miRkwood <em>ab initio</em></h1>
        </div>

        <div id="center_sup">
            <div id="link_home" style="display:inline-block"><a href="../index.php" class="text_onglet"><img src="/Style/icon/home_w.png" alt="home_general"/></a></div>
            <div class="tabs" id="menu_central" style="display:inline-block">
                <?php include("./header_menu.txt") ?>
            </div>
            <div id="arborescence"></div>
        </div>

        <div id="main">

            <br /> <br />

            <p>miRkwood <em>ab initio</em> allows to scan a genomic sequence (up to 100,000 nt) and find all potential microRNAs.</p>

            <p>miRkwood is able to face the diversity of plant miRNAs and allows the prediction of precursors of both conserved and non-conserved microRNAs. It first identifies potential precursors of  microRNAs on the basis of the secondary structure, and then refines this result through a variety of additional complementary features that bring new evidence to the prediction: Thermodynamical stability, nature of complementarity of the miRNA:miRNA* duplex...</p>

            <ul>
                <li><a href='help.php'>read user manual</a></li>
                <li><a href='/cgi-bin/mirkwood/web_scripts/interface.pl'>use web server</a></li>
            </ul>

            <br /><br />

            <img style='width:660px; display: block; margin: 0 auto;' src='../style/pipeline_abinitio.jpg' alt='' />

        </div><!-- main-->

        <?php require("/bio1/www/html/lib.inc")?>
        <?php footer("miRkwood","miRkwood", "mirkwood@univ-lille1.fr","2013"); ?>

    </body>

</html>
