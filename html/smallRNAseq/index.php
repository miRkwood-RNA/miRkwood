<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../../Style/css/bioinfo.css' />
        <link type='text/css' rel='stylesheet' href='../style/mirkwood.css' />
<!--
        <link href="/Style/css/page_theme.css" rel="stylesheet" style="text/css"/>
-->
        <script type="text/javascript" src="../../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="../../libs/jquery.history.js"></script>
        <script type="text/javascript" src="../../scripts/bioinfo.js"></script>
        <title>miRkwood small RNA-seq</title>
    </head>
    <body>
        <div class="frametitle">
            <h1 id="title">miRkwood small RNA-seq</h1>
        </div>

        <div id="center_sup">
            <div id="link_home" style="display:inline-block"><a href="../index.php" class="text_onglet"><img src="/Style/icon/home_w.png" alt="home_general"/></a></div>
            <div class="tabs" id="menu_central" style="display:inline-block">
                <? include("./header_menu.txt") ?>
            </div>
            <div id="arborescence"></div>
        </div>

        <div id="main">

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

        </div>
        <?php require("../lib.inc")?>
        <?php footer("miRkwood","miRkwood", "mirkwood@univ-lille1.fr","2013"); ?>
    </body>

</html>
