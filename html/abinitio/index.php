<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../style/mirkwood.css' />
        <link type='text/css' rel='stylesheet' href='../../Style/css/bioinfo.css' />
        <link href="/Style/css/page_theme.css" rel="stylesheet" style="text/css"/>
	<script type="text/javascript" src="/scripts/bioinfo.js"></script>
        <script type="text/javascript" src="../../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="../js/header.js"></script>
        <title>miRkwood ab initio</title>
    </head>
    <body>

        <div class="frametitle">
            <h1 id="title" onclick="loadLink('/mirkwood/abinitio/index.php');">miRkwood <em>ab initio</em></h1>
        </div>

        <div id="center_sup">
                       <div class="tabs" id="menu_central" style="display:inline-block">
                <?php include("./header_menu.txt") ?>
            </div>
            <div id="arborescence"></div>
        </div>

        <div id="main">

            <br /> 

<p>
<b>miRkwood ab initio</b> is a web application for the discovery of plant microRNAs in raw
genomic sequences. It is able to face the diversity of plant microRNAs
and searches for both conserved and non-conserved
microRNAs and their precursor. Mirkwood can scan  sequences up to
                                               100,000 nt.
</p>

<br />

<h3>Web server</h3>
<ul>
 <li><a href='/cgi-bin/mirkwood/web_scripts/interface.pl'>click here </a></li>
</ul>

<br />
<h3>Documentation</h3>

            <ul>
                <li><a href='help.php'>read user manual</a></li>
            </ul>

            <br /><br />

      
        </div><!-- main-->

        <?php require("/bio1/www/html/lib.inc")?>
        <?php footer("miRkwood","miRkwood", "mirkwood@univ-lille1.fr","2013"); ?>

    </body>

</html>
