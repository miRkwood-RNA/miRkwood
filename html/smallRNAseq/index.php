<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../style/mirkwood.css' />
        <link type='text/css' rel='stylesheet' href='../../Style/css/bioinfo.css' />
    <script type="text/javascript" src="/scripts/bioinfo.js"></script>
        <script type="text/javascript" src="../../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="../../libs/jquery.history.js"></script>
        <script type="text/javascript" src="../js/header.js"></script>
        <title>miRkwood small RNA-seq</title>
    </head>
    <body>
        <div class="frametitle">
            <h1 id="title" onclick="loadLink('/mirkwood/smallRNAseq/index.php');">miRkwood small RNA-seq</h1>
        </div>

        <div class="tabs" id="menu_central" style="display:inline-block">
            <?php include("./header_menu.txt") ?>
        </div>
        <div id="arborescence"></div>

        <div id="main">

            <br />

<p>
<b>miRkwood small RNA-seq</b> is a comprehensive toolbox  for the analysis of small RNA deep
sequencing data in plants.
It takes as input a set of short expressed reads (from 15 to 35 nt)
that have been previously mapped on a reference genome, and searches for all microRNAs present in the data.
Furthermore, it provides additional features that help the user to make the most out of the data:
thermodynamic stability, quality of the miRNA-miRNA* duplex, read coverage, read visualisation, miRNA conservation,...
</p>

<br />

            <h3> Web server</h3>

            <ul>
                <li> <a href='/cgi-bin/mirkwood/web_scripts/BAMinterface.pl'>click here</a></li>
            </ul>

            <h3> Documentation </h3>

            <ul>
                <li><a href='example.php'>quick start with an example</a></li>
                <li><a href='BED_file.php'>how to prepare the input file</a></li>
                <li><a href='help.php'>user manual</a></li>
            </ul>

            <h3> Local install </h3>

            <ul>
                <li>see <a href='https://github.com/miRkwood-RNA/miRkwood'> GitHub</a></li>
            </ul>

    
            <br /><br />

        </div>
        <?php include('../static/footer.php') ?>
    </body>

</html>
