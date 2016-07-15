<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../style/mirkwood.css' />
        <link href="/Style/css/page_theme.css" rel="stylesheet" style="text/css"/>
        <link type='text/css' rel='stylesheet' href='../style/help.css' />
        <link type='text/css' rel='stylesheet' href='../style/rna.css' />
        <link type='text/css' rel='stylesheet' href='/Style/css/bioinfo.css' />
	<script type="text/javascript" src="/scripts/bioinfo.js"></script>
        <script type="text/javascript" src="../../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="../../libs/jquery.history.js"></script>
        <script type="text/javascript" src="../js/header.js"></script>
        <title>miRkwood small RNA-seq - Run an example</title>
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

           
<div id="mirkwood_title">Quick start</div>

            <p>This page covers the basics on how to get started
                                               with miRkwood small RNA-seq. A comprehensive user guide is to be found <a href="./help.php">here</a>.</p>

<p>The main input of
            miRkwood is a set of reads (ranging from 15 nt to 35
            nt approximately) that have been previously mapped on
            a reference genome and that are stored in a file in
            BED format. We advise you to first practice with the
            sample BED file below. </p>

            <ul><li> <a href="sample.bed">Download the sample BED file</a> </li></ul> 

            <p>This file contains Illumina reads for <i> Arabidopsis thaliana</i> from SRA <a href="http://www.ncbi.nlm.nih.gov/sra/?term=SRR960237 ">SRR960237 </a>. We have mapped these reads on <i>TAIR 10</i>, and have selected 2 portions of reads, one from chromosome 1 and one from chromosome 4. This gives a total of 319,103 reads. </p>

            <p>If you want to know more on how to create this BED file for your data, check out the <a href="BED_file.php">detailed instruction</a>.</p> 


            <h2> Input form</h2> 

            <p>You can now open the <a href="/cgi-bin/mirkwood/web_scripts/BAMinterface.pl"> input form</a> on the web server.</p>

            <p><b>Upload your set of reads: </b> choose the downloaded sample BED file on your local disk.</p>

            <p><b>Select a species: </b> choose the <i> Arabidospsis thaliana TAIR 10</i> assembly.</p>

            <p><b>Parameters: </b> keep all default options. </p>

            <p>Finally, click the <i>Run miRkwood</i> button. Results will be displayed in a couple of minutes. </p>


            <h2> Results page </h2>

            <p>The top of the page displays the job ID, that can
            be saved for future usage (with the link <a href="../id.php">retrieve a
            result with an ID</a>, on the main menu).</p> 

            <p>The result page has then two main parts. The first one
            (<i>Options summary</i>) is simply a summary of your
            job parameters. The other one (<i>Results summary</i>)
            provides the detailed results.</p> 


            <img style='width:700px; display: block; margin: 0
            auto;'src='../style/reads_repartition.png' alt='read graphics' />

            <p> This diagram indicates the proportion of reads
            found for each category. This is a graphical
            visualisation of the results displayed below.</p>

            <ul>
                <li>
                    <i>Total number of reads:</i> 319,103 (29,795 unique reads) <br> <br>
                    The initial BED file contains 319,103 reads, forming 29,795 unique reads.
                    <br> <br>
                </li>

                <li>
                    <i>CoDing Sequences:</i> 2,599 reads <br> <br>
                    2,599 reads fall within a coding region (annotated as CDS) and are
                    discarded from the analysis.  You can list them by clicking on
                    the <i>download</i> link (GFF file).
                    <br> <br>
                </li>

                <li>
                    <i>rRNA/tRNA/snoRNA:</i> 6,275 reads <br> <br>
                    6,275 reads fall within a region annotated either as ribosomal RNA, or
                    transfer RNA or snoRNA, and are discarded from the analysis. 
                    You can list them by clicking on
                    the <i>download</i> link (GFF file).
                    <br> <br>
                </li>

                <li>
                    <i>Multiply mapped reads:</i> 1,665 reads <br> <br> 
                    1,665 reads map to more than five locations, and are discarded from the
                    analysis.  You can list them by clicking on
                    the <i>download</i> link (BED file).
                    <br> <br>
                </li>

                <li>
                    <i>Orphan cluster of reads:</i> 26,826 reads<br><br>
                    A cluster of reads is a short region in the genome
                    that has been enriched with aligned reads. Here we
                    report the number of reads that occur in a cluster not classified as
                    miRNA by miRkwood, due to the secondary structure.  You
                    can obtain the list of the corresponding orphan clusters by clicking on
                    the <i>download</i> link (BED file).
                    <br> <br>
                </li>

                <li>
                    <i>Orphan hairpins:</i> 277 reads<br><br>
                    An orphan hairpin is a candidate with a global score of 0
                    and showing no conservation with miRBase.
                    If you select the option "filter out low quality hairpins",
                    low quality candidates will be discarded and you can obtain the list 
                    of the corresponding orphan hairpins by clicking on the 
                    <i>download</i> link (BED file).
                    <br> <br>
                </li>

                <li>
                    <i> Unclassified reads:</i> 49,706 reads<br><br>
                    49,706 reads do not belong to any cluster, or do not
                    fall in any annotated region.
                    <br> <br>
                </li>

                <li> 
                    <i>Known miRNAs:</i> 5 <br> <br>
                    5 precursors of miRNA present in miRbase
                    intersect with a total of 224,134 reads
                    present in the input BED file.
                    <br> <br> 
                    You can display detailed results by clicking on the link <i>see results</i>.
                    From this new page, you will be able to access a short report for each miRNA found :
                    miRBase accession number, positions of the
                    precursor, sequence of the miRNA, secondary structure, read
                    distribution. Further information can then be retrieved from
                    the miRBase website.
                    <br> <br>
                </li>

                <li> 
                    <i> Novel miRNAs: </i> 30 <br> <br>
                    miRkwood finds 30 new miRNAs that have not
                    been previously  reported in miRbase and that are supported by a significant
                    read coverage  and a stemloop secondary
                    structure. These 30 miRNAs represent a total
                    of 7,621 reads.<br> <br> 
                    You can display detailed results by clicking on the link <i>see results</i>.
                    From this new page, you will be able to access a comprehensive report
                    for each predicted miRNA :
                    positions of the precursor, sequence of the
                    miRNA, secondary structure, read
                    distribution, thermodynamic stability, precision of the duplex
                    processing, conservation, ... 
                    Predictions are ranked according to a
                    quality score.
                    <br /> <br />
                </li> 

            </ul>

            <p>If you want to know more about all parameter options, export
            formats, visualisation tools,  please visit our main 
            <a href="help.php">help page</a>.</p>

            <center>

                <img style='width:700px; display: block; margin: 0 auto;'src='../style/results_novelmirna.png' alt='results table' />

                <br />

<pre class='example'>
Locus  : 1:234009-234092
Strand : -

GAAAUGAUGCGCAAAUGCGGAUAUCAAUGUAAAUCAGGGAGAAGGCAUGAUAUACCUUUAUAUCCGCAUUUGCGCAUCAUCUCU
((.(((((((((((((((((((((.((.(((.((((.(.......).)))).))).)).))))))))))))))))))))).)).
         <------miRBase------>                          <------miRBase------>
*********************............................................................... length=21 depth=5
.......*********************........................................................ length=21 depth=2
........................................................*********************....... length=21 depth=1
............................................................*********************... length=21 depth=16
</pre>

                <br />

                <img style='width:600px; display: block; margin: 0 auto;' src='../style/hairpin_with_mature.png' alt='hairpin with mature' />		

            </center>

            <br />

        </div>

        <?php require("/bio1/www/html/lib.inc")?>
        <?php footer("miRkwood","miRkwood", "mirkwood@univ-lille1.fr","2013"); ?>  
    </body>

</html>
