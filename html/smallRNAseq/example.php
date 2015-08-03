<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../../Style/bioinfo.css' />
        <link type='text/css' rel='stylesheet' href='../style/help.css' />
        <link type='text/css' rel='stylesheet' href='../style/rna.css' />
        <script type='text/javascript' src='../js/miARN.js'></script>
        <title>miRkwood small RNA-seq - Run an example</title>
    </head>
    <body>
        <div class="theme-border"></div>
        <a href="/">
            <div class="logo"></div>
        </a>
        <div class="bloc_droit">
            <?php include("./header_menu.txt") ?>
            
            <div class="main-full">
        
                <p class='title'>Quick start</p>
	      
                <p>This page covers the basics on how to get started
                with miRkwood small RNA-seq. The main input of
                miRkwood is a set of reads (ranging from 15 nt to 35
                nt approximately) that have been previously mapped on
                a reference genome and that are stored in a file in
                BED format. We advise you to first practice with the
                sample BED file below. </p>
                
                <ul><li> <a href="sample.bed">Download the sample BED file</a> </li></ul> 
                 
                <p>This file contains Illumina reads for <i> Arabidopsis thaliana</i> from SRA <a href="http://www.ncbi.nlm.nih.gov/sra/?term=SRR051927">SRR051927</a>. We have mapped these reads on <i>TAIR 10</i>, and have selected the portion of reads corresponding to the first 2,200,00 nt of chromosome 1. This gives a total of 19,756 reads. </p>
                
                <p>If you want to know more on how to create this BED file for your data, check out the <a href="BED_file.php">detailed instruction</a>.</p> 

                <h2> Input form</h2> 

                <p>You can now open the <a href="/cgi-bin/mirkwood/web_scripts/BAMinterface.pl"> input form</a> on the web server.</p>

                <p><b>Upload your set of reads: </b> choose the downloaded sample BED file on your local disk.</p>
                
                <p><b>Select a species: </b> choose the <i> Arabidospsis thaliana TAIR 10</i> assembly.</p>
                
                <p><b>Parameters: </b> keep all default options. </p>
                
                <p>Finally, click the <i>Run miRkwood</i> button. Results will be displayed in a couple of minutes. </p>

                <h2> Results page </h2>

		<p>The top of the page displays the job ID, that can
		be saved for future usage (with the link <i>retrieve a
		result with an ID</i>, on the main menu).</p> 
		
                <p>The result page has then two main parts. The first one
                (<i>Options summary</i>) is simply a summary of your
                job parameters. The other one (<i>Results summary</i>)
                provides the detailed results.</p> 


	<img style='width:700px; display: block; margin: 0
	auto;'src='../style/repartition.png' alt='read graphics' />

		<p> This diagram indicates the proportion of reads
		found for each category. This is a graphical
		visualisation of the results displayed below
		</p>	
                <ul>
                    <li>
                        <i>Total number of reads:</i> 19756 (4075 unique reads) <br> <br>
                        The initial BED file contains 19756 reads, forming 4075 unique reads.
                        <br> <br> 
                    </li> 

                    <li>
                        <i>CoDing Sequences:</i> 680 reads <br> <br>
                        680 reads fall within a coding region (annotated as CDS) and are
                        discarded from the analysis.  You can list them by clicking on
                        the <i>download</i> link (GFF file).
                        <br> <br>
                    </li>
                    
                    <li>
                        <i>rRNA/tRNA/snoRNA:</i> 11568 reads <br> <br>
                        11568 reads fall within a region annotated either as ribosomal RNA, or
                        transfer RNA or snoRNA, and are discarded from the analysis. 
                        You can list them by clicking on
                        the <i>download</i> link (GFF file).
                        <br> <br>
                    </li> 

                    <li>
                        <i>Multiply mapped reads:</i> 16 reads <br> <br> 
                        16 reads map to more than five locations, and are discarded from the
                        analysis.  You can list them by clicking on
                        the <i>download</i> link (BED file).
                        <br> <br>
                    </li> 


		    <li>
		    <i>Orphan cluster of reads:</i> 3562 reads<br><br>
		    A cluster of reads is a short region in the genome
		    that has been enriched with aligned reads. Here we
		    report the number of reads that occur in a cluster not classified as
		    miRNA by miRkwood, due to the secondary structure.  You
		    can obtain the list of the corresponding orphan clusters by clicking on
                    the <i>download</i> link (BED file).
		    </li>
            <br />
		    <li><i> Unclassified reads:</i> 2464 reads<br><br>
		    2464 reads do not belong to any cluster, or do not
		    fall in any annotated region.
		    </li>
		    <br />
                    <li> 
                        <i>Known miRNAs:</i> 5 <br> <br>
                        5 precursors of miRNA present in miRbase
                        intersect with a total of 346 reads
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
                        <i> Novel miRNAs: </i> 21 <br> <br>
                        miRkwood finds 21 new miRNAs that have not
                        been previously  reported in miRbase and that are supported by a significant
                        read coverage  and a stemloop secondary
                        structure. These 21 miRNAs represent a total
                        of 1120 reads.<br> <br> 
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
 formats, visualisation tools,  please visit our  main <a
 href="help.php">help page</a>.</p>
		
		<center>

		<img style='width:700px; display: block; margin: 0 auto;'src='../style/table_novelmirna.png' alt='results table' />
        
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

            </div><!-- main-->
        </div><!-- bloc droit-->

        <? include("../static/footer.txt") ?>
    </body>

</html>
