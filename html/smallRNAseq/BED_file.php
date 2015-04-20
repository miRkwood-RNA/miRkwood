<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../../Style/bioinfo.css' />
        <link type='text/css' rel='stylesheet' href='../style/help.css' />
        <link type='text/css' rel='stylesheet' href='../style/rna.css' />
        <script type='text/javascript' src='/js/miARN.js'></script>
        <title>miRkwood - MicroRNA identification - Help</title>
    </head>
    <body>
        <div class="theme-border"></div>
        <a href="/">
        <div class="logo"></div>
        </a>
        <div class="bloc_droit">
            <?php include("./header_menu.txt") ?>
            
            <div class="main-full">

                <p class='title'>How to prepare my input files ?</p>

                <p>The input of miRkwood is a set of reads produced by
                deep sequencing of small RNAs and mapped to a
                reference genome. Typically, length of the reads
                should range between 16nt and 35nt. The user is
                required to upload a BED file that contains all
                positions of mapped sequence tags. This file can be
                obtained from the raw sequencing data by taking three
                easy steps on your computer. If you are new to
                miRkwood, you might also want to test it with a sample
                BED file.</p>
                
                <ul><li> <a href="sample.bed">Download the sample BED file</a> </li></ul>

	      
                <h2 id='adapter'>Remove adapter sequences</h2>
                
                <p>This can be performed with cutadapt, for example.</p>
                
                <pre class='code'>cutadapt -a AACCGGTT -o output.fastq input.fastq</pre>
                
                
                <h2 id='map'>Map the trimmed reads onto the reference genome</h2>
                
                <p>The goal of this step is to generate a BAM file that contains the alignments of the expressed reads with the reference genome.</p>
                
                <p>We recommend to perform exact matching. For that, you can use <a href="http://bowtie-bio.sourceforge.net/index.shtml">Bowtie2</a> with the following parameters. Any other read mapper can also do the job.</p>
                
                <pre class='code'>bowtie -v 0 -f --all --best --strata -S &lt;genome&gt; &lt;reads&gt;   </pre>
                
                <p>Reads file must be in FASTA, FASTQ, or colorspace-fasta format. Genome file must be in FASTA format.</p>
      
                <p>The list of assemblies accepted by miRkwood is given in <a href="./help.php#reference">Section "Select an assembly"</a> on the help page. 
<!--
                Should you work with any other organism or assembly, then you can upload part of the reference genome (<a href="./help.php#reference">Section "Select an assembly"</a> on the help page).
-->
                </p>
                
                
                <h2 id='convert'>Convert the BAM file into a BED file</h2>
                
                <p>For this step, you should use our custom script mirkwood-bam2bed.pl (<a href="../../cgi-bin/mirkwood/web_scripts/getScript.pl?file=mirkwood-bam2bed.pl">download the script</a>).  mirkwood-bam2bed.pl is a perl script dependent upon the installation of SAMtools. In practice, the BED file is up to 10 times smaller than the BAM file and up to XXX times smaller than the set of raw reads, while retaining all information needed to conduct the analysis. This allows to reduce significantly the bandwidth necessary to upload the data to miRkwood server. </p>
                
                <pre class='code'>mirkwood-bam2bed.pl -bam /path/to/your/bam/file -out /output/directory/</pre>
                
                <br />
                
                <p>The generated BED file has the following syntax.</p>
                
                <pre class='example'>
1    18092    18112    SRR051927.5475072    1    -
1    18094    18118    SRR051927.2544175    2    +
1    18096    18119    SRR051927.3033336    1    +
1    18100    18124    SRR051927.172198     9    +
</pre>
                
                <p>In this file, each line is a unique read. The fields are, from left to right:  name of the chromosome, starting position, ending position, read identifier, number of occurrences of the read in the data, strand. Positions follow the BED numbering convention: the first base of the chromosome is considered position 0 (0-based position) and the feature does not include the stop position. </p>
                
                <br />
                
                <p><b>You are now ready to use miRkwood small RNA-seq
                on your data. </b></p>
            
            </div> <!-- main full -->
        
        </div><!-- bloc droit-->
        <?php include("../static/footer.txt") ?>
    </body>

</html>
