<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../style/mirkwood.css' />
        <link type='text/css' rel='stylesheet' href='../style/help.css' />
        <link type='text/css' rel='stylesheet' href='../style/rna.css' />
        <link type='text/css' rel='stylesheet' href='/Style/css/bioinfo.css' />
        <script type="text/javascript" src="/scripts/bioinfo.js"></script>
        <script type="text/javascript" src="../../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="../../libs/jquery.history.js"></script>
        <script type="text/javascript" src="../js/header.js"></script>
        <title>miRkwood small RNA-seq - input file</title>
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

            <div id="mirkwood_title">How to prepare my input files ?</div>

            <p>The input of miRkwood is a set of reads produced by
            deep sequencing of small RNAs and mapped to a
            reference genome. Typically, length of the reads
            should range between 15nt and 35nt. The user is
            required to upload a BED file that contains all
            positions of mapped sequence tags. This file can be
            obtained from the raw sequencing data by taking three
            easy steps on your computer. If you are new to
            miRkwood, you might also want to test it with the sample
            BED file provided below.</p>

            <ul><li> <a href="sample.bed">Download the sample BED file</a> </li></ul>

            <h2 id='adapter'>Remove adapter sequences</h2>

            <p>This can be performed with cutadapt, for example.</p>

            <pre class='code'>cutadapt -a AACCGGTT -o output.fastq input.fastq</pre>

            <h2 id='prinseq'>Run quality control</h2>

            <p>The aim of this step is to filter too short or too
            long sequences and to remove or to trim the low
            quality sequences. This can be achieved using prinseq,
            with this command line as example.</p>

<pre class='code'> prinseq-lite.pl -fastq &lt;short_reads_file.fastq&gt; -min_len 18 -max_len 25 -noniupac
-min_qual_mean 25 -trim_qual_right 20 -ns_max_n 0</pre>

            <p>In this command line we keep only the sequences between 18 and 25 nt with a mean
            quality of at least 25 (phred score) and composed of nucleotides ACGT.
            The sequences are trimmed by quality score from the 3'-end with
            a value of 20 as threshold.</p>


            <h2 id='map'>Map the trimmed reads on the reference genome</h2>

            <p>The goal of this step is to generate a BAM file that contains the alignments of the expressed reads with the reference genome.
            For that, you can use <a href="http://bowtie-bio.sourceforge.net/index.shtml">Bowtie</a> with the following parameters (exact matching).
            Any other read mapper can also do the job.</p>

            <pre class='code'>bowtie -v 0 -f/q --all --best --strata -S &lt;genome&gt; &lt;reads&gt; > output.sam  </pre>

            <p>Reads file must be in FASTA, FASTQ, or colorspace-fasta format. Genome file must be in FASTA format.
            The list of genome assemblies accepted by miRkwood is given in Section <a href="./help.php#reference"> "Select an assembly"</a>
            on the help page.
<!--
            Should you work with any other organism or assembly, then you can upload part of the reference genome (<a href="./help.php#reference">Section "Select an assembly"</a> on the help page).
-->
            </p>


            <h2 id='convert'>Convert the BAM file into a BED file</h2>

            <p>For this step, you should use our custom script mirkwood-bam2bed.pl
            (<a href="../../cgi-bin/mirkwood/web_scripts/getScript.pl?file=mirkwood-bam2bed.pl">download the script</a>).
            mirkwood-bam2bed.pl is a perl script dependent upon the installation of SAMtools.
            In practice, the BED file is up to 10 times smaller than the BAM file, while retaining all information needed to conduct the analysis.
            This allows to reduce significantly the bandwidth necessary to upload the data to miRkwood server. </p>

            <pre class='code'>mirkwood-bam2bed.pl --in /input/file --bed /output/file/ --min X --max Y</pre>

<!--
            <br />
-->
            
            <p>
                --in : path to your input file (format BAM or SAM)<br />
                --bed : path to your output BED file<br />
                --min : keep only reads with length &ge; min (default 18)<br />
                --max  : keep only reads with length &le; max (default 25)
            </p>

            <p>The generated BED file has the following syntax.</p>

<pre class='example'>
1    18092    18112    AAACGTGTAGAGAGAGACTCA    1    -
1    18094    18118    GATTCTTTTGTTTGCCACT    2    +
1    18096    18119    TCGATAGGATCAAGTACATCT    1    +
1    18100    18124    AAGAAGAAAAAGAAGAAGAAGAAG     9    +
</pre>

            <p>In this file, each line is a unique read. The fields are, from left to right: name of the chromosome,
            starting position, ending position, read sequence, number of occurrences of the read in the data, strand.
            Positions follow the BED numbering convention: the first base of the chromosome is considered position 0
            (0-based position) and the feature does not include the stop position. </p>

            <br />

            <p>You are now ready to use miRkwood small RNA-seq
            on your data.
            Go to <a href='/cgi-bin/mirkwood/web_scripts/BAMinterface.pl'>web server</a>.
            </p>

            <br /><br /><br />

        </div><!-- bloc droit-->
        <?php require("/bio1/www/html/lib.inc")?>
        <?php footer("miRkwood","miRkwood", "mirkwood@univ-lille1.fr","2013"); ?>
    </body>

</html>
