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
        <title>miRkwood small RNA-seq - Help</title>
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

            <div id="mirkwood_title">User manual</div>

            <br />

            <p>This page is a user manual for the <a href='/cgi-bin/mirkwood/web_scripts/BAMinterface.pl'>miRkwood small RNA-seq web application</a>.</p>

            <p>If this is your first time using miRkwood, we would
            suggest visiting our <a href="./example.php">quick start</a>
            guide before.</p>

            <br /> <br />

            <div>
                <ol>
                    <li><a href="#input_form"> Input form</a>
                        <ol>
                            <li><a href="#reads"> Upload your set of reads</a> </li>
                            <li><a href="#reference"> Select an assembly</a></li>
                            <li><a href="#parameters"> Parameters</a></li>
                                <ol>
                                    <li><a href="#reads_distribution"> Parameters for the processing of the read data</a></li>
                                    <li><a href="#secondary_structure"> Parameters for the secondary structure of the hairpin precursor</a></li>
                                </ol>
                            <li><a href="#submission"> Submit the job</a></li>
                        </ol>
                    </li>
                    <li><a href="#results_page"> Results page</a>
                        <ol>
                            <li><a href="#overview"> Overview</a></li> 
                            <li><a href="#known_mirnas"> Known miRNAs</a></li>
                            <li><a href="#novel_mirnas"> Novel miRNAs</a></li>
                        </ol>
                    </li>
                    <li><a href="#export"> Export</a>
                        <ol>
                            <li><a href="#gff"> GFF format</a></li>
                            <li><a href="#fasta"> FASTA format</a></li>
                            <li><a href="#dot_bracket"> Dot-bracket format</a></li>
                            <li><a href="#csv"> Tabular format (CSV) </a></li>
                            <li><a href="#org"> Full report in ORG format</a></li>
                            <li><a href="#read_cloud">Reads cloud</a></li>
                        </ol>
                    </li>
                    <li><a href="#html_report"> HTML report</a>
                        <ol>
                            <li><a href="#html_known_mirnas"> HTML report for known miRNAs</a></li>
                            <li><a href="#html_novel_mirnas"> HTML report for novel miRNAs</a></li>
                        </ol>
                    </li>
                </ol>
            </div>

            <br />


            <h2 id='input_form'>Input form</h2>

            <h3 id='reads'>Upload your set of reads</h3>

            <p>The input is a set of reads produced by deep sequencing
            of small RNAs and then mapped to a reference genome. For
            that, you should organize your data in a BED file. See the
            <a href="./BED_file.php">detailed instructions</a> on how
            to build this file, or grab a sample file below.</p>

            <ul><li> <a href="sample.bed" download>Download the sample BED file</a> </li></ul> 

            <p>You will find more information on this sample file in our
            <a href="./example.php">quick start</a> guide.</p>


            <h3 id='reference'>Select an assembly</h3>

            <p>miRkwood currently proposes 12 assemblies, which are listed below.</p>

            <ul>
                <li><em>Arabidopsis lyrata</em> - v1.0 - GFF: <a href="http://phytozome.jgi.doe.gov/pz/portal.html">Phytozome</a></li>
                <li><em>Arabidopsis thaliana</em> - TAIR 10 - GFF: <a href="http://phytozome.jgi.doe.gov/pz/portal.html">Phytozome</a></li>
                <li><em>Brassica napus</em> - Bnapus 4.1 - GFF: <a href="http://brassicadb.org/brad/">brassicadb</a></li>
                <li><em>Brassica rapa</em> - Brapa 1.1 - GFF: <a href="http://brassicadb.org/brad/">brassicadb</a></li>
                <li><em>Glycine max</em> - v1.0.23 - GFF: <a href="http://www.ensembl.org/index.html">Ensembl</a></li>
                <li><em>Lotus japonicus</em> - Lj2.5 - GFF: <a href="http://www.kazusa.or.jp/lotus/">kazusa.or.jp</a></li>
                <li><em>Medicago truncatula</em> - JCVI_Mt3.5.2 - GFF: <a href="http://www.jcvi.org/cms/home/">JCVI</a></li>
                <li><em>Oryza sativa</em> - MSU7 - GFF: <a href="http://phytozome.jgi.doe.gov/pz/portal.html">Phytozome</a></li>
                <li><em>Populus trichocarpa</em> - JGI_Poptr2.0 - GFF: <a href="http://phytozome.jgi.doe.gov/pz/portal.html">Phytozome</a></li>
                <li><em>Solanum lycopersicum</em> - SL2.4 - GFF: <a href="http://phytozome.jgi.doe.gov/pz/portal.html">Phytozome</a></li>
                <li><em>Sorghum bicolor</em> - Sorbi1 - GFF: <a href="http://www.plantgdb.org/">PlantGDB</a></li>
                <li><em>Vitis vinifera</em> - Genoscope-20100122 - GFF: <a href="http://phytozome.jgi.doe.gov/pz/portal.html">Phytozome</a></li>
            </ul>

            <p>When available, the assembly is supplemented by two GFF
            files. The first one contains the genome coordinates of
            annotated  CDSs, tRNAs, rRNAs and snoRNAs,  and is
            used to apply masking options described in <a
            href="#parameters">Section 1.3 - Parameters</a>.  The
            source of this file is indicated above, case by case. The
            other GFF file compiles all miRNAs and precursors of miRNAs
            available in MiRBase 21 and is used to detect known miRNAs that are expressed in the sequencing data. </p>


            <h3 id='parameters'>Parameters</h3>

            <p>miRkwood comes with a series of options, that allow to customise the search and enhance the results.
            These options can be divided in two main types: parameters concerning the processing of the sequence reads
            and parameters concerning the secondary structure of the miRNA precursor.
            Note that they only apply to novel miRNAs. For known miRNAs, the user can rely on the additional
            information delivered by miRBase and make up his/her own mind.</p>


            <h4 id='reads_distribution'>Parameters for the processing of the read data</h4>

            <p>The first step of miRkwood is to locate signals into
            the set of mapped reads.  This is performed by scanning
            the set of reads and detecting statistically significant
            clusters of reads. To this end, it is advised to filter out the data
            beforehand with the following options. </p>

            <p id='mask_coding_regions'><b>Mask coding regions:</b>
            This option allows selecting reads that are aligned to non-coding sequences.
            The selection is performed with the GFF annotation file.
            All reads that intersect a CDS feature are removed. Default: checked.</p>

            <p id='filter_tRNA_rRNA'><b>Filter out tRNA/rRNA/snoRNA:</b>
            When this option is checked, products from tRNA, rRNA and snoRNA degradation are filtered out
            from the input reads. This task is performed based on the existing annotation provided
            in the GFF annotation file. All reads that intersect a tRNA, rRNA, or snoRNA feature are removed.
            Default: checked.</p>

            <p id='filter_multimapped'><b>Remove multiply mapped reads:</b>
            All reads that are mapped to more than 5 loci on the reference sequence are discarded.
            This allows to avoid spurious predictions due to transposons. Default: checked.</p>


            <h4 id='secondary_structure'>Parameters for the secondary structure of the hairpin precursor</h4>

            <p>After cluster detection, miRkwood aims at determining
            which sequences can fold into a stemloop structure. This
            gives a set of candidate precursors of miRNAs.
            For each candidate, it is possible to calculate additional criteria
            that help to bring further evidence to the quality of the prediction and
            to distinguish accurate miRNA precursors from pseudo-hairpins. </p>

            <p id='filter_mfei'><b>Select only sequences with MFEI
            &lt; -0.6:</b> MFEI is the minimum folding free energy
            index, and expresses the thermodynamic stability of the precursor.
            It is calculated by the following equation: </p>

            <p class='equation'>MFEI = [MFE / sequence length x 100] / (G+C%)</p>

            <p>where MFE (minimum free energy) denotes the negative folding free energies of a secondary structure,
            and is calculated using the Matthews-Turner nearest neighbor model implemented in RNAeval.
            When checked, this option removes all candidate pre-miRNAs with an MFEI greater than or equal to -0.6.
            Indeed, more than 96% of miRBase precursors have an MFEI smaller than -0.6,
            whereas pseudo-hairpins show significantly larger values of MFEI. Default: checked.</p>

            <p id='thermodynamic_stability'><b>Compute thermodynamic stability:</b>
            The significance of the stability of the sequence can also be measured
            by comparison with other equivalent sequences. <a href="http://www.ncbi.nlm.nih.gov/pubmed/15217813">Bonnet <em>et al</em></a>
            have established that the majority of the pre-miRNA sequences exhibit a MFE that is lower than that for shuffled sequences.
            We compute the probability that, for a given sequence, the MFE of the secondary structure is different from a distribution
            of MFE computed with 300 random sequences with the same length and the same dinucleotide frequency.
            Default: unchecked.</p>

            <p id='flag_conserved_mirnas'><b>Flag conserved mature miRNAs:</b>
            This option permits to check if the predicted miRNA belongs to some known miRNA family.
            For that, we compare the sequence of the precursor with the database of mature miRNAs of plant
            (<em>Viridiplantae</em>) deposited in <a href="http://www.mirbase.org/ftp.shtml">miRBase</a> (Release 21).
            We select alignments with at most three errors (mismatch, deletion or insertion) against
            the full-length mature miRNA and that occur in one of the two arms of the stemloop.
            Moreover, this alignment allows to infer a putative location for the miRNA within the precursor.
            This location is then validated with <a href="http://www.cs.mcgill.ca/~blanchem/mirdup/">miRdup</a>,
            that assesses the stability of the miRNA:miRNA* duplex. Here, it was trained on miRbase <em>Viridiplantae</em> V21.
            Default: checked.</p>

            <p id='filter_bad_hairpins'><b>Filter out low quality hairpins:</b>
            When this option is checked, hairpins found by miRkwood with a global score of 0
            and no alignment with miRbase are discarded. Default: checked.
            </p>

            <h3 id='submission'>Submission</h3>

            <p>Each job is automatically assigned an ID. </p>

            <p><b>Job title.</b> It's possible to identify the tool result by giving it a name. </p>

            <p><b>Email Address.</b> You can enter your email address to be notified when the job is
            finished. The email contains a link to access the results for 2 weeks.</p>


            <h2 id='results_page'>Results page</h2>

            <h3 id="overview">Results overview</h3>
            <p>This page has two main parts. The first one (<i>Options
            summary</i>) is simply a summary of your job parameters. The other one
            (<i>Results summary</i>) provides the detailed results.</p> 


            <p><b>Total number of reads (unique reads):</b>
            This is the total number of reads in your initial file.
            The number of unique reads (obtained after merging identical reads) is indicated in parentheses. </p>

            <p> <b>CoDing Sequences:</b> This is the number of reads that have been discarded by the option <i>Mask coding regions</i>.
            You can list them by clicking on the <i>download</i> link.</p>

            <p><b>rRNA/tRNA/snoRNA:</b> This is the number of reads that have been discarded by the option <i>Filter out tRNA/rRNA/snoRNA</i>.
            You can list them by clicking on the <i>download</i> link.</p>

            <p><b>Multiply mapped reads:</b> This is the number of reads that have been discarded by the option <i> Remove multiply mapped reads</i>.
            You can list them by clicking on the <i>download</i> link.</p>

            <p><b>Orphan cluster of reads:</b> An orphan cluster is a short region in the genome
            that is enriched with aligned reads but that shows no secondary structure compatible
            with a hairpin. You can obtain the list of the such clusters by clicking on the <i>download</i> link (BED file). </p>

            <p><b>Orphan hairpins:</b>
            An orphan hairpin is a candidate with a global score of 0 and showing
            no conservation with miRBase.
            By default, if you select the option "filter out low quality hairpins",
            such hairpins will be discarded automatically and you can obtain the list
            by clicking on the download link (BED file).
            </p>

            <p><b> Unclassified reads:</b> Unclassified reads
            are isolated reads, that do not belong to any
            cluster or any orphan hairpin, or do not fall in any annotated region.
            </p>

            <p><b>Known miRNAs:</b> This is the number of loci annotated as microRNA precursors
            in miRBase that intersect with reads from the BED file.
            You can display detailed results by clicking on the link <i>see results</i>.
            See <a href="#known_mirna">Section 2.2</a>.
            </p>

            <p><b>Novel miRNAs:</b> This is the number of miRNAs found by  miRkwood that have not
            been previously  reported in miRbase.
            You can display detailed results by clicking on the link <i>see results</i>.
            See <a href="#novel_mirna">Section 2.3</a>.</p>


            <h3 id="known_mirnas">Known miRNAs</h3>

            <p>Known miRNAs are miRNAs that are already present with their precursor in the
            miRBase database (version 21). We consider that a known
            microRNA is found in the data as soon as there is at least
            one read  on the precursor sequence. The quality score (see definition below)
            helps to determine which are the best candidates.</p>

            <br />
            <img style='width:660px; display: block; margin: 0 auto;'src='../style/results_knownmirna.png' alt='results table' />
            <br />

            <p>The list of all known miRNAs found is displayed in a two-way table.
            Each row corresponds to a pre-miRNA, and each column to a feature.
            By default, results are sorted by sequence and then by position.
            It is possible to have them sorted by quality (see definition below).
            You can view all information related to a given prediction by clicking on the row
            (see <a href="#html_report">section HTML Report</a>).</p>

            <p><b>Chr:</b> Number of the chromosome.</p>

            <p><b>Position:</b> Start and end positions of the miRNA precursor, as documented in miRBase.</p>

            <p><b>+/- :</b> Strand, forward (+) or reverse (-).</p>

            <p id='known_mirna_quality'><b>Quality:</b> This score measures the consistency between the distribution of reads
            along the locus and the annotation provided in miRbase.
            It ranges between 0 and 2 stars, and  is calculated as follows. </p>

            <ul>
                <li>the locus contains more than 10 reads: add one star.</li>
                <li>more than half of the reads intersect either with the miRNA or the miRNA*: add one star. </li>
            </ul>

            <p><b>miRBAse name:</b> miRBase identifier.</p>

            <p><b>Reads:</b> Number of reads included in the locus.</p>

            <p><b>miRNA Sequence:</b> Sequence of the miRNA.</p>

            <p><b>miRNA Length:</b> Length of the miRNA.


            <h3 id="novel_mirnas">Novel miRNAs</h3> 

            <p>Novel miRNAs are miRNAS that are not reported in miRBase. The prediction is supported
            by the presence of a stemloop secondary structure, a
            significant read coverage and read distribution.</p>

            <br />
            <img style='width:660px; display: block; margin: 0 auto;'src='../style/results_novelmirna.png' alt='results table' />
            <br />

            <p>Each row corresponds to miRNA precursor, and each
            column to a feature. By default, results are sorted by
            sequence and then by position.</p>

            <p>It is possible to have the results sorted by quality.
            The <b id='novel_mirna_quality'>quality</b> is the sum of four
            values: </p>
            <ul>
                <li>the existence of a miRNA sequence, </li>
                <li>the validation of the miRNA/miRNA* duplex by <a href="http://www.cs.mcgill.ca/~blanchem/mirdup/">miRdup</a>, </li>
                <li>the score of reads distribution, </li>
                <li>the value of the MFEI of the precursor secondary structure (&lt;-0.8)
                    (see definitions below). </li>
            </ul>

            <p>You can view all information related
            to a given prediction by clicking on the row (see <a
            href="#html_report">section HTML Report</a>).</p>

            <p><b>Chr:</b> Number of the chromosome.</p>

            <p><b>Position:</b> Start and end positions of the
            putative miRNA precursor in the original sequence in 1 based notation (consistently to the GFF format).</p>

            <p><b>+/-:</b> Strand, forward (+) or reverse (-). </p>

            <p><b>Reads:</b> Total number of reads included in the locus.</p>

            <p><b>Reads distribution:</b> This score, ranging from 0
            to 3-stars, allows to qualify the pattern of reads mapping
            to a putative microRNA precursor. It aims at determining
            if this distribution of reads presents a typical 2-peaks
            profile, corresponding to the guide miRNA and the miRNA* respectively. </p>

            <ul>
                <li><em>Number of reads:</em> The  locus has either
                at least 10 reads mapping to each arm, or at least 100
                reads mapping in total. </li>
                <li><em>Precision of the precursor processing :</em>
                At least 75% of reads start in a window [-3,+3]
                centered around the start position of the miRNA, or
                [-5,+5] on the opposite arm of the stemloop.  </li>
                <li><em>Presence of the miRNA:miRNA* duplex:</em>
                There is at least one read in the window [-5,+5] on
                the strand of the miRNA*. </li>

            </ul>

            <p>Each criterion contributes equally to the overall ranking, and adds one star.</p>

            <p><b>MFEI:</b> This is the minimum folding free energy
            index. This value expresses the thermodynamic stability of
            the precursor and is calculated by the following equation: </p>

            <p class='equation'>MFEI = [MFE / sequence length x 100] / (G+C%)</p>

            <p>where MFE is the  minimum free energy of the secondary structure 
            (computed with <a href="http://www.tbi.univie.ac.at/RNA/RNAeval.html">RNAeval</a>)
            When the MFEI is  &lt; -0.8, then the value is displayed
            in purple, indicating a significantly stable hairpin. This
            MFEI threshold covers 83% of miRBase miRNA precursors,
            whereas it is observed in less than 13% of pseudo hairpins.</p>

            <p><b>Shuffles (option):</b> proportion of shuffled sequences whose MFE is lower
            than the MFE of the candidate miRNA precursor (see <a href="#thermodynamic_stability">Compute thermodynamic stability</a>).
            This value ranges between 0 and 1. The smaller it is, the more significant is the MFE.
            We report pre-miRNA stemloops for which the value is smaller than 0.01,
            which covers more than 89% of miRBase sequences. Otherwise, if the P-value is greater than 0.01,
            we say that it is non significant, and do not report any value.</p>

            <p><b>miRNA Sequence:</b> Sequence of the miRNA. It is the sequence
            of the most common read, with a frequency of at least 33%.</p>

            <p><b>miRNA Length:</b> Length of the miRNA (when existing).</p>

            <p><b>miRNA Weight:</b> Depth of the miRNA read (when existing)
            divided by the number of possible alignments of this read in the genome.</p>

            <p><b>Conserved miRNA (option):</b> This cell is checked <img src='../style/check.png' alt='arobas' style='width:15px; height:15px;' />
            when an alignment between the miRNA precursor sequence and plant mature miRNA database of miRBase is found
            (see <a href="#flag_conserved_mirnas">Flag conserved mature miRNAs</a>).
            It is doubled checked <img src='../style/check.png' alt='arobas' style='width:15px; height:15px;' />
            <img src='../style/check.png' alt='arobas' style='width:15px; height:15px;' />
            if at least 40% of reads overlap with the miRBase sequence.
            The alignments are visible in the HTML report.</p>


            <h2 id='export'>Export </h2>

            <p>Results, or a selection of them, can be exported to a variety of formats, and saved to a local folder for further analyses. </p>

            <p id='gff'><b>GFF:</b> General annotation format, that displays the list of positions of pre-miRNA found
            (see more explanation on <a href="http://www.ensembl.org/info/website/upload/gff.html">Ensembl documentation</a>).</p>

            <p id='fasta'><b>FASTA:</b> This is the compilation of all pre-miRNA sequences found. </p>

            <p id='dot_bracket'><b>Dot-bracket notation:</b> This is the compilation of all pre-miRNA sequences found,
            together with the predicted secondary structure. The secondary structure is given as a set of matching parentheses
            (see more explanation on <a href="https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/rna_structure_notations.html">Vienna website</a>).  </p>

            <p id='csv'><b>CSV <em>(comma separated value)</em>:</b>
            It contains the same information as the result table, plus the FASTA sequences
            and the dot-bracket secondary structures. This tabular format is supported by spreadsheets like Excel. </p>

            <p id='org'><b>ORG:</b>
            This is an equivalent of the <a href="#html_export">HTML report</a>, and contains the full report of the predictions.
            This file can be easily edited by the user. </p>

            <p id='read_cloud'><b> Reads cloud:</b>  This archive is a
            compilation of all reads clouds. Each reads cloud is a text file that
            summarizes all information available for a potential precursor:
            positions, sequence, secondary structure, existence of an alignment
            with miRbase, distribution of mapped reads. It can easily be parsed.</p>

            <h2 id='html_export'>HTML report </h2>

            <p>The HTML report contains all information related to a given predicted pre-miRNA.</p>

            <h3 id='html_known_mirnas'>HTML report for known miRNAs</h3>

            <ul>
                <li><b>miRBase name:</b> miRBase identifier, and link to access the miRBase entry.</li>
                <li><b>Chromosome:</b> Name of the chromosome</li>
                <li><b>Position:</b> Start and end positions of the pre-miRNA, such as indicated in miRBase
                    in 1-based notation (consistently with the GFF format). The length is indicated in parentheses.
                </li>
                <li><b>Strand:</b> + (forward) or - (reverse).</li>
                <li><b>G+C content:</b> Percentage of bases that are either guanine or cytosine.</li>
                <li><b>miRNA sequence:</b> Sequence of the predicted miRNA.</li>
                <li><b>Sequence (FASTA format):</b> Link to download the sequence of the precursor.</li>
                <li><b>Stemloop structure:</b> Link to download the secondary structure in <em>dot-bracket format</em>.
                    The first line contains a FASTA-like header. The second line contains the nucleic acid sequence.
                    The last line contains the set of associated pairings encoded by brackets and dots.
                    A base pair between bases <em>i</em> and <em>j</em> is represented by a '(' at position <em>i</em> and a ')'
                    at position <em>j</em>. Unpaired bases are represented by dots.
                </li>
                <li><b>Alternative candidates:</b>
                    This is the set of stemloop sequences that overlap the current prediction.
                    The choice between several alternative overlapping candidate pre-miRNAs
                    is made according to the best MFEI.
                </li>
            </ul>
<pre class='example'>
>  1:234009-234092,-, stemloop structure
GAAAUGAUGCGCAAAUGCGGAUAUCAAUGUAAAUCAGGGAGAAGGCAUGAUAUACCUUUAUAUCCGCAUUUGCGCAUCAUCUCU
((.(((((((((((((((((((((.((.(((.((((.(.......).)))).))).)).))))))))))))))))))))).)).
</pre>

            <p><b> Quality:</b> It details the computation of the quality score as explained in
                <a href="#known_mirnas">Section 2.2 Known miRNAs</a>.</p>

            <p><b> Reads</b></p>

            <ul>
                <li><b>Number of reads:</b> The total number of reads included in the locus.</li>
                <li><b>Reads length distribution:</b> This plot gives the sequence length distribution
                for all reads in the precursor locus.</li>
                <li><b>Reads cloud:</b> This is a visual representation of all reads included in the locus.
                <br />
                Each <tt>**********</tt> string is a unique read. Its length and its depth 
                (its number of occurrences in the set of reads) are reported at the end of the dotted line.
                <tt><------miRBase------> </tt>indicates the positions of the miRNA referenced in miRBase.
                </li>
            </ul>

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

            <p><b>Thermodynamics stability</b></p>

            <ul>
                <li><b>MFE:</b> Value of the Minimum Free Energy (computed by <a href="http://www.tbi.univie.ac.at/RNA/RNAeval.html">RNAeval</a>) </li>
                <li><b>AMFE:</b> Value of the adjusted MFE : MFE/(sequence length) x 100</li>
                <li><b>MFEI:</b> Value of the minimum folding energy index (as defined <a href="#filter_mfei">here</a>)</li>
                <li><b>Shuffles (option):</b> Proportion of shuffled sequences whose MFE is lower than the MFE of the candidate miRNA precursor
                    (see <a href="#thermodynamic_stability">Compute thermodynamic stability</a>).
                    This value ranges between 0 and 1. The smaller it is, the more significant is the MFE.
                    We report pre-miRNA stemloops for which the value is smaller than 0.01, which covers more than 89% of miRBase sequences.
                    Otherwise, if the P-value is greater than 0.01, we say that it is non significant, and do not report any value.</li>
            </ul>



            <h3 id='html_novel_mirnas'>HTML report for novel miRNAs</h3>

            <ul>
                <li><b>Chromosome:</b> Name of the chromosome</li>
                <li><b>Position:</b> Start and end positions of the pre-miRNA, such as indicated in miRBase
                    in 1-based notation (consistently with the GFF format). The length is indicated in parentheses.
                </li>
                <li><b>Strand:</b> + (forward) or - (reverse).</li>
                <li><b>G+C content:</b> Percentage of bases that are either guanine or cytosine.</li>
                <li><b>miRNA sequence:</b> Sequence of the predicted miRNA.</li>
                <li><b>miRNA depth:</b> Depth of the miRNA read (weight: miRNA depth
                    divided by the number of possible alignments of this read in the genome.)</li>
                <li><b>Candidates with the same miRNA:</b> List of other miRkwood predictions
                    that involve the same miRNA sequence.</li>
                <li><b>Sequence (FASTA format):</b> Link to download the sequence of the precursor.</li>
                <li><b>Stemloop structure:</b> Link to download the secondary structure in <em>dot-bracket format</em>.
                    The first line contains a FASTA-like header. The second line contains the nucleic acid sequence.
                    The last line contains the set of associated pairings encoded by brackets and dots.
                    A base pair between bases <em>i</em> and <em>j</em> is represented by a '(' at position <em>i</em> and a ')'
                    at position <em>j</em>. Unpaired bases are represented by dots.
                </li>
            </ul>

<pre class='example'>
>  1:234009-234092,-, stemloop structure
GAAAUGAUGCGCAAAUGCGGAUAUCAAUGUAAAUCAGGGAGAAGGCAUGAUAUACCUUUAUAUCCGCAUUUGCGCAUCAUCUCU
((.(((((((((((((((((((((.((.(((.((((.(.......).)))).))).)).))))))))))))))))))))).)).
</pre>

            <ul> 
                <li><b>Optimal MFE secondary structure:</b> If the stemloop structure is not the MFE structure,
                    we also provide a link to download the MFE structure.
                </li>
                <li><b>Alternative candidates (dot-bracket format):</b>
                    This is the set of stemloop sequences that overlap the current prediction.
                    The choice between several alternative overlapping candidate pre-miRNAs is made according to the best MFEI.
                </li>
            </ul>

            <p>The stemloop structure of the miRNA precursor is also displayed with <a href="http://varna.lri.fr/">Varna.</a></p>

            <img style='width:400px; display: block; margin: 0 auto;' src='../style/varna.png' alt='Varna image' />


            <p><b>Quality:</b></p>
            This score includes the reads ditribution score (as defined <a href="#novel_mirnas">here</a>)
            and extends it with some additional criteria.
            It ranges from 0 to 6.

            <ul>
                <li><b>MFEI &lt; -0.8:</b> This MFEI threshold covers 83% of miRBase pre-miRNAs,
                    whereas it is observed in less than 13% of pseudo hairpins.</li>
                <li><b>Criteria number of reads:</b> The locus has either
                    at least 10 reads mapping to each arm, or at least 100
                    reads mapping in total. </li>
                <li><b>Existence of a miRNA:</b> It is the sequence
                    of the most common read, if the frequency is at least 33%. </li>
                <li><b>The miRNA is validated by miRdup:</b>  This measures the stability of the miRNA/miRNA* duplex
                    (<a href="http://www.cs.mcgill.ca/~blanchem/mirdup/">see more information</a>), </li>
                <li><b>Criteria presence of the miRNA:miRNA* duplex:</b> There is at least one read in the window [-5,+5] on
                the strand of the miRNA*. </li>
                <li><b>Criteria precision of the precursor processing:</b> At least 75% of reads start in a window [-3,+3]
                    centered around the start position of the miRNA, or
                    [-5,+5] on the opposite arm of the stemloop. </li>
            </ul>


            <p><b> Reads</b></p>

            <ul>
                <li><b>Number of reads:</b> The total number of reads included in the locus.</li>
                <li><b>Reads length distribution:</b> This plot gives the sequence length distribution
                for all reads in the precursor locus.</li>
                <li><b>Reads cloud:</b> This is a visual representation of all reads included in the locus.
                    <br />
                    Each <tt>**********</tt> string is a unique read. Its length and its depth
                    (the number of occurrences of this read  in the total set of reads) are reported at the end of the dotted line.
                    Each <tt><------miRBase------></tt> string indicates the existence of an alignmnent with a sequence from miRBase.
                    More information on this alignment is provided in the sequel of the report, in Section <i>Conservation of the mature miRNA</i>.
                </li>
            </ul>

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

            <p><b>Thermodynamics stability</b></p>

            <ul>
                <li><b>MFE:</b> Value of the Minimum Free Energy (computed by <a href="http://www.tbi.univie.ac.at/RNA/RNAeval.html">RNAeval</a>) </li>
                <li><b>AMFE:</b> Value of the adjusted MFE : MFE/(sequence length) x 100</li>
                <li><b>MFEI:</b> Value of the minimum folding energy index (as defined <a href="#filter_mfei">here</a>)</li>
                <li><b>Shuffles (option):</b> Proportion of shuffled sequences whose MFE is lower than the MFE of the candidate miRNA precursor
                    (see <a href="#thermodynamic_stability">Compute thermodynamic stability</a>).
                    This value ranges between 0 and 1. The smaller it is, the more significant is the MFE.
                    We report pre-miRNA stemloops for which the value is smaller than 0.01, which covers more than 89% of miRBase sequences.
                    Otherwise, if the P-value is greater than 0.01, we say that it is non significant, and do not report any value.</li>
            </ul>

            <p><b>Conserved mature miRNA</b></p>

            <p>All alignments with miRBase are reported and gathered according to their positions. </p>

            <img style='width:610px; display: block; margin: 0 auto;' src='../style/alignment.png'' alt='alignment' />

            <p><i>query</i> is the user sequence, and <i>miRBase</i> designates the mature miRNA found in miRBase.
            It is possible to access the corresponding miRBase entry by clicking on the link under the alignment.
            Finally, we provide an ASCII representation of the putative miRNA within the stemloop  precursor.</p>

            <img style='width:600px; display: block; margin: 0 auto;' src='../style/hairpin_with_mature.png'' alt='hairpin with mature' />

        </div>
        <?php require("/bio1/www/html/lib.inc")?>
        <?php footer("miRkwood","miRkwood", "mirkwood@univ-lille1.fr","2018"); ?>
    </body>

</html>
