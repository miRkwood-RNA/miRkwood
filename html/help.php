<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
  <head>
        <link type='text/css' rel='stylesheet' href='../Style/bioinfo.css' />
		<link type='text/css' rel='stylesheet' href='./style/help.css' />
		<link type='text/css' rel='stylesheet' href='./style/rna.css' />
        <script type='text/javascript' src='/js/miARN.js'></script>
        <title>miRkwood - MicroRNA identification - Help</title>
    </head>
 <body>
        <div class="theme-border"></div>
        <a href="/">
            <div class="logo"></div>
        </a>
        <div class="bloc_droit">
        <?php include("./static/header_menu.txt") ?>

<div class="main-full">


  
<p><a href='/mirkwood/'>miRkwood</a> is a computational pipeline for
  the identification of plant miRNAs and their hairpin precursors.</p>
  <p>
  This page is a user manual for <a href='/cgi-bin/mirkwood/web_scripts/interface.pl'>miRkwood website</a>.
For the full detail of the method implemented in miRkwood, see <a href="method.php">miRkwood method</a>.</p>

<br> 

<hr class='section'/>

<div class="table-of-contents">
  <ol>
    <li><a href="#input_form">Input form</a>
    <ol>
      <li> Enter query sequence </li>
      <li> Parameters</li>
      <li> Submit the job</li> 
    </ol>
   </li>
   <li><a href="#results_page">Result page</a></li>
   <li><a href="#export">Export</a>
   <ol>
     <li> Tabular format (CSV) </li>
     <li> FASTA</li>
     <li> dot-bracket format (plain sequence + secondary structure) </li>
     <li> Full report in document format (ODF)</li>
     <li> GFF format</li>
   </ol>
 </li> 
 <li><a href="#html_report">HTML report</a>
  <ol>
    <li> Header</li>
    <li> Thermodynamic stability</li>
    <li> Conservation of the mature miRNA</li> 
  </ol>
 </li>
</ol>
</div>

<hr class='section'/>


<h2 id='input_form'>Input form</h2>


<p>The input form of miRkwood has three main sections:</p>
<ol>
	<li>Enter query sequence,</li>
	<li>Parameters: select additional annotation criteria for the miRNA precursors,</li>
	<li>Submit the job.</li>
</ol>
<p>Results are displayed on a new page.</p>

<h3>Enter query sequence</h3>

<p>miRkwood input is (Multi) FASTA  format.  Lower-case and upper-case letters are both accepted, as well as T/U.  Other characters, such as N, R, Y,&#8230; are prohibited.</p>

<p>You can either paste or upload a file. The maximum size for a submission is 100 000 nt.</p>

<dl>
<dt>Scan both strands</dt>
<dd>miRkwood normally analyses data in forward direction only. Checking this option will cause the program to search both the forward and reverse complement strands.</dd>

<dt>Mask coding regions</dt>
<dd>This option allows selecting non-coding sequences in the input
  data by masking putative coding regions. It consists in a BlastX
  search against the protein sequences from the chosen species
  (E-value=1E-5).  Currently available species are: <i>Arabidopsis
  thaliana</i> (<a href="http://www.arabidopsis.org/">TAIR</a>, V10), <i>Medicago truncatula</i> (<a href='http://www.jcvi.org/medicago/'>Medicago truncatula genome project</a>, Mt4.0) and <i>Oriza sativa</i> (<a href='http://rice.plantbiology.msu.edu/'>Rice genome annotation project</a>, V7.0).</dd>
</dl>
<h3>Parameters</h3>

<p>miRkwood folds the input sequence to identify miRNA precursor secondary structures (for more explanation on this step, see <a href="method.php">miRkwood method</a>). This gives a set of candidate pre-miRNAs. For each candidate pre-miRNA, it is possible to calculate additional criteria that help to bring further evidence to the quality of the prediction and to distinguish accurate miRNA precursors from pseudo-hairpins.</p>


<p>When all options are unchecked, no criteria is applied,
and the list of predictions is exactly the list of all  candidate
pre-miRNAs obtained with the folding step.</p>

<dl>
<dt id='mfei_definition'>Select only sequences with MFEI &lt; -0.6</dt>
<dd>MFEI is the minimal folding free energy index. It is calculated by the following equation:

<p class='equation'>MFEI = [MFE / sequence length x 100] / (G+C%)</p>

where MFE (minimal free energy) denotes the negative folding free
  energies of a secondary structure, and is calculated using the
  Matthews-Turner nearest neighbor model implemented in <a
  href='http://www.tbi.univie.ac.at/~ronny/RNA/RNAeval.html'>RNAeval</a>. When
  checked, this option removes all candidate pre-miRNAs with an MFEI
  greater than or equal to -0.6. Indeed, more than 96% of plant miRBase precursors have an MFEI smaller than -0.6, whereas pseudo-hairpins show significantly larger values of MFEI (for more details, see <a href='method.php'>miRkwood method</a>).</dd>

<dt id='compute-thermodynamic-stability'>Compute thermodynamic stability</dt>
<dd>The significance of the stability of the sequence can also be measured by comparison with other equivalent sequences. <em><a href='http://www.ncbi.nlm.nih.gov/pubmed/15217813'>Bonnet et al</a></em> have established that the majority of the pre-miRNA sequences exhibit a MFE that is lower than that for shuffled sequences.  We compute the probability that, for a given sequence, the MFE of the secondary structure is different from a distribution of MFE computed with 300 random sequences with the same length and the same dinucleotide frequency. </dd>

<dt id='flag-conserved-mature-mirnas'>Flag conserved mature miRNAs</dt>
<dd>Some families of mature miRNAs are highly conserved through evolution. In this case, it is possible to localize  the mature miRNA within the pre-miRNA  by similarity. For that,  we compare each sequence with  the mature miRNAs of plant (<i>Viridiplantae</i>) deposited in <a href='http://www.mirbase.org/ftp.shtml'>miRBase</a> (Release 20). Alignments are performed with <a href='https://www.ebi.ac.uk/~guy/exonerate/'>Exonerate</a>, which implements an exact model for pairwise alignment. We select alignments with at most three errors (mismatch, deletion or insertion) against the full-length mature miRNA and that occur in one of the two arms of the stem-loop.  The putative location obtained is then validated with <a href='http://www.cs.mcgill.ca/~blanchem/mirdup/'>miRdup</a>, that assesses the stability of the miRNA-miRNA* duplex. Here, it was trained on miRBase Viridiplantae v20.</dd>
</dl>

<h3>Submission</h3>

<p>Each job is automatically assigned an ID.</p>
<dl>
<dt>Job title</dt>
<dd>It is possible to identify the tool result by giving it a name.</dd>
<dt>Email address</dt>
<dd>You can enter your email address to be notified when the job is finished. The email contains a link to access the results for 2 weeks.</dd>
</dl>

<hr class='section'/>

<h2 id='results_page'>Result page</h2>

SCREENSHOT

<p>Results are summarized in a two-way table. Each row corresponds to a pre-miRNA, and each column to a feature. By default, results are sorted by sequence and then by position. It is possible to have them sorted by quality (<a href='#definition_quality'>see definition</a>). You can view all information related to a given prediction by clicking on the row (<a href='#html_report'>see section HTML Report</a>).</p>
<dl>
<dt>Name</dt>
<dd>Name of the original sequence, as specified in the heading of the FASTA format.</dd>

<dt>Position</dt>
<dd>Start and end positions of the putative pre-miRNA in the original sequence.</dd>

<dt>+/- (option)</dt>
<dd>Strand, forward or reverse complement. </dd>

<dt id='definition_quality'>Quality</dt>
<dd>The quality is a distinctive feature of miRkwood. It is a combination of all other criteria described afterwards, and allows to rank the predictions according to the significance, from zero- to three- stars. It is calculated as follows.
<ul>
<li><em>MFEI &lt; -0.8:</em> add one star. This MFEI threshold covers 83% of miRBase pre-miRNAs, whereas it is observed in less than 13% of pseudo hairpins (see <a href="method.php">miRkwood method</a>).</li>

<li><em>Existence of a conserved miRNA in miRBase (alignment):</em> add one star. We allow up to three errors in the alignment with mature miRBase, which corresponds to an estimated P-value of  3E-2 for each pre-miRNA. Alignments with 2 errors or less have an estimated P-value of 4E-3.</li>

<li><em>The location of the mature miRNA obtained by alignment is validated by miRdup:</em> add one star</li>.
</ul>
</dd>

<dt>MFE</dt>
<dd>Value of the minimal free energy (computed with <a href='http://www.tbi.univie.ac.at/~ronny/RNA/RNAeval.html'>RNAeval</a>).</dd>

<dt>MFEI</dt>
<dd>Value of the MFEI (<a href='#mfei_definition'>see definition</a>).</dd>

<dt>Shuffles (option)</dt>
<dd>Proportion of shuffled sequences whose MFE is lower than the MFE of the candidate miRNA precursor (<a href='#compute-thermodynamic-stability'>see Compute thermodynamic stability</a>).  This value ranges between 0 and 1. The smaller it is, the more significant is the MFE.  We report pre-miRNA stem-loops for which the value is smaller than 0.01, which covers more than 89% of miRBase sequences. Otherwise, if the P-value is greater than 0.01, we say that it is non significant, and do not report any value.</dd>

<dt>miRBase alignment (option)</dt>
<dd>This cell is checked when an alignment between the candidate sequence and miRBase is found (<a href='#flag-conserved-mature-mirnas'>see Flag conserved mature miRNAs</a>). It is doubled checked when the location of the candidate mature miRNA is validated by <a href='http://www.cs.mcgill.ca/~blanchem/mirdup/'>miRdup</a>. The alignments are visible in the HTML or ODF report.</dd>

<dt>2D structure</dt>
<dd>You can drag the mouse over the zoom icon to visualize the stem-loop structure of the pre-miRNA. The image is generated with <a href='http://varna.lri.fr/'>Varna</a>.</dd>
</dl>

<hr class='section' />
<h2 id='export'>Export</h2>

<p>Results, or a selection of them, can also be exported to a variety of formats, and saved to a local folder for further analyses.</p>


<h3>Tabular format (CSV)</h3>
<p>
It contains the same information as the result table, plus the FASTA
sequences and the dot-bracket secondary structures. The CSV  
format is supported by spreadsheets like Excel. </p>



<h3>FASTA</h3>

<p>
This is the compilation of all pre-miRNA sequences found in FASTA
format. The header of the FASTA format contains the initial name of the
sequence, as well as the positions and the strand of the predicted pre-miRNA.
</p>

<h3 id='dot_bracket'>Dot-bracket format</h3>
<p>
This is the compilation of all pre-miRNA sequences found, together
with the predicted secondary structure.  The first line contains a
FASTA-like header. The second line contains the nucleic acid
sequence. The last line contains the secondary structure, that is
given as a set of matching brackets.  A base pair between bases
<em>i</em> and <em>j</em> is represented by a "(" at position
<em>i</em> and a ")" at position <em>j</em>. Unpaired bases are
represented by dots (see more explanation on <a
href='http://www.tbi.univie.ac.at/RNA/bracket.html'>Vienna
website</a>).</p>
<pre class='example'>
> Sample_1001-1085, stemloop structure
cugagauacugccauagacgacuagccaucccucuggcucuuagauagccggauacagugauuuugaaagguuugugggguacag
(((...((((.((((((((........(((.((((((((.......)))))))....).)))........)))))))))))))))
</pre>


<h3>Full report in document format (ODF)</h3>

<p>This is an equivalent of the <a href='#html_report'>HTML report</a>, and contains the full report of the predictions. This document format is compatible with Word or OpenOffice.</p>

<h3>GFF format</h3>

<p>General annotation format, that displays the list of positions of pre-miRNA found (see more explanation on <a href='http://www.ensembl.org/info/website/upload/gff.html'>Ensembl documentation</a>)</p>

<pre class="example">
##gff-version 3
# miRNA precursor sequences found by miRkwood have type 'miRNA_primary_transcript'.
# Note, these sequences do not represent the full primary transcript,
# rather a predicted stem-loop portion that includes the precursor.

sample  miRkwood  miRNA_primary_transcript  1001 1085	.  +  . Name=preMir_sample_1001-1085	
</pre>


<p>
Following the convention of miRBase, we consider that hairpin
precursors have type <tt>miRNA_primary_transcript</tt> even if they
are shorter than the primary transcript.
</p> 


<hr class='section' />
<h2 id='html_report'>HTML Report</h2>

<p>The HTML report contains all information related to a given predicted pre-miRNA.</p>

<h3>Header</h3>

<p>The report begins with the following information.</p>

<dl>
<dt>Name</dt>
<dd>Name of the initial sequence, as specified in the heading of the FASTA format.</dd>

<dt>Position</dt>
<dd>Start and end positions of the putative pre-miRNA in the original sequence. The length is indicated in parentheses.</dd>

<dt>Strand</dt>
<dd>+ (forward) or - (reverse complement).</dd>

<dt>GC content</dt>
<dd>Percentage of bases that are either guanine or cytosine.</dd>

<dt>Sequence (FASTA format)</dt>
<dd>Link to download the sequence.</dd>

<dt>Stem-loop structure</dt>
<dd>Link to download the secondary structure in dot-bracket format
  (<a href='#dot_bracket'>see definition</a>). 
</dd>

<dt>Optimal MFE secondary structure</dt>
<dd>If the stem-loop structure is not the MFE structure, we also provide a link to download the MFE structure.</dd>

<dt>Alternative candidates (dot-bracket format)</dt>
<dd>This is the set of stem-loop sequences that overlap the current prediction. The choice between several alternative overlapping candidate pre-miRNAs is made according to the best MFEI.</dd>
</dl>

<p>The stem-loop structure of the miRNA precursor is also displayed with  <a href='http://varna.lri.fr/'>Varna</a>.</p>
<img style='width:400px; display: block; margin: 0 auto;' src='style/structure.png' alt='The stem-loop structure of the miRNA precursor' />

<h3>Thermodynamics stability</h3>
<dl>
<dt>MFE</dt>
<dd>Value of the Minimum Free Energy (computed by <a href='http://www.tbi.univie.ac.at/~ronny/RNA/RNAeval.html'>RNAeval</a>).</dd>

<dt>AMFE</dt>
<dd>Value of the adjusted MFE : MFE &divide; (sequence length) &times; 100</dd>

<dt>MFEI</dt>
<dd>Value of the minimum folding energy index (<a href='#mfei_definition'>see definition</a>).</dd>

<dt>Shuffles</dt>
<dd>Proportion of shuffled sequences whose MFE is lower than the MFE of the candidate miRNA precursor (<a href='#compute-thermodynamic-stability'>see Compute thermodynamic stability</a>).  This value ranges between 0 and 1. The smaller it is, the more significant is the MFE.  We report pre-miRNA stem-loops for which the value is smaller than 0.1, which covers more than 89% of miRBase sequences. Otherwise, if the P-value is greater than 0.1, we say that it is non significant, and do not report any value.</dd>
</dl>

<h3>Conservation of the miRNA</h3>

<p>All alignments with miRBase are reported and gathered according to their positions.</p>
<div class='example'>
<pre class='alignment'>
query            19 ucgcuuggugcaggucggga- 38
                    ||||||||||||| ||||||  
miRBase           1 ucgcuuggugcagaucgggac 21
</pre>
<span class="others">miRBase sequences: <a href='http://mirbase.org/cgi-bin/mirna_entry.pl?acc=MIMAT0001045'>osa-miR168a-5p</a>, <a href='http://mirbase.org/cgi-bin/mirna_entry.pl?acc=MIMAT0001452'>sbi-miR168</a>, <a href='http://mirbase.org/cgi-bin/mirna_entry.pl?acc=MIMAT0001665'>sof-miR168a</a>, <a href='http://mirbase.org/cgi-bin/mirna_entry.pl?acc=MIMAT0001726'>zma-miR168a-5p</a>, <a href='http://mirbase.org/cgi-bin/mirna_entry.pl?acc=MIMAT0001727'>zma-miR168b-5p</a>, <a href='http://mirbase.org/cgi-bin/mirna_entry.pl?acc=MIMAT0018215'>hvu-miR168-5p</a></span>
</div>

<p><tt>query</tt> is the user sequence, and <tt>miRBase</tt>
designates the mature miRNA  (or the miRNA*) found in miRBase. It is possible to access the corresponding miRBAse entry by clicking on the link under the alignment. The report also indicates whether the location is validated by <a href='http://www.cs.mcgill.ca/~blanchem/mirdup/'>miRdup</a>. Finally, we provide an ASCII representation of the putative miRNA within the stem-loop  precursor.</p>

<p>
TO ADD
</p>
  
</div> <!-- main full -->

</div><!-- bloc droit-->
       <?php include("./static/footer.txt") ?>
    </body>
    
</html>
