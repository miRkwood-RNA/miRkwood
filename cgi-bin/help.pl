#!/usr/bin/perl -w




print <<DATA
Content-type: text/html


<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<meta name="keywords" content="RNA, ARN, mfold, fold, structure, prediction, secondary structure" />
<link title="test" type="text/css" rel="stylesheet" href="/arn/css/script.css" />

		<script src="/arn/js/miARN.js" type="text/javascript" LANGUAGE="JavaScript"></script>
	
		
		<title>MicroRNA identification</title>


</head>
<body>
<div class="theme-border"></div>
<div class="logo"></div>
<!-- debut du fichier qui represente le menu colore vertical de gauche -->

<div class="bloc_gauche">
  <div class="nav-top"></div>
  <div class="menu">
    
    <!-- SEQUOIA -->
    <ul class="menu_partie">
      <li class ="case_titre Research_titre"> 
	<a href="/index.php">bioinfo.lifl.fr</a>
      </li>
      <li class ="case Research"> 
	<a href="http://www.lifl.fr/bonsai">Bonsai</a>
      </li>
    </ul>
    
    <!-- Mreps -->
    <ul class="menu_partie"> 
      <li class ="case_titre Mreps_titre">
	<a href="/mreps/index.php">mreps</a>
      </li>
    </ul>
    
    <!-- Yass -->
    <ul class="menu_partie">
      <li class ="case_titre Yass_titre">
	<a href="/yass/index.php">YASS</a>
      </li>
    </ul> 

    <!-- Magnolia -->
    <ul class="menu_partie">
      <li class="case_titre Magnolia_titre">
	<a href="/magnolia/index.php">Magnolia</a>
      </li>
    </ul>
    

    <ul class="menu_partie">
      <li class="case_titre proteins_titre">
	<a href="/proteins/index.php">Proteins</a>
      </li>
    <!-- Path -->
      <li class ="case proteins">
        <a href="/path/index.php">Path</a>
      </li>
    <!-- Protea -->
      <li class="case proteins">
	<a href="/protea/index.php">Protea</a>
      </li>    
    <!-- Reblosum -->
      <li class="case proteins">
	<a href="/reblosum/index.php">ReBLOSUM</a>
      </li>
    </ul>
        
    <!-- RNA -->
    <ul class="menu_partie"> 
      <li class="case_titre Rna_titre">
	<a href="/RNA/index.php">RNA</a>
      </li>
      <li class="case Rna">
	<a href="/RNA/carnac/index.php">Carnac</a>
      </li>
      <li class="case Rna">
	<a href="/RNA/RNAfamily/rnafamily.php">RNAfamily</a>
      </li>
      <li class="case Rna">
	<a href="/RNA/gardenia/index.php">Gardenia</a>
      </li>
      <li class="case Rna">
	<a href="/RNA/regliss/index.php">Regliss</a>
      </li>
      <li class="case Rna">
	<a href="/CGseq/index.php">CG-seq</a>
      </li>

    </ul>
    
    <!-- TFM -->
    <ul class="menu_partie"> 
      <li class="case_titre Tfm_titre">
	<a href="/TFM/">TFM</a>
      </li>
      <li class ="case Tfm">
	<a href="/TFM/TFME/">TFM-Explorer</a>
      </li>
      <li class ="case Tfm">
	<a href="/TFM/TFMscan/">TFM-Scan</a>
      </li>
      <li class="case Tfm">
	<a href="/TFM/TFMpvalue/">TFM-Pvalue</a>
      </li>      
      <li class="case Tfm">
	<a href="/TFM/TFM-CUDA/">TFM-CUDA</a>
      </li>      
    </ul>

    <!-- NRPS -->
    <ul class="menu_partie no-margin">
      <li class ="case_titre Nrps_titre">
	<a href="/norine/">NRPS</a>
      </li>
      <li class ="case Nrps"> 
	<a href="/norine/form.jsp">Norine</a>
      </li>

    </ul>
    
  </div> <!-- left menu -->
  <div class="nav-bottom"></div>  
</div> <!-- bloc gauche-->


<script type="text/javascript">
var pkBaseURL = (("https:" == document.location.protocol) ? "https://bioinfo.lifl.fr/piwik/" : "http://bioinfo.lifl.fr/piwik/");
document.write(unescape("%3Cscript src='" + pkBaseURL + "piwik.js' type='text/javascript'%3E%3C/script%3E"));
</script>
<script type="text/javascript">
try {
var piwikTracker = Piwik.getTracker(pkBaseURL + "piwik.php", 1);
piwikTracker.trackPageView();
piwikTracker.enableLinkTracking();
} catch( err ) {}
</script>
<noscript><p><img src="http://bioinfo.lifl.fr/piwik/piwik.php?idsite=1" style="border:0" alt=""/></p></noscript>


<!-- End Piwik Tag -->

<!-- fin du fichier menu.txt -->


<div class="bloc_droit">
<div class="frametitle">
 <h1>miREST :: identification of miRNA/miRNA hairpins in plants</h1>
   <div class="menu_central">
     <a href="./interface.pl">home</a>
 	<a href="./help.pl">help</a>
     <a href="./id.pl">retrieve result with an ID</a>
     </div>
</div>

<div class="main">

<h2>What is miREST</h2>

<p>

<a name="fasta"></a> 
miREST is a pipeline for the identification of miRNA and their hairpins precursors from assembled EST (Sanger or 454) or from short (<10 kb) genomic sequences. It is dedicated to plants. Fonctionnement du programme selon choix final et paramétrage des outils.
</p>

<h2>Input</h2>
<p>
miREST input is (Multi)Fasta format. A sequence in FASTA format consists of a single-line description, followed by lines of sequence data. The first character of the description line is a greater-than (">") symbol in the first column. All lines should be shorter than 80 characters:
Fasta file example 
</p>
<div class="exemple">
<pre>&gt; Name of the sequence 1
ctgcgagcgcgcgatgatagcgcggcgagcatgtagcatgctagctgtcgcgagcact
cggccgagatcaggcgatgcatgcgcagggagcagcgagcgacgagcacagcatgcta
gctagatgcatgctgtaggcagc
cgccgagagacgatggagctgc

</pre>
</div>

<p>
You can either paste or upload a file. Lower-case and upper-case letters are both accepted. Note that even if the full standard IUPAC nucleic acid code, numerical digits 0, ..., 9, - and dot . symbols are accepted, they are ignored by miREST. 
</p>  


<h2>Mask coding regions</h2>
<p>
This option allows selecting non-coding sequences in the input data by masking putative coding regions. It consists in a BlastX search against the translated coding sequences from the chosen species [A VERIFIER: 'translated coding sequences' ou 'proteins']. Available species are currently: A. thaliana, XXX It may be slow according the number and the length of sequences in input.
</p>


<h2>Additional features</h2>
<p>

For each candidate sequence, it is possible to calculate additional features that help to bring furhter evidence to the quality of the prediction and to distinguish accurate miRNA precursors from random stemloops. These features rely mainly on statistical features based on thermodynamics and conservation of the mature miRNA sequence.
</p>

<p>
<b>MFE/AMFE/AMFEI</b><br>
<li><b> MFE </b>(minimal free energy)<br><p> denotes the negative folding free energies (ΔG) calculated using RNAfold. MFE is estimated by considering the minimum energy values obtained by complementary base pairs decreased by the stacking energy of successive base pairs or increased by the destabilizing energy associated with non- complementary bases (Zuker and Stiegler, 1981; Schuster et al., 1994) [from Bonnet 2004, phrase à modifier]. It may be slow with a large number of sequences.</p><br>
<br>
<li><b> AMFE</b><br><p> is an adjusted MFE. It is calculated using the following equation: AMFE = MFE / sequence length x 100.  Plant pre-miRNAs are known to ….</p>
<br><br>
<li><b> MFEI</b><br><p> is the minimal folding free energy index. It is calculated using the following equation: MFEI = [MFE / sequence length x 100] / GC%. Plant pre-miRNAs are known to ….</p>
<br><br>
<b>Thermodynamic stability probability (Randfold): </b><br><p>Randfold compute the probability that, for a given sequence, the MFE of the secondary structure is different from a distribution of MFE computed with 999 random sequences with the same dinucleotide frequency (Bonnet et al. 2004). <br>
<b>Warning</b>: It may be slow with a large number of candidate sequences.</p>
<br><br>
<b>Compute self-containment index:</b>
<br>
<p>for each sequence of interest, a set of 1000 random sequences of the same length are generated... A compléter
<b>Warning</b>: It may be slow with a large number of candidate sequences.
</p><br>
<b>Align against mature microRNAs from miRBase</b>
<p>
It allows the identification of conserved miRNAs by comparing each sequence against the mature miRNAs of plant (Viridiplantae) deposited in miRBase (Realease 18). Alignments are performed with Exonerate (Slater et al. 2005). We select all  alignments with at most three errors (mismatch, deletion or insertion) against the full-length mature miRNA. 
</p>



<h2>Ouput</h2>
<br>
<p>
Results are available in a two-way table. Each row corresponds to a single prediction, and each column to a given property of the prediction. 
By default, the following information is displayed:
<li><b>name</b>:<p> the name of the initial sequence, as specified in the headiong of the Fasta format</p>
<li><b>position</b>:<p> start and end position of the putative pre-miRNA on the initial sequence</p>
<li>XXX:<p> sequence and predicted secondary structure of the putative pre-miRNA (format ? vienna, ct ?)</p>
<li><b>2D drawing</b>:<p> planar representation of the sequence and predicted secondary structure of the putative pre-miRNA. The figure is generated with XXX<p>
<p>
Then there is one additional column for each selected additional feature.</p>
<li><b>MFEI</b>:<p> value of the MFEI</p>
<li><b>miRBase alignment</b>:   
<li><b>randfold</b>:<p> P-value returned by ranfold. This P-value is the proportion of equivalent randomized sequences with the same dinucleotide frequency whose MFE is lower than the MFE of the candidate miRNA precursor. Its value ranges between 0 and 1. O indicates that all random sequences have a MFE larger than the one observed with the </p>
<li>


It is possible to have more information on each value by clicking on the corresponding entry. It is also possible to view all information related to a given prediction by clicking on the first entry of the row (column name). Finally, it is possible to view all possible values of a given column by clicking on the heading of the column.
</p> 

 
</div>

 <div class="copyright">
          2004 - <a href="http://www.lifl.fr/bonsai">Bonsai bioinformatics</a>
         </div>
</div><!-- bloc droit-->

</body>
</html>

DATA
###End###
