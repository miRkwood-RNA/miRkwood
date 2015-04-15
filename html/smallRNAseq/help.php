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
            
            <p>This page is a user manual for <a href='/cgi-bin/mirkwood/web_scripts/BAMinterface.pl'>miRkwood small RNA-seq website</a>. This web application allows to search for microRNAs in small RNA sequencing data in plants.</p>
            
            <p>The method implemented in miRkwood is described in full detail in this other page.</p>

            <br /> 
            
            <div class="table-of-contents">
                <ol>
                    <li><a href="#input_form"> Input form</a>
                        <ol>
                            <li><a href="#reads"> Enter the set of expressed reads</a> </li>
                            <li><a href="#reference"> Enter the reference sequence</a></li>
                            <li><a href="#parameters"> Parameters</a></li>
                                <ol>
                                    <li><a href="#reads_distribution"> Parameters concerning the distribution of reads</a></li>
                                    <li><a href="#secondary_structure"> Parameters concerning the secondary structure of the hairpin precursor</a></li>
                                </ol>
                            <li><a href="#submission"> Submit the job</a></li> 
                        </ol>
                    </li>
                    <li><a href="#results_page"> Result page</a>
                        <ol>
                            <li><a href="#known_mirnas"> Known miRNAs</a></li>
                            <li><a href="#novel_mirnas"> Novel miRNAs</a></li>
                        </ol>
                    </li>    
                    <li><a href="#export"> Export</a>
                        <ol>
                            <li><a href="#gff"> GFF format</a></li>
                            <li><a href="#fasta"> FASTA</a></li>
                            <li><a href="#dot_bracket"> dot-bracket format<br/>(plain sequence + secondary structure) </a></li>
                            <li><a href="#csv"> Tabular format (CSV) </a></li>
                            <li><a href="#odf"> Full report in document format (ODF)</a></li>
                            
                        </ol>
                    </li> 
                    <li><a href="#html_report"> HTML report</a>
                        <ol>
                            <li><a href="#html_known_mirnas"> HTML Report for known miRNAs</a></li>
                            <li><a href="#html_novel_mirnas"> HTML Report for novel miRNAs</a></li>
                        </ol>
                    </li>
                </ol>
            </div>
            
            
            <h2 id='input_form'>Input form</h2>
            
            <h3 id='reads'>Enter the set of expressed sequence tags</h3>
            
            <p>The input is a set of reads produced by deep sequencing of small RNAs and then mapped to a reference genome. For that, you should organize your data in a BED file. <a href="./prepare_file.php">See the detailed instructions.</a> </p>
            
            
            <h3 id='reference'>Select an assembly</h3>
            
            <p>miRkwood proposes 7 assemblies, which are listed below.</p>
            
            <ul>
                <li><em>Arabidopsis lyrata</em> - v1.0</li>
                <li><em>Arabidopsis thaliana</em> - TAIR 10</li> 
                <li><em>Brassica rapa</em> - Brapa 1.1</li>
                <li><em>Glycine max</em> - v1.0.23</li>
                <li><em>Lotus japonicus</em> - Lj2.5</li>
                <li><em>Oryza sativa</em> - MSU7</li>
                <li><em>Solanum lycopersicum</em> - SL2.4</li>
            </ul>
            
            <p>Each assembly is supplemented by two GFF format files. The first one contains the genome coordinates of  annotated mRNAs, CDS, tRNAs and rRNAs,  and is used to apply masking options described in <a href="#parameters">Section 1.3 - Parameters</a>.  The source of this file is indicated above, case by case. The other GFF file compiles all miRNAs and pre-miRNAs available in MiRBase 21 and is used to detect known miRNAs that are expressed in the sequencing data. </p>
            
            
            <h3 id='parameters'>Parameters</h3>
            
            <p>miRkwood comes with a series of options, that allow to customize the search and enhance the results. These options can be divided in  two main types: parameters concerning the processing of the sequence reads and parameters concerning the secondary structure of the miRNA precursor. Note that they only apply  to novel miRNAs. For known miRNAs, the user can rely on the additional information delivered by miRBase and make up his/her own mind.</p>
            
            
            <h4 id='reads_distribution'>Parameters concerning the distribution of reads</h4>
            
            <p>The first step of miRkwood is to locate signals into the set of mapped reads.  This is performed by scanning the set of reads and detecting peaks (see <a href="../method.php">miRkwood method</a> for more explanation on this step).</p>
            
            <p><b>Remove multiply mapped reads:</b> All reads that are mapped to more than 5 loci on the reference sequence are discarded.  This allows to XXXX transposons. Default: checked.</p>
            
            <p><b>Filter out tRNA/rRNA:</b> When this option is checked, products from tRNA and rRNA degradation are filtered out from the input reads.  This task is performed based on the existing annotation provided in the GFF annotation file. All reads that intersect a tRNA or rRNA feature are removed. Default: checked.</p>
            
            <p><b>Mask coding regions:</b> This option allows selecting reads that are aligned to non-coding sequences. The selection is performed with the GFF annotation file. All reads that intersect  a CDS feature are removed. Default: checked.</p>
            
            
            <h4 id='secondary_structure'>Parameters concerning the secondary structure of the hairpin precursor</h4>
            
            <p>After peak detection,  miRkwood aims at determining which sequences can fold into a MIR stem-loop structure (see <a href="../method.php">miRkwood method</a> for more explanation on this step). This gives a set of candidate pre-miRNAs. For each candidate pre-miRNA, it is possible to calculate additional criteria that help to bring further evidence to the quality of the prediction and to distinguish accurate miRNA precursors from pseudo-hairpins. </p>
            
            <p><b>Select only sequences with MFEI &lt; -0.6:</b> MFEI is the minimal folding free energy index. It is calculated by the following equation: </p>
            
            <p class='equation'>MFEI = [MFE / sequence length x 100] / (G+C%)</p>
            
            <p>where MFE (minimal free energy) denotes the negative folding free energies of a secondary structure, and is calculated using the Matthews-Turner nearest neighbor model implemented in RNAeval. When checked, this option removes all candidate pre-miRNAs with an MFEI greater than or equal to -0.6.  Indeed, more than 96% of miRBase precursors have an MFEI smaller than -0.6, whereas pseudo-hairpins show significantly larger values of MFEI (see <a href="../method.php">miRkwood method</a> for more details). Default: checked.</p>
            
            <p id='thermodynamic_stability'><b>Compute thermodynamic stability:</b> The significance of the stability of the sequence can also be measured by comparison with other equivalent sequences. <a href="http://www.ncbi.nlm.nih.gov/pubmed/15217813">Bonnet <em>et al</em></a> have established that the majority of the pre-miRNA sequences exhibit a MFE that is lower than that for shuffled sequences.  We compute the probability that, for a given sequence, the MFE of the secondary structure is different from a distribution of MFE computed with 300 random sequences with the same length and the same dinucleotide frequency. </p>
            
            <p><b>Flag conserved mature miRNAs:</b> This option permits to check if the predicted miRNA belongs to some known miRNA family.  For that, we compare the sequence of the precursor with the database of mature miRNAs of plant (<em>Viridiplantae</em>) deposited in <a href="http://www.mirbase.org/ftp.shtml">miRBase</a> (Release 20). Alignments are performed with <a href="http://www.ebi.ac.uk/~guy/exonerate/">Exonerate</a>, which implements an exact model for pairwise alignment. We select alignments with at most three errors (mismatch, deletion or insertion) against the full-length mature miRNA and that occur in one of the two arms of the stem-loop. Moreover, this alignment allows to infer a putative location for the miRNA within the precursor.  This location is then validated with <a href="http://www.cs.mcgill.ca/~blanchem/mirdup/">miRdup</a>, that assesses the stability of the miRNA-miRNA* duplex. Here, it was trained on miRbase <em>Viridiplantae</em> V20.</p>


            <h3 id='submission'>Submission</h3>

            <p>Each job is automatically assigned an ID. </p>

            <p><b>Job title.</b> It's possible to identify the tool result by giving it a name. </p>

            <p><b>Email Address.</b> You can enter your email address to be notified when the job is finished. The email contains a link to access the results for 2 weeks.</p>
            
            
            <h2 id='results_page'>Result page</h2>
            
            <p>The page starts with a heading that summarizes the main results of the computation. Full results for known miRNAs (<a href="#known_mirnas">section 2.1</a>) and novel miRNAs (<a href="#new_mirnas">section 2.2</a>) are displayed in two new pages. </p>
            
            
            <h3 id="known_mirnas">Known miRNAs</h3>
            
            <p>Known miRNAs are miRNAs that are already present in the miRBase database (version 21). We consider that a known microRNA is found in the data as soon as there is at least one read  on the precursor sequence.</p>
            
            <p>The list of all known miRNAs found is displayed in a two-way table.</p>
            
            <p>Each row corresponds to a pre-miRNA, and each column to a feature. By default, results are sorted by sequence and then by position. It is possible to have them sorted by quality (see definition). You can view all information related to a given prediction by clicking on the row (see <a href="#html_report">section HTML Report</a>).</p> 
            
            <p><b>Name:</b> miRBase identifier</p>
            
            <p><b>Position:</b> Start and end positions of the pre-miRNA, as provided in miRBase</p>
            
            <p><b>+/- :</b> Strand, forward or reverse complement.</p>
            
            <p><b>Quality:</b> This score measures the consistency between the distribution of reads along the locus and the annotation provided in miRbase. It ranges between 0 and 2 stars, and  is calculated as follows: </p>
            
            <ul>
                <li>the locus contains more than 10 reads: add one star.</li>
                <li>more than half of the reads intersect either with the miRNA or the miRNA*: add one star. </li>
            </ul>

            <p><b>Reads:</b> number of reads included in the locus.</p>
                        
            <p><b>2D structure:</b> You can drag the mouse over the zoom icon to visualize the stem-loop structure of the pre-miRNA. The image is generated with <a href="http://varna.lri.fr/">Varna.</a></p>
                        
                        
            <h3 id="novel_mirnas">Novel miRNAs</h3> 
            
            <p>Novel miRNAs are miRNAS that are not reported in miRBase.</p> 
            
            <p>Each row corresponds to a pre-miRNA, and each column to a feature. By default, results are sorted by sequence and then by position. It is possible to have them sorted by quality (see definition) You can view all information related to a given prediction by clicking on the row (see <a href="#html_report">section HTML Report</a>).</p>
            
            <p><b>Name:</b> Name of the original sequence, as specified in the heading of the FASTA format.</p>
            
            <p><b>Position:</b> Start and end positions of the putative pre-miRNA in the original sequence in 1 based notation (consistently to the GFF format)</p>
            
            <p><b>+/- :</b> Strand, forward or reverse complement. </p>
            
            <p><b>Quality:</b> It is a combination of all other criteria described afterwards, and allows to rank the predictions according to the significance, from zero- to five- stars. It is calculated as follows.</p>
            
             <ul>
                <li><em>Thermodynamic stability of the hairpin precursor: MFEI &lt; -0.8.</em>  This MFEI threshold covers 83% of miRBase pre-miRNAs, whereas it is observed in less than 13% of pseudo hairpins (see XXX).</li>
                <li><em>Number of reads:</em> The  sequence has either at least 10 reads mapping to each arm, or  at least 5 reads mapping to each arm *and* at least 100 reads mapping in total. </li>
                <li><em>Presence of the miRNA:miRNA* duplex:</em> The most abundant reads from each arm of the precursor pair in the mature microRNA duplex with 0-4 nt overhang at their 3' ends.</li>
                <li><em>Precision of the precursor processing:</em>  At least 50% of reads mapping to each arm of the hairpin precursor have the same 5' end.</li>
                <li><em>Stability of the miRNA:miRNA* duplex:</em> The location of the mature miRNA is validated by miRdup.</li>
            </ul>           
            
            <p>Each criterion contributes equally to the overall ranking, and adds one star. Criteria on the number of reads, the presence of the miRNA:miRNA* duplex  and the precision of the precursor processing  are taken from miRBase 21, where they are used to select “high qulity” sequences (see <a href="http://nar.oxfordjournals.org/content/39/suppl_1/D152">http://nar.oxfordjournals.org/content/39/suppl_1/D152</a> and <a href="http://www.mirbase.org/blog/2014/03/high-confidence-micrornas/">http://www.mirbase.org/blog/2014/03/high-confidence-micrornas/</a>).</p>
            
            <p><b>Reads:</b> number of reads included in the locus</p>
            
            <p><b>MFE:</b> value of the minimal free energy of the secondary structure (computed with <a href="http://www.tbi.univie.ac.at/RNA/RNAeval.html">RNAeval</a>).</p>
            
            <p><b>MFEI:</b> value of the MFEI, as defined <a href="">here.</a></p>
            
            <p><b>Shuffles (option):</b> proportion of shuffled sequences whose MFE is lower than the MFE of the candidate miRNA precursor (see <a href="">Compute thermodynamic stability</a>).  This value ranges between 0 and 1. The smaller it is, the more significant is the MFE.  We report pre-miRNA stem-loops for which the value is smaller than 0.01, which covers more than 89% of miRBase sequences. Otherwise, if the P-value is greater than 0.01, we say that it is non significant, and do not report any value.</p>
            
            <p><b>Conservation  (option):</b> This cell is checked <img src='../style/check.png' alt='arobas' style='width:15px; height:15px;' /> when an alignment between the candidate sequence and miRBase is found (see <a href="">Flag conserved mature miRNAs</a>). It is doubled checked <img src='../style/check.png' alt='arobas' style='width:15px; height:15px;' /><img src='../style/check.png' alt='arobas' style='width:15px; height:15px;' /> when the location of the candidate mature miRNA is validated by miRdup. The alignments are visible in the HTML or ODF report.</p>
            
            <p><b>2D structure:</b> You can drag the mouse over the zoom icon to visualize the stem-loop structure of the pre-miRNA. The image is generated with <a href="http://varna.lri.fr/">Varna.</a></p>
            
            
            <h2 id='export'>Export </h2>
            
            
            <p>Results, or a selection of them, can be exported to a variety of formats, and saved to a local folder for further analyses. </p>

            <p id='gff'><b>GFF:</b> General annotation format, that displays the list of positions of pre-miRNA found (see more explanation on <a href="http://www.ensembl.org/info/website/upload/gff.html">Ensembl documentation</a>)</p>

            <p id='fasta'><b>FASTA:</b> This is the compilation of all pre-miRNA sequences found </p>

            <p id='dot_bracket'><b>Dot-bracket notation:</b> This is the compilation of all pre-miRNA sequences found, together with the predicted secondary structure. The secondary structure is given as a set of matching parentheses (see more explanation on <a href="">Vienna website</a>).  </p>

            <p id='csv'><b>CSV <em>(comma separated value)</em>:</b> It contains the same information as the result table, plus the FASTA sequences and the dot-bracket secondary structures. This tabular format is supported by spreadsheets like Excel. </p>

            <p id='odf'><b>ODF:</b> This is an equivalent of the <a href="#html_export">HTML report</a>, and contains the full report of the predictions. This document format is compatible with Word or OpenOffice. </p>    
                        

            <h2 id='html_export'>HTML Export </h2>
            
            
            <p>The HTML report contains all information related to a given predicted pre-miRNA.</p>
            
            
            <h3 id='html_known_mirnas'>HTML Report for known miRNAs</h3>
            
            <ul>
                <li><b>Name:</b> miRBase identifier, and link to access the miRBase entry.</li>
                <li><b>Position:</b> Start and end positions of the pre-miRNA, such as indicated in miRBase in 1-based notation (consistently with the GFF format). The length is indicated in parentheses.</li>
                <li><b>Strand:</b> + (forward) or - (reverse).</li>
                <li><b>GC content:</b> Percentage of bases that are either guanine or cytosine.</li>
                <li><b>Sequence (FASTA format):</b> Link to download the sequence.</li>
                <li><b>Stem-loop structure:</b> Link to download the secondary structure in <em>dot-bracket format</em>.  The first line contains a FASTA-like header. The second line contains the nucleic acid sequence. The last line contains the set of associated pairings encoded by brackets and dots. A base pair between bases <em>i</em> and <em>j</em> is represented by a '(' at position <em>i</em> and a ')' at position <em>j</em>. Unpaired bases are represented by dots. 
                    <p class='code'>
                    > Sample_1001-1085, stemloop structure <br />
                    cugagauacugccauagacgacuagccaucccucuggcucuuagauagccggauacagugauuuugaaagguuugugggguacag<br />
                    (((...((((.((((((((........(((.((((((((.......)))))))....).)))........)))))))))))))))<br />
                    </p>
                
                </li>
                <li><b>Number of reads:</b> The total number of reads included in the locus.</li>
                <li><b>Mapping profile:</b> Visual representation of all reads included in the locus.</li>
            </ul>
            
            
            <h3 id='html_novel_mirnas'>HTML Report for novel miRNAs</h3>
            
            <ul>
                <li><b>Name:</b> miRBase identifier, and link to access the miRBase entry.</li>
                <li><b>Position:</b> Start and end positions of the pre-miRNA, such as indicated in miRBase in 1-based notation (consistently with the GFF format). The length is indicated in parentheses.</li>
                <li><b>Strand:</b> + (forward) or - (reverse).</li>
                <li><b>GC content:</b> Percentage of bases that are either guanine or cytosine.</li>
                <li><b>Sequence (FASTA format):</b> Link to download the sequence.</li>
                <li><b>Stem-loop structure:</b> Link to download the secondary structure in <em>dot-bracket format</em>. The first line contains a FASTA-like header. The second line contains the nucleic acid sequence. The last line contains the set of associated pairings encoded by brackets and dots. A base pair between bases <em>i</em> and <em>j</em> is represented by a '(' at position <em>i</em> and a ')' at position <em>j</em>. Unpaired bases are represented by dots. 
                    <p class='code'>
                    > Sample_1001-1085, stemloop structure<br />
                    cugagauacugccauagacgacuagccaucccucuggcucuuagauagccggauacagugauuuugaaagguuugugggguacag<br />
                    (((...((((.((((((((........(((.((((((((.......)))))))....).)))........)))))))))))))))<br />
                    </p>
                </li>
                <li><b>Number of reads:</b> The total number of reads included in the locus.</li>
                <li><b>Mapping profile:</b> Visual representation of all reads included in the locus.</li>
                <li><b>Optimal MFE secondary structure:</b> If the stem-loop structure is not the MFE structure, we also provide a link to download the MFE structure.</li>
                <li><b>Alternative candidates (dot-bracket format):</b> This is the set of stem-loop sequences that overlap the current prediction. The choice between several alternative overlapping candidate pre-miRNAs is made according to the best MFEI.</li>
            </ul>
            
            <p>The stem-loop structure of the miRNA precursor is also displayed with <a href="http://varna.lri.fr/">Varna.</a></p>
            
            <img style='width:400px; display: block; margin: 0 auto;' src='../style/varna.png'' alt='Varna image' />
            
            <p><b>Thermodynamics stability</b></p>
            
            <ul>
                <li><b>MFE:</b> Value of the Minimum Free Energy (computed by <a href="http://www.tbi.univie.ac.at/RNA/RNAeval.html">RNAeval</a>) </li>
                <li><b>AMFE:</b> Value of the adjusted MFE : MFE/(sequence length) x 100</li>
                <li><b>MFEI:</b> Value of the minimum folding energy index (as defined <a href="">here</a>)</li>
                <li><b>Shuffles:</b> Proportion of shuffled sequences whose MFE is lower than the MFE of the candidate miRNA precursor (see <a href="#thermodynamic_stability">Compute thermodynamic stability</a>).  This value ranges between 0 and 1. The smaller it is, the more significant is the MFE.  We report pre-miRNA stem-loops for which the value is smaller than 0.01, which covers more than 89% of miRBase sequences. Otherwise, if the P-value is greater than 0.01, we say that it is non significant, and do not report any value.</li>
            </ul>
            
            <p><b>Conservation of the mature miRNA</b></p>
            
            <p>All alignments with miRBase are reported and gathered according to their positions. </p>
            
            <img style='width:610px; display: block; margin: 0 auto;' src='../style/alignment.png'' alt='alignment' />
            
            <p>query is the user sequence, and miRBase designates the mature miRNA found in miRBase. It is possible to access the corresponding mirBase entry by clicking on the link under the alignment. The report also indicates whether the location is validated with <a href="http://www.cs.mcgill.ca/~blanchem/mirdup/">miRdup</a>.  Finally, we provide an ASCII representation of the putative miRNA within the stem-loop  precursor.</p>
            
            <img style='width:600px; display: block; margin: 0 auto;' src='../style/hairpin_with_mature.png'' alt='hairpin with mature' />
            
            </div> <!-- main full -->
        
        </div><!-- bloc droit-->
        <?php include("../static/footer.txt") ?>
    </body>

</html>
