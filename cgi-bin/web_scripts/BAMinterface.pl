#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Spec;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::WebTemplate;

my @css = (
    miRkwood::WebTemplate->get_server_css_file(),
    miRkwood::WebTemplate->get_css_file(),
    miRkwood::WebTemplate->get_mirkwood_css_file()
);
my @js  = (
    miRkwood::WebTemplate->get_js_file(),
    miRkwood::WebTemplate->get_bioinfo_js_file()
);

my $help_page = File::Spec->catfile( File::Spec->catdir( miRkwood::WebPaths->get_html_path(), 'smallRNAseq'), 'help.php');

my $page = <<"END_TXT";
<div id="main">
    <br />
    <form id="form" onsubmit="return verifyBEDForm();" method="post" action="./BAMpipeline.pl" enctype="multipart/form-data">
        <fieldset id="fieldset">
            <div class="forms">
                <label for='job'><b>Job title</b> (optional)</label>
                <input type="text" id ='job' name="job" size="20"/>
            </div> 
            
            <div class="forms">
                <label for='bedFile'><b>Upload your set of reads </b>&nbsp;</label>
                <br /><br />
                This must be a BED file created by our script <i>mirkwood-bam2bed.pl [<a href="$help_page#input_form">?</a>]</i>
                <p><input type="file" name="bedFile" id="bedFile" /></p>
            </div> 
            
            <div class="forms">
                <b>Select an assembly </b>
                <br /><br />
                <label for="species">Select a reference genome in the list below [<a href="$help_page#reference">?</a>]</label><br />
                <p><select name="species" id="species">
                    <option value="" selected>--Choose assembly--</option>
                    <option value="Arabidopsis_lyrata">Arabidopsis lyrata (v1.0)</option>
                    <option value="Arabidopsis_thaliana">Arabidopsis thaliana (TAIR10)</option>
                    <option value="Brassica_napus">Brassica napus (Bnapus_V4.1)</option>
                    <option value="Brassica_rapa">Brassica rapa (Brapa_1.1)</option>
                    <option value="Glycine_max">Glycine max (v1.0.23)</option>
                    <option value="Lotus_japonicus">Lotus japonicus (Lj2.5)</option>
                    <option value="Medicago_truncatula">Medicago truncatula (JCVI_Mt3.5.2)</option>
                    <option value="Oryza_sativa">Oryza sativa (MSU7)</option>
                    <option value="Populus_trichocarpa">Populus trichocarpa (JGI_Poptr2.0)</option>
                    <option value="Solanum_lycopersicum">Solanum_lycopersicum (SL2.4)</option>
                    <option value="Sorghum_bicolor">Sorghum bicolor (Sorbi1)</option>
                    <option value="Vitis_vinifera">Vitis vinifera (Genoscope-20100122)</option>
                </select>  </p>
                <!-- <br /><br />
                <label for='seqArea'>Or paste your reference sequence in FASTA format &nbsp;[<a href="$help_page">?</a>]</label>
                <br /><br />
                <textarea id='seqArea' name="seqArea"  rows="10" cols="150" ></textarea>
                <br />
                <a id="area_clear" onclick="resetElementByID('seqArea')">Clear</a>
                <br /><br />
                <label for='seqFile'>Or upload a file: </label><input type="file" name="seqFile" id="seqFile" />
                <br />
                <a id="area_clear" onclick="resetElementByID('seqFile');">Clear Fasta</a> -->
            </div>
            
            <div class="forms">
                    <p>
                        <b>Parameters</b>
                    </p>
                    <div id='listParam'> 
                    <p>
                        <input class="checkbox" type="checkbox" checked="checked" name='CDS' id ="CDS"/>
                        &#160;<label for='CDS'>Mask coding regions</label>  [<a href="$help_page#mask_coding_regions">?</a>]
                    </p>
                    <p>
                        <input class="checkbox" type="checkbox" checked="checked" name='filter-tRNA-rRNA' id ="filter-tRNA-rRNA"/>
                        &#160;<label for='filter-tRNA-rRNA'>Filter out tRNA/rRNA/snoRNA</label>  [<a href="$help_page#filter_tRNA_rRNA">?</a>]
                    </p>
                    <p>
                        <input class="checkbox" type="checkbox" checked="checked" name="filter_multimapped" id="filter_multimapped" value="filter_multimappedChecked" />
                        &#160;<label for='filter_multimapped'>Remove multiply mapped reads ( &gt; 5 positions) </label>   [<a href="$help_page#filter_multimapped">?</a>]                       
                    </p>
                    <p>
                        <input class="checkbox" type="checkbox" checked="checked" name="mfei" id="mfei" value="mfeiChecked" />
                        &#160;<label for='mfei'>Select only sequences with MFEI &lt; -0.6</label>  [<a href="$help_page#filter_mfei">?</a>]
                    </p>
                    <p>
                        <input class="checkbox" type="checkbox" name="randfold" id="randfold" value="randfoldChecked" />
                        &#160;<label for='randfold'>Compute thermodynamic stability <i>(shuffled sequences)</i></label>  [<a href="$help_page#thermodynamic_stability">?</a>]
                    </p>
                    <p>
                        <input class="checkbox" type="checkbox" checked="checked" name="align" id="align" value="alignChecked" />
                        &#160;<label for='align'>Flag conserved mature miRNAs <i>(alignment with miRBase + miRdup)</i></label>  [<a href="$help_page#flag_conserved_mirnas">?</a>]
                    </p>
                    <p>
                        <input class="checkbox" type="checkbox" checked="checked" name='filter-bad-hairpins' id ="filter-bad-hairpins"/>
                        &#160;<label for='CDS'>Filter out low quality hairpins</label>  [<a href="$help_page#filter_bad_hairpins">?</a>]
                    </p>
                </div>
            </div>
            
            <div class="forms">
                <label for='mail'><b>E-mail address</b> (optional)</label>
                <input type="text" id='mail' name="mail" size="20"/>
            </div>
            
            <div class="center">
                <input type="reset" name="clear" id="clear" value="Clear form"/>
                <input type="submit" name="upload" id="upload" value="Run miRkwood"/>
            </div> 
        </fieldset>       
    </form>

</div><!-- main -->
END_TXT

my $title = 'miRkwood small RNA-seq';
my $html = miRkwood::WebTemplate::get_HTML_page_for_content('smallrnaseq', $page, \@css, \@js, $title);
print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
