#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Spec;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::WebTemplate;

my @css = (PipelineMiRNA::WebTemplate->get_server_css_file(), PipelineMiRNA::WebTemplate->get_css_file());
my @js  = (PipelineMiRNA::WebTemplate->get_js_file());

my $help_page = File::Spec->catfile( PipelineMiRNA::WebPaths->get_css_path(), '..', 'help.php');

my $page = <<"END_TXT";
<div class="main">
  <form id='form' onsubmit="return verifySequence();" method="post" action="./results.pl" enctype="multipart/form-data">
    <fieldset id="fieldset">
      <div class="forms">
          <label for='job'>&nbsp;<b>Job title</b> (optional)</label>
          <input type="text" id ='job' name="job" size="20"/>
      </div>
      <div class="forms">
        <p>
          <label for='seqArea'><b>Enter query sequence</b>: Paste your RNA sequence(s) in FASTA format &nbsp;[<a href="$help_page#query_sequence">?</a>]</label>
        </p>
        <textarea id='seqArea' name="seqArea"  rows="10" cols="150" ></textarea>

        <p>
            <label for='seqFile'>or, upload a file</label><input type="file" name="seqFile" id="seqFile" />
        </p>
        <br/>
        <p><input class="checkbox" type="checkbox" name="strand" id="strand" value="strand"/>
            &#160;<label for='strand'>Scan both strands</label> [<a href="$help_page#scan-both-strands">?</a>]</p>
        <p><input class="checkbox" type="checkbox" name='CDS' id ="CDS" onclick="showHideBlock()"/>
            &#160;<label for='CDS'>Mask coding regions <i>(BlastX)</i></label>  [<a href="$help_page#mask-coding-regions">?</a>]</p>
        <div id="menuDb">
          <label class="choixDiv selectdb" for="db">Choose organism database:</label>
          <select class="db" name="db" id='db'>
            <option class="db" selected="selected">Arabidopsis_thaliana</option>
            <option class="db">Oryza_sativa</option>
            <option class="db">Medicago_truncatula</option> 
          </select>
        </div>
        <p><input class="checkbox" type="checkbox" name='filter-tRNA-rRNA' id ="filter-tRNA-rRNA"/>
            &#160;<label for='filter-tRNA-rRNA'>Filter out tRNA/rRNAs <i>(tRNAscan-SE / RNAmmer)</i></label>  [<a href="$help_page#filter_tRNA_rRNA">?</a>]</p>
        <p id='exempleClear'>
        <a id="area_clear" onclick="ResetForm();">clear</a> | <a id="seq_button"  onclick="generateExample();">run with an example</a>
      </p>
      </div>
      <div class="forms">
        <p><b>Parameters</b>: Choose the annotation criteria for the miRNA precursors [<a href="$help_page#parameters">?</a>]</p>
        <div id='listParam'>  <p><input class="checkbox" type="checkbox" checked="checked" name="mfei" id="mfei" value="mfeiChecked" />
            &#160;<label for='mfei'>Select only sequences with MFEI &lt; -0.6</label></p>
          <p><input class="checkbox" type="checkbox" name="randfold" id="randfold" value="randfoldChecked" />
            &#160;<label for='randfold'>Compute thermodynamic stability <i>(shuffled sequences)</i></label></p>
          <p><input class="checkbox" type="checkbox" checked="checked" name="align" id="align" value="alignChecked" />
            &#160;<label for='align'>Flag conserved mature miRNAs <i>(alignment with miRBase + miRdup)</i></label></p>
        </div></div>
       <div class="forms">
         <label for='mail'>&nbsp;<b>E-mail address</b> (optional)</label>
         <input type="text" id='mail' name="mail" size="20"/>
        </div>
        <div class="center">
           <input type="submit" name="upload" id="upload" value="Run miRkwood"/>
        </div>
    </fieldset>
  </form>
</div><!-- main -->
END_TXT

my $html = PipelineMiRNA::WebTemplate::get_HTML_page_for_content($page, \@css, \@js);
print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
###End###
