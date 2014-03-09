#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::WebTemplate;

my @css = (PipelineMiRNA::WebTemplate->get_server_css_file(), PipelineMiRNA::WebTemplate->get_css_file());
my @js  = (PipelineMiRNA::WebTemplate->get_js_file());

my $page = <<"END_TXT";
<div class="main">
  <form  name="form" onsubmit="return verifySequence();" onsubmit="wainting()" method="post" action="./results.pl" enctype="multipart/form-data">
    <fieldset id="fieldset">
      <div class="forms">
        <tr>
          <td class="label">
          &nbsp;<b> Job title </b> (optional)
          <input type="text" name="job" size="20">
          </td>
        </tr>
      </div>
      <div class="forms">
        <p class="label">
          <b>Enter query sequence</b>: Paste your RNA sequence(s) in FASTA format &nbsp;[<a href="./help.pl">?</a>]
        </p>
        <textarea id='seqArea' name="seqArea"  rows="10" cols="150" ></textarea>
      
        <p class="label">
            <p>or, upload a file <input type="file" name="seqFile" id="file" /> </p>
        </p>
        <br>
        <p><input class="checkbox" type="checkbox" name="strand" id="strand" value="strand"/>&#160;<label for='strand'>Scan both strands</label></p>
        <p><input class="checkbox" type="checkbox" name='CDS' id ="CDS" onclick="showHideBlock()"/>&#160;<label for='CDS'>Mask coding regions <i>(BlastX)</i></label>
        
        </p>
        <div id="menuDb">
          <p class="choixDiv" for="db">Choose organism database :</p>
          <p class="selectdb" ><select class="db" name="db">
            <option class="db" selected>Arabidopsis_thaliana</option>
            <option class="db">Oryza_sativa</option>
            <option class="db">Medicago_truncatula</option> 
          </select></p>
        </div>
        <p id='exempleClear'>
        <a id="area_clear" onclick="ResetForm();">clear</a> | <a id="seq_button"  onclick="generateExample();" />run with an example</a>
      </p>
      </div>
      <div class="forms">
        <p><b>Parameters</b>: Choose the annotation criteria for the miRNA precursors</p>
        <br>
        <P>
          
          <P><input class="checkbox" type="checkbox" checked="checked" name="mfei" id="mfei" value="mfeiChecked" />&#160;<label for='mfei'>Select only sequences with MFEI < -0.6</label></P>
          <P><input class="checkbox" type="checkbox" name="randfold" id="randfold" value="randfoldChecked" />&#160;<label for='randfold'>Compute thermodynamic stability <i>(shuffled sequences)</i></label></P>
          <P><input class="checkbox" type="checkbox" checked="checked" name="align" id="align" value="alignChecked" />&#160;<label for='align'>Flag conserved mature miRNAs <i>(alignment with miRBase + miRdup)</i></label></P>
       </div>
       <div class="forms">
         <tr>
           <td class="label">&nbsp; <b>E-mail address</b> (optional)
             <input type="text" name="mail" size="20">
           </td>
          </tr>
        </div>
        <div class="center">
           <input type="submit" name="upload" id="upload" value="Run miRkwood">
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
