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
          Enter a <b> name </b>for your job  <i>(optional)</i></i>:
          <input type="text" name="job" size="20">
          </td>
        </tr>
      </div>
      <div class="forms">
        <p class="label">
          <b>Paste</b> your RNA sequences in FASTA format &nbsp;&nbsp;[<a href="./help.pl">?</a>]
        </p>
        <textarea id='seqArea' name="seqArea"  rows="10" cols="150" ></textarea>
        <p>or</p>
        <p class="label">
          <b>upload</b> a file
          <input type="file" name="seqFile" id="file" />
        </p>
        <input id="seq_button" type="button" value="Example" onclick="generateExample();" />
        <p><input class="checkbox" type="checkbox" checked="checked" name="strand" value="strand">Process both strand </input></p>
        <p class="label">
          <b>Mask coding regions [<a href="./help.pl">?</a>] <i><small>(may be slow) </small></i></b> :
          <input  id ="CDS" type="checkbox" name="check" value="checked" onclick="showHideBlock()">
        </p>
        <div id="menuDb">
          <p class="choixDiv" for="db">Choose organism database :</p>
          <select class="db" name="db">
            <option class="db" selected>ATpepTAIR10</option>
            <option class="db">plante</option>
          </select>
        </div>
      </div>
      <div class="forms">
        <p><b>Select additional features</b>:</p>
        <P>
          <P><input class="checkbox" type="checkbox" checked="checked" name="randfold" value="randfoldChecked">Compute thermodynamic stability [<a href="./help.pl">?</a>]</input></P>
          <P><input class="checkbox" type="checkbox" checked="checked" name="mfei" value="mfeiChecked">Compute MFE/MFEI/AMFE (minimal folding energy)[<a href="./help.pl">?</a>]</input></P>
          <P><input class="checkbox" type="checkbox" checked="checked" name="align" value="alignChecked">Align against mature microRNAs miRBase [<a href="./help.pl">?</a>]</input></P>
       </div>
       <div class="forms">
         <tr>
           <td class="label">Enter your <b>E-mail</b> address <i>(optional)</i>:
             <input type="text" name="mail" size="20">
           </td>
          </tr>
        </div>
        <div class="center">
           <input type="submit" name="upload" id="upload" value="Run MiRNA">
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
