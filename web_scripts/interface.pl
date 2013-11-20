#!/usr/bin/perl -w
use strict;
use warnings;

use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::WebTemplate;

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();

my $css = PipelineMiRNA::WebTemplate->get_css_file();
my $js  = PipelineMiRNA::WebTemplate->get_js_file();

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<meta name="keywords" content="RNA, ARN, mfold, fold, structure, prediction, secondary structure" />
<link title="test" type="text/css" rel="stylesheet" href="$css" />

        <script src="$js" type="text/javascript" LANGUAGE="JavaScript"></script>
       
        <title>miREST</title>


</head>
<body>
<div class="theme-border"></div>
<div class="logo"></div>

$bioinfo_menu

<div class="bloc_droit">

$header_menu

<div class="main">

        <form  name="form" onsubmit="verifySequence()" onsubmit="wainting()" method="post" action="./results.pl" enctype="multipart/form-data">
            <fieldset id="fieldset">    
				 <div class="forms">
				<tr>
					<td class="label"> 
						Enter a name for your job (optional)</i>: 
						<input type="text" name="job" size="20">
					</td>
				  </tr>
				</div>
				<div class="forms">
					<label class="textDiv" for="seqArea">
					   <legend><h2>Paste your sequence(s)[<a href="./help.pl">?</a>]:</h2></legend>
				    </label><br />
					<textarea id='seqArea' name="seqArea"  rows="10" cols="150" ></textarea>
					<label class="textDiv" style="color:#333" for="seq">Or, upload file :</label><br />
					<input type="file" name="seqFile" id="file" />
					<input id="seq_button" type="button" value="Example" onclick="generateExample();" />
				</div>
				  <div class="forms">
					 <label class="checkbox" for="db"><h2>Mask coding regions [<a href="./help.pl">?</a>] <i><small>(may be slow) </small></i> :
					<input  id ="CDS" type="checkbox" name="check" value="checked" onclick="showHideBlock()"> </h2>               
					<div id="menuDb"> 
					<label class="choixDiv" for="db">Choose organism database :</label>
						<select class="db" name="db">
							<option class="db" selected>ATpepTAIR10</option>
							<option class="db">plante</option>
						</select>
					</div>
               </div>
               <div class="forms">
				<h2>Select additional features:</h2>
				<P>
               <P > <input class="checkbox" type="checkbox" checked="checked" name="randfold" value="randfoldChecked">Compute thermodynamic stability [<a href="./help.pl">?</a>]</input></P>
                  
               <P > <input class="checkbox" type="checkbox" checked="checked" name="mfei" value="mfeiChecked">Compute MFE/MFEI/AMFE (minimal folding energy)[<a href="./help.pl">?</a>]</input></P>
             
               <P > <input  class="checkbox" type="checkbox" checked="checked" name="align" value="alignChecked">Align against mature microRNAs miRBase [<a href="./help.pl">?</a>]</input></P>
                
                </div>
               </P>
                
                </label>
				<div class="forms">
				<tr>
					<td class="label"> 
						Enter your <b>E-mail</b> address <i>(optional)</i>: 
						<input type="text" name="mail" size="20">
					</td>
				  </tr>
				</div>
              <input style="margin-left:632px" type="submit" name="upload" id="upload" value="Submit JOB">
                
            
            
            </fieldset>
        
        </form>
        
        
        
        
        
        
</div><!-- main -->

$footer

</div><!-- bloc droit-->

</body>
</html>

DATA
###End###
