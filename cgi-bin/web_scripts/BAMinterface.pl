#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Spec;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::WebTemplate;

my @css = (miRkwood::WebTemplate->get_server_css_file(), miRkwood::WebTemplate->get_css_file());
my @js  = (miRkwood::WebTemplate->get_js_file());

my $help_page = File::Spec->catfile( miRkwood::WebPaths->get_html_path(), 'help.php');

my $page = <<"END_TXT";
<div class="main">
    <form id="form" onsubmit="return verifyBEDForm();" method="post" action="./BAMpipeline.pl" enctype="multipart/form-data">
        <fieldset id="fieldset">
            <div class="forms">
                <label for='job'><b>Job title</b> (optional)</label>
                <input type="text" id ='job' name="job" size="20"/>
            </div> 
            
            <div class="forms">
                <label for='bedFile'><b>Upload your BED file: </b>&nbsp;[<a href="$help_page">?</a>]</label>
                <input type="file" name="bedFile" id="bedFile" />
            </div> 
            
            <div class="forms">
                <label for="species"><b>Model organism: </b>Choose an organism in the list below: </label><br />
                <select name="species" id="species">
                    <option value="" selected>--Choose species--</option>
                    <option value="Arabidopsis_thaliana">Arabidopsis thaliana</option>
                </select>  
                <br />
                <p>
                    <label for='seqArea'>Or paste your reference sequence in FASTA format &nbsp;[<a href="$help_page">?</a>]</label>
                </p>
                <textarea id='seqArea' name="seqArea"  rows="10" cols="150" ></textarea>
                <br /><a id="area_clear" onclick="resetTextarea('seqArea')">Clear</a> 
            </div>
            
            <div class="forms">
                <b>Parameters:</b>
                <p>
                    <input class="checkbox" type="checkbox" name="strand" id="strand" value="strand"/>
                    <label for='strand'>Scan both strands</label> [<a href="$help_page#scan-both-strands">?</a>]
                </p>
                <p>
                    <input class="checkbox" type="checkbox" name='filter-tRNA-rRNA' id ="filter-tRNA-rRNA"/>
                    <label for='filter-tRNA-rRNA'>Filter out tRNA/rRNA <i>(tRNAscan-SE / RNAmmer)</i></label>  [<a href="$help_page#filter_tRNA_rRNA">?</a>]
                </p>
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

my $html = miRkwood::WebTemplate::get_HTML_page_for_content($page, \@css, \@js);
print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
