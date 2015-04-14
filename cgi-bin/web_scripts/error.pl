#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::WebTemplate;
use miRkwood::WebPaths;

my @css = (miRkwood::WebTemplate->get_server_css_file(), miRkwood::WebTemplate->get_css_file());
my @js  = (miRkwood::WebTemplate->get_js_file());

my $index_page = File::Spec->catfile( miRkwood::WebPaths->get_html_path(), 'index.php');

my $cgi    = CGI->new();
my $errorType   = $cgi->param('type');

my $errorMessage = "Unknown error";
if ( $errorType eq "noFasta" ){
    $errorMessage = <<"END_MSG";
    <p>An error occured when processing your data: miRkwood could not recognize the format of one or more sequences.</p>
    <p> Sequences should be in <b>FASTA</b> format. Lower-case and upper-case letters are both accepted.
    The full  standard I UPAC nucleic acid code is not supported : only <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt> and <tt>U</tt> symbols are recognized.</p>
    <br /><br />
END_MSG
}
elsif ( $errorType eq "severalSequences" ){
    $errorMessage = <<"END_MSG";
    <p>An error occured when processing your data:</p>
    <p> You must type only one sequence, in <b>FASTA</b> format. </p>
    <br /><br />    
    
    
END_MSG
}
elsif ( $errorType eq "tooLongSequence" ){
    $errorMessage = <<"END_MSG";
    <p>An error occured when processing your data:</p>
    <p>Your sequence must be shorter than 100000 nucleotides. </p>
    <br /><br />    
    
    
END_MSG
}
elsif ( $errorType eq "noBED" ){
    $errorMessage = <<"END_MSG";
    <p>An error occured when processing your data:</p>
    <p>Your BED file is not valid. Make sure you correctly used our provided script to convert BAM into BED file.</p>
    <br /><br />    
    
    
END_MSG
}
elsif ( $errorType eq "noGenome" ){
    $errorMessage = <<"END_MSG";
    <p>An error occured when processing your data:</p>
    <p>No genome file were found for the chosen species.</p>
    <br /><br />    
    
    
END_MSG
}
   
my $page = <<"END_TXT";
<div class="main">
  <div id="page">
    <h2>Error</h2>
      $errorMessage
    <div class="center">
      <input type="button" name="nom" value="Submit a new job" onclick="window.location.href='$index_page'" />
    </div>
  </div>
</div>
END_TXT

my $html = miRkwood::WebTemplate::get_HTML_page_for_content( 'static/', $page, \@css, \@js);

print <<"DATA" or miRkwood::WebTemplate::web_die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
###End###
