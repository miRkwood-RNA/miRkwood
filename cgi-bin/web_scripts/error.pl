#!/usr/bin/perl -w
use strict;
use warnings;

use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::WebTemplate;

my @css = (PipelineMiRNA::WebTemplate->get_server_css_file(), PipelineMiRNA::WebTemplate->get_css_file());
my @js  = (PipelineMiRNA::WebTemplate->get_js_file());

my $page = <<"END_TXT";
<div class="main">
  <div id="page">
    <h2>Error</h2>
      <p>An error occured when processing your data: miRkwood could not recognize the format of one or more sequences.</p>
      <p> Sequences should be in <b>FASTA</b> format. Lower-case and upper-case letters are both accepted.
      The full  standard I UPAC nucleic acid code is not supported : only <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt> and <tt>U</tt> symbols are recognized.</p>
    <br /><br />
    <div class="center">
      <input type="button" name="nom" value="Submit a new job" onclick="window.location.href='./interface.pl'" />
    </div>
  </div>
</div>
END_TXT

my $html = PipelineMiRNA::WebTemplate::get_HTML_page_for_content($page, \@css, \@js);

print <<"DATA" or PipelineMiRNA::WebTemplate::web_die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
###End###
