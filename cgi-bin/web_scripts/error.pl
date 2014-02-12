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
    <h2>Error!</h2>
    <div class="section">
      miRkwood could not recognize the format of one or more sequences.
    </div>
    <ul>
      <li>Your data set should include <b>at least two</b> distinct sequences.</li>
      <li> Sequences should be in <b>FASTA</b> format.
        A sequence in FASTA format consists of a
        single-line description, followed by lines of sequence
        data. The first character of the description line is a
        greater-than (&quot;&gt;&quot;) symbol in the first
        column. All lines should be shorter than 80 characters. An
        example  in FASTA format is:
<pre>
> Name of the sequence 1
ctgcgagcgcgcgatgatagcgcggcgagcatgtagcatgctagctgtcgcgagcact
cggccgagatcaggcgatgcatgcgcagggagcagcgagcgacgagcacagcatgctagctagatgcatgctgtaggcagc
cgccgagagacgatggagctgc
> Name of the sequence 2
gacagatacgataagaggacgggatagaacgtagacatcgccgagagacgatggagctgc
cggccgagatcaggcgatgcatgcgcagggagcaggcgagcatgtagcatgctagctgtcgcgagcact
</pre>
      </li>
      <li><b>Lower-case and upper-case</b> letters are both accepted.</li>
      <li>The full  standard <b>IUPAC</b> nucleic acid code is not supported : only <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt> and <tt>U</tt> symbols are recognized.</li>
      <li><b>Numerical digits <tt>0</tt>, ..., <tt>9</tt>, <tt>-</tt> and dot <tt>.</tt></b> symbols are accepted. They are simply ignored by miRkwood.</li>
    </ul>
    <br /><br />
    <div class="home-made-center">
      <input type="button" name="nom" value="Submit a new job" onclick="window.location.href='./interface.pl'" />
    </div>
  </div>
</div>
END_TXT

my $html = PipelineMiRNA::WebTemplate::get_HTML_page_for_content($page, \@css, \@js);

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
###End###
