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
  <h2>Retrieve result with an ID</h2>
  <div class="forms">
    <form method="post" action="./resultsWithID.pl">
      <p>The ID remains valid 15 days after sequence submission.</p>
      <label for='run_id'><b>Enter the ID:</b></label>
        <input type="text" name="run_id" id='run_id' size="20">
        <input type="hidden" name="command" value="result">
        <input type="submit" value="Go">
      </p>
    </form>
  </div> <!-- form -->
</div>
END_TXT

my $html = PipelineMiRNA::WebTemplate::get_HTML_page_for_content($page, \@css, \@js);

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
###End###
