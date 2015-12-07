#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin;
use File::Spec;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::WebTemplate;

my @css = (miRkwood::WebTemplate->get_server_css_file(), miRkwood::WebTemplate->get_css_file());
my @js  = ();

my $page = <<"END_TXT";
<div class="main">
    Be patient, I'm currently working on it ;-)

</div><!-- main -->
END_TXT

my $title = 'miRkwood - under construction';
my $html = miRkwood::WebTemplate::get_HTML_page_for_content( 'static', $page, \@css, \@js, $title);
print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
