#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::Components;
use PipelineMiRNA::WebTemplate;
use PipelineMiRNA::Utils;

my $cgi = CGI->new();

my $name     = $cgi->param('nameSeq');
my $factor   = $cgi->param('factor');
my $value    = $cgi->param('value');
my $position = $cgi->param('position');
my $typePage = $cgi->param('typePage');
my $url      = $cgi->param('url');

my @css = (PipelineMiRNA::WebTemplate->get_server_css_file(), PipelineMiRNA::WebTemplate->get_css_file());
my @js  = (File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'results.js'));

my $page;
if ( $typePage eq 'simpleCell' ) {
    $page = <<"DATA";
    <div class="main main-full">
        <div class="titreDiv"> MicroRNA identification results:</div>
        <div id = 'showInfo'>
            <h2 ><u>Sequence Informations </u></h2><br/>
            <li><b>Name sequence:</b> ${name}$position
            </li>
            <li><b>$factor:</b> $value
            </li>
        </div>
    </div>
DATA
}
elsif ( $typePage eq 'alignement' ) {

    my $error = "Deprecated feature.";
    print PipelineMiRNA::WebTemplate::get_error_page($error);
    die($error);

}
elsif ( $typePage eq 'image' ) {
    my @url = split( /\//xms, $url );
    my $image_name = $url[5];
    $page = <<"DATA";
<div class="main main-full">
  <div class="titreDiv"> MicroRNA identification results:</div>
  <h2>Structure:</h2>
  <div class="figure" >
    <img src='$url' border=0 alt='image'/>
	<p>Fig : $image_name sequence</p>
  </div>
</div>
DATA
}

my $html = PipelineMiRNA::WebTemplate::get_HTML_page_for_content($page, \@css, \@js, 1);

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
###End###
