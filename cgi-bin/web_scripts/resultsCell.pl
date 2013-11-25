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

my $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
my $footer       = PipelineMiRNA::WebTemplate::get_footer();

my $bioinfo_css = PipelineMiRNA::WebTemplate->get_server_css_file();
my $project_css = PipelineMiRNA::WebTemplate->get_css_file();
my $js  = File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'results.js');

=method make_HTML

Returns the HTML page with the given <body>

=cut

sub make_HTML {
    my @args = @_;
    my $body = shift @args;
    my $html = <<"DATA";
Content-type: text/html

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html>
    <head>
        <link title="test" type="text/css" rel="stylesheet" href="$project_css" />
        <link title="test" type="text/css" rel="stylesheet" href="$bioinfo_css" />
        <LINK rel="stylesheet" type="text/css" href="$css" />
        <script type="text/javascript" language="Javascript" src="$js"> </script>
        <title>MicroRNA identification</title>
    </head>
    <body>
        <div class="theme-border"></div>
        <div class="logo"></div>
	$bioinfo_menu
        <div class="bloc_droit">

            $header_menu
            <div class="main main-full">
                $body
            </div><!-- main -->
        </div><!-- bloc droit-->
	$footer
    </body>
</html>
DATA
    return $html;
}

# Script code

my $body;
if ( $typePage eq 'simpleCell' ) {
    $body = <<"DATA";
        <div class="titreDiv"> MicroRNA identification results:</div>
        <div id = 'showInfo'>
            <h2 ><u>Sequence Informations </u></h2><br/>
            <li><b>Name sequence:</b> ${name}$position
            </li>
            <li><b>$factor:</b> $value
            </li>
        </div>
DATA
    my $html = make_HTML($body);
}
elsif ( $typePage eq 'alignement' ) {

    my $error = "Deprecated feature.";
    print PipelineMiRNA::WebTemplate::get_error_page($error);
    die($error);

}
elsif ( $typePage eq 'image' ) {
    my @url = split( /\//xms, $url );
    my $image_name = $url[5];
    $body = <<"DATA";
		<div class="titreDiv"> MicroRNA identification results:</div>
		<h2> Structure :</h2>
		<div class="figure" >
		  <img src='$url' border=0 alt='image'>
		  <p>Fig : $image_name sequence
    	</div>	
DATA
}

print make_HTML($body) or die("Error when displaying HTML: $!");

###End###