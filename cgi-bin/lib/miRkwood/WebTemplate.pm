package miRkwood::WebTemplate;

# ABSTRACT: The HTML templates (or bits of) used by the web interface

use strict;
use warnings;

use File::Spec;
use miRkwood::Paths;
use miRkwood::WebPaths;

=method get_static_file

Return the contents of a given file in the stati directory

=cut

sub get_static_file {
    my @args = @_;
    my $file_name = shift @args;
    my $file = File::Spec->catfile(miRkwood::WebPaths->get_static_path(),
                                   $file_name);
    open my $FILE, '<', $file;
    my $contents = do { local $/; <$FILE> };
    close $FILE;
    return $contents;
}

=method get_bioinfo_menu

Return the HTML BioInfo left menu

=cut

sub get_bioinfo_menu {
    my @args = @_;
    return get_static_file('bioinfo_menu.txt');
}

=method get_header_menu

Return the HTML BioInfo header menu

=cut

sub get_header_menu {
    my @args = @_;
    my $pipeline = shift @args;
    if ( $pipeline eq 'smallrnaseq' ){
        return get_static_file('../smallRNAseq/header_menu.txt');
    }
    elsif ( $pipeline eq 'abinitio' ){
        return get_static_file('../abinitio/header_menu.txt');
    }
    return get_static_file('header_menu.txt');
}

=method get_footer

Return the HTML BioInfo footer

=cut

sub get_footer {
    my @args = @_;
    return get_static_file('footer.txt');
}

=method get_link_back_to_results

Return the link back to the main results page

=cut

sub get_link_back_to_results {
    my @args = @_;
    my $jobId = shift @args;
    return ("./resultsWithID.pl?run_id=$jobId");
}

=method get_link_back_to_BAM_results

Return the link back to the main results page
for BAM pipeline

=cut

sub get_link_back_to_BAM_results {
    my @args = @_;
    my $jobId = shift @args;
    return ("./BAMresults.pl?run_id=$jobId");
}

=method get_link_back_to_BED_known_results

Return the link back to the known results page
for BED pipeline

=cut

sub get_link_back_to_BED_known_results {
    my @args = @_;
    my $jobId = shift @args;
    return ("./BAMresults_for_mirnas.pl?jobID=$jobId&type=Known");
}

=method get_link_back_to_BED_new_results

Return the link back to the new results page
for BED pipeline

=cut

sub get_link_back_to_BED_new_results {
    my @args = @_;
    my $jobId = shift @args;
    return ("./BAMresults_for_mirnas.pl?jobID=$jobId&type=New");
}

=method get_css_file

Return the project CSS file

=cut

sub get_css_file {
    my @args = @_;
    return File::Spec->catfile(miRkwood::WebPaths->get_css_path(), 'script.css');
}

=method get_mirkwood_css_file

Return the miRkwood CSS file

=cut

sub get_mirkwood_css_file {
    my @args = @_;
    return File::Spec->catfile(miRkwood::WebPaths->get_css_path(), 'mirkwood.css');
}

=method get_server_css_file

Return the server CSS file

=cut

sub get_server_css_file {
    my @args = @_;
    return File::Spec->catfile(miRkwood::WebPaths->get_server_css_path, 'bioinfo.css');
}


=method get_js_file

Return the main JavaScript file

=cut

sub get_js_file {
    my @args = @_;
    return File::Spec->catfile(miRkwood::WebPaths->get_js_path(), 'miARN.js');
}


=method get_bioinfo_js_file

Return the bioinfo JavaScript file

=cut

sub get_bioinfo_js_file {
    my @args = @_;
    return File::Spec->catfile(miRkwood::WebPaths->get_js_path(), 'header.js');
}


=method get_jquery_lib

Return the path to jquery lib on Bioinfo

=cut
sub get_jquery_lib {
    my @args = @_;
    return File::Spec->catfile(miRkwood::WebPaths->get_server_libs_path(), 'jquery-1.11.3.min.js');
}


=method get_error_page

Return a generic error page

=cut

sub get_error_page {
    my @args = @_;
    my $error_message = shift @args;
    my @css = (get_server_css_file(), get_css_file());
    my @js  = (get_js_file());
    my $header = 'Sorry, something went wrong with miRkwood';
    my $explanation = 'The error which occured is:';
    my $footer = 'Please send this to the miRkwood team, at the address in the footer.';
    my $contents = "<br/><br/>$header<br/><br/>$explanation<br/><br/>$error_message<br/><br/><br/>$footer";
    my $title = 'miRkwood - Error';
    my $html = miRkwood::WebTemplate::get_HTML_page_for_content( 'static', $contents, \@css, \@js, $title);
    my $res = <<"HTML";
Content-type: text/html

$html
HTML
    return $res;
}

=method web_die

Die in a web context

=cut

sub web_die {
    my @args = @_;
    my $error_message = shift @args;
    print get_error_page($error_message);
    die($error_message);
    return;
}


=method get_cgi_url

Make an URL to the given CGI page

=cut

sub get_cgi_url {
    my @args = @_;
    my $page = shift @args;
    # dirname( $ENV{HTTP_REFERER} );
    my $path = File::Spec->catfile($ENV{SERVER_NAME}, miRkwood::WebPaths->get_web_scripts(), $page);
    my $url  = 'http://'. $path;
    return $url;
}


sub get_HTML_page_for_content {
    my @args      = @_;
    my $pipeline  = shift @args;    # should be 'abinitio' or 'smallrnaseq'
    my $page      = shift @args;
    my $css_files = shift @args;
    my $js_files  = shift @args;
    my $title     = shift @args;
    my $no_menu   = shift @args;

    #~ my $bioinfo_menu = '';
    #~ if (! $no_menu){
        #~ $bioinfo_menu = miRkwood::WebTemplate::get_bioinfo_menu();
    #~ }

    my $header_menu  = miRkwood::WebTemplate::get_header_menu($pipeline);
    my $footer       = miRkwood::WebTemplate::get_footer();

    my $index_link;
    if ( $pipeline eq 'static' ){
        $index_link = '/mirkwood/index.php';
    }
    elsif ( $pipeline eq 'abinitio' ){
        $index_link = '/mirkwood/abinitio/index.php';
    }
    else{
        $index_link = '/mirkwood/smallRNAseq/index.php';
    }

    my $body = <<"END_TXT";
    <body>
        <div class="frametitle">
            <h1 id="title">$title</h1>
        </div>

        <div id="center_sup">
            <div id="link_home" style="display:inline-block"><a href="$index_link" class="text_onglet"><img src="/Style/icon/home_w.png" alt="home_general"/></a></div>
            <div class="tabs" id="menu_central" style="display:inline-block"> 
                $header_menu
            </div>
            <div id="arborescence"></div>
        </div>
        <div class="bloc_droit">
        
            $page
        </div><!-- bloc droit-->
        $footer
    </body>
END_TXT
    my $HTML = get_HTML_page_for_body($body, $css_files, $js_files, $title);
    return $HTML;
}

sub get_HTML_page_for_body {
    my @args      = @_;
    my $body      = shift @args;
    my @css_files = @{shift @args};
    my @js_files  = @{shift @args};
    my $title     = shift @args;

    if ( ( ! defined($title) ) || $title eq '' ){
        $title = 'miRkwood - MicroRNA identification';
    }

    my $css_html = '';
    foreach my $css (@css_files){
        $css_html .= "<link type='text/css' rel='stylesheet' href='$css' />\n";
    }

    my $js_html = '';
    foreach my $js (@js_files){
        $js_html .= "<script type='text/javascript' src='$js'></script>\n";
    }

    my $jquery_html = "<script type='text/javascript' src='" . get_jquery_lib() . "'></script>\n";

    my $HTML = <<"END_TXT";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
        <meta name="keywords" content="microRNA, miRNA, premir, plant, arabidopsis thaliana, mirkwood, RNAfold" />
        $jquery_html
        <title>$title</title>
        $css_html        $js_html    </head>
    $body
</html>
END_TXT

    return $HTML;
}

sub send_email {
    my @args  = @_;
    my $mode  = shift @args;    # should be 'abinitio' or 'smallRNAseq'
    my $to    = shift @args;
    my $jobId = shift @args;
    my $title = ' ';
    if (@args) {
        $title = shift @args;
    }
    require MIME::Lite;

    my $results_page  = '';
    if ( $mode eq 'smallRNAseq' ){
        $results_page = 'BAMresults.pl';
    }
    else{
        $results_page  = 'resultsWithID.pl';
    }
    my $results_baseurl = get_cgi_url($results_page);

    my $from = 'mirkwood@univ-lille1.fr';
    my $subject = "[miRkwood] Results for job $jobId";

    my $res_arguments = "?run_id=$jobId";
    my $results_url   = $results_baseurl . $res_arguments;

    my $msg;
    if ($title){
        $msg = "Your job \"$title\" is completed.";
    } else {
        $msg = 'Your job is completed.';
    }


    my $message = <<"DATA";
Dear miRkwood user,

$msg

Results are available at
$results_url

Thank you for using miRkwood!

-- 
The miRkwood team
DATA

    my $email = MIME::Lite->new(
                 From     => $from,
                 To       => $to,
                 Subject  => $subject,
                 Data     => $message
                 );
    $email->send('sendmail');
    return;
}

1;
