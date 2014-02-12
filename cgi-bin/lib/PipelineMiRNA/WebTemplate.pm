package PipelineMiRNA::WebTemplate;

# ABSTRACT: The HTML templates (or bits of) used by the web interface

use strict;
use warnings;

use File::Spec;
use PipelineMiRNA::Paths;

=method get_static_file

Return the contents of a given file in the stati directory

=cut

sub get_static_file {
    my @args = @_;
    my $file_name = shift @args;
    my $file = File::Spec->catfile(PipelineMiRNA::Paths->get_static_path(),
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

=method get_css_file

Return the project CSS file

=cut

sub get_css_file {
    my @args = @_;
    return File::Spec->catfile(PipelineMiRNA::Paths->get_css_path(), 'script.css');
}

=method get_server_css_file

Return the server CSS file

=cut

sub get_server_css_file {
    my @args = @_;
    return File::Spec->catfile(PipelineMiRNA::Paths->get_server_css_path, 'bioinfo.css');
}


=method get_js_file

Return the main JavaScript file

=cut

sub get_js_file {
    my @args = @_;
    return File::Spec->catfile(PipelineMiRNA::Paths->get_js_path(), 'miARN.js');
}

=method get_error_page

Return a generic error page

=cut

sub get_error_page {
    my @args = @_;
    my $error_message = shift @args;
    my $css = get_css_file();
    my $js  = get_js_file();
    my $html = <<"HTML";
Content-type: text/html

<html>
    <head>
        <LINK rel="stylesheet" type="text/css" href="$css" />
        <script src="$js" type="text/javascript" LANGUAGE="JavaScript"></script>
        <title>MicroRNA identification</title>
    </head>
    <body>
        $error_message
    </body>
HTML
}

=method make_mirbase_link

Return the URL to MirBase given the identifier

=cut

sub make_mirbase_link {
    my @args = @_;
    my $id   = shift @args;
    my $url  = 'http://mirbase.org/cgi-bin/mirna_entry.pl?acc=';
    return $url . $id;
}

=method make_url

Make an URL to the given page based on the server root

=cut

sub make_url {
    my @args = @_;
    my $page = shift @args;
    # dirname( $ENV{HTTP_REFERER} );
    my $path = File::Spec->catfile($ENV{SERVER_NAME}, PipelineMiRNA::Paths->get_web_root(), $page);
    my $url  = 'http://'. $path;
    return $url;
}

sub get_simple_results_page {
    my @args = @_;
    my $page = shift @args;
    my $css  = shift @args;
    my $HTML = <<"END_TXT";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <title>MiRest - MicroRNA identification</title>
        <STYLE type="text/css">$css</STYLE>
   </head>
    <body>
        $page
    </body>
</html>
END_TXT
    return $HTML
}

sub get_HTML_page_for_content {
    my @args      = @_;
    my $page      = shift @args;
    my $css_files = shift @args;
    my $js_files  = shift @args;
    my $no_menu   = shift @args;

    my $bioinfo_menu = '';
    if (! $no_menu){
        $bioinfo_menu = PipelineMiRNA::WebTemplate::get_bioinfo_menu();
    }

    my $header_menu  = PipelineMiRNA::WebTemplate::get_header_menu();
    my $footer       = PipelineMiRNA::WebTemplate::get_footer();

    my $body = <<"END_TXT";
    <body>
        <div class="theme-border"></div>
        <div class="logo"></div>
        $bioinfo_menu
        <div class="bloc_droit">
        $header_menu
            $page
        $footer
        </div><!-- bloc droit-->
    </body>
END_TXT
    my $HTML = get_HTML_page_for_body($body, $css_files, $js_files);
    return $HTML;
}

sub get_HTML_page_for_body {
    my @args      = @_;
    my $body      = shift @args;
    my @css_files = @{shift @args};
    my @js_files  = @{shift @args};

    my $css_html = '';
    foreach my $css (@css_files){
        $css_html .= "<link type='text/css' rel='stylesheet' href='$css' />";
    }

    my $js_html = '';
    foreach my $js (@js_files){
        $js_html .= "<script type='text/javascript' src='$js'></script>";
    }

    my $HTML = <<"END_TXT";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
        <meta name="keywords" content="RNA, ARN, mfold, fold, structure, prediction, secondary structure" />
        <title>miREST - MicroRNA identification</title>
        $css_html
        $js_html
    </head>
    $body
</html>
END_TXT

    return $HTML;
}

1;