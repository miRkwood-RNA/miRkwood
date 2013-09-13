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
    my $file = PipelineMiRNA::Paths->get_absolute_path('static', $file_name);
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

=method get_error_page

Return a generic error page

=cut

sub get_error_page {
    my @args = @_;
    my $error_message = shift @args;
    my $html = <<"HTML";
Content-type: text/html

<html>
    <head>
        <LINK rel="stylesheet" type="text/css" href="/arn/style/script.css" />
        <script src="/arn/js/miARN.js" type="text/javascript" LANGUAGE="JavaScript"></script>
        <title>MicroRNA identification</title>
    </head>
    <body>
        $error_message
    </body>
HTML
}

1;