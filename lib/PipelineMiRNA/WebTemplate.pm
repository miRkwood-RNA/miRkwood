package PipelineMiRNA::WebTemplate;

use strict;
use warnings;

use File::Spec;
use PipelineMiRNA::Paths;

sub get_static_file {
    my @args = @_;
    my $file_name = shift @args;
    my $file = PipelineMiRNA::Paths->get_absolute_path('static', $file_name);
    open my $FILE, '<', $file;
    my $contents = do { local $/; <$FILE> };
    close $FILE;
    return $contents;
}

sub get_bioinfo_menu {
    my @args = @_;
    return get_static_file('bioinfo_menu.txt');
}

sub get_header_menu {
    my @args = @_;
    return get_static_file('header_menu.txt');
}

sub get_footer {
    my @args = @_;
    return get_static_file('footer.txt');
}

1;