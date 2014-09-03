package miRkwood::FileUtils;

# ABSTRACT: library of utility methods around files

use strict;
use warnings;

use File::Spec;

=method is_directory

Return whether the given file object (name and parent)
is a directory or not.

=cut

sub is_directory {
    my @args = @_;
    my ($directory_name, $parent_directory) = @args;
    my $directory = File::Spec->catdir( $parent_directory, $directory_name );
    return ( $directory_name ne '.'
             && $directory_name ne '..'
             && -d $directory );
}

=method get_dirs_from_directory

Get the sub-directories from the given directory,
avoiding self and parent.

=cut

sub get_dirs_from_directory {
    my @args = @_;
    my $parent_directory = shift @args;
    opendir DIR, $parent_directory;
    my @dirs = grep { is_directory($_, $parent_directory) }  readdir DIR;;
    closedir DIR;
    return sort {$a <=> $b} @dirs;
}

1;