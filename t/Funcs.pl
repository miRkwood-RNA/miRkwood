use strict;
use warnings;

sub input_file {
    my @args = @_;
    return File::Spec->catfile( $FindBin::Bin, 'data', @args );
}

sub slurp_file {
    my @args = @_;
    my $file = shift @args;
    open my $fh, '<', $file or die $!;
    my $contents = do { local $/; <$fh> };
    close $fh or die $!;
    return $contents;
}

1;
