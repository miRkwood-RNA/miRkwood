package PipelineMiRNA::Utils;

use strict;
use warnings;

sub parse_multi_fasta {
    my ( $INPUT_FH ) = @_ ;
    my %tab      = ();
    my $nameSeq;
    while ( my $line = <$INPUT_FH> ) {
        if ( grep { /^>/ } $line ) {
            chomp $line;
            $nameSeq = (split(' ', $line))[0];
            $tab{$nameSeq} = '';
        }
        else {
            chomp $line;
            $tab{$nameSeq} = $tab{$nameSeq}.$line;
        }
    }
    return %tab;
}

1;

