package PipelineMiRNA::Utils;

use strict;
use warnings;

sub parse_multi_fasta {
    my ($INPUT_FH) = @_;
    my %tab = ();
    my $nameSeq;
    while ( my $line = <$INPUT_FH> ) {
        if ( grep { /^>/ } $line ) {
            chomp $line;
            $nameSeq = ( split( ' ', $line ) )[0];
            $nameSeq =~ s/\|/-/g;
            $tab{$nameSeq} = '';
        }
        else {
            chomp $line;
            $tab{$nameSeq} = $tab{$nameSeq} . $line;
        }
    }
    return %tab;
}

sub filter_fasta {
    my $FASTA_FH = shift;
    my $RES_FH = shift;
    my $sequence_to_filter = shift;
    my %sequence_to_filter = %$sequence_to_filter;

    my $lineSeq;
    my $RESULT = "";
    while ( my $line = <$FASTA_FH> ) {
        if ( grep { /^>/msx } $line ) {
            $lineSeq = substr $line, 1, -1;
        }
        if(!exists($sequence_to_filter{$lineSeq})) {
            printf $RES_FH $line;
        }
    }
    return;
}

1;