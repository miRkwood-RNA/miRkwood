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

sub find_matching_count {
    my $structure = shift;
    my $stop = 1;
    my $parenthesisCounter = 0;
    my $iteration = -1;
    my $element;
    while($stop){
        $iteration += 1;
        $element = substr($structure, $iteration, 1);
        if($element eq '.'){
        }
        elsif($element eq '('){
            $parenthesisCounter += 1;
        }elsif($element eq ')'){
            $stop = 0;
        }
    }
    return $parenthesisCounter
}

sub make_loop {
    my $sequence = shift;
    my (@upper, @middle, @lower);
    my $quotient = int(length($sequence) / 2);
    my $modulo = int(length($sequence) % 2);
    push(@middle, (" ") x $quotient);
    if ($modulo != 0){
        push(@middle, substr($sequence, $quotient, 1));
    }
    push(@upper, split('', substr($sequence, 0, $quotient)));
    push(@lower, split('', substr($sequence, $quotient, $quotient)));
    my @AOA;
    $AOA[0] = [ @upper ];
    $AOA[1] = [ @middle ];
    $AOA[2] = [ @lower ];
    return [ @AOA ];
}

1;

