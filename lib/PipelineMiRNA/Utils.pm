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

=method make_loop

Return an ASCII representation of the final loop.

=cut

sub make_loop {
    my $sequence = shift;
    my (@top, @upper, @middle, @lower, @bottom);
    my $quotient = int( (length($sequence) - 2) / 2);
    print $sequence, " ", length($sequence), " ", $quotient;
    my $modulo = int(length($sequence) % 2);
    push(@middle, (" ") x $quotient);
    push(@upper, (" ") x $quotient);
    push(@lower, (" ") x $quotient);
    if ($modulo != 0){
        push(@middle, substr($sequence, $quotient + 1, 1));
    }
    push(@top, split('', substr($sequence, 0, $quotient)));
    push(@upper, split('', substr($sequence, $quotient, 1)));
    push(@lower, split('', substr($sequence, $quotient + 2, 1)));
    push(@bottom, split('', substr($sequence, $quotient + 2, $quotient)));
    my @AOA;
    $AOA[0] = [ @top ];
    $AOA[1] = [ @upper ];
    $AOA[2] = [ @middle ];
    $AOA[3] = [ @lower ];
    $AOA[4] = [ @bottom ];
    print "\n============\n";
    print @upper;
    print "\n============\n";
    return @AOA;
}


=method make_ASCII_viz

Return an ASCII representation of a stemloop.

=cut

sub make_ASCII_viz {
    print "make_ASCII_viz()\n";
    my ($sequence, $structure) = @_;

    # Variables initialisation
    my @top;
    my @upper;
    my @middle;
    my @lower;
    my @bottom;
    my $left = 0;
    my $right = length($sequence) - 1;
    my $parenthesis_counter = 0;
    my $stop = 0;
    my $parenthesis_number = find_matching_count($structure);

    while (! $stop){
        my $element_left = substr($structure, $left, 1);
        my $letter_left = substr($sequence, $left, 1);
        my $element_right = substr( $structure, $right, 1);
        my $letter_right = substr( $sequence, $right, 1);

        if ($element_left eq '.' and $element_right eq ')'){
            push(@top, $letter_left);
            push(@upper, ' ');
            push(@middle, ' ');
            push(@lower, ' ');
            push(@bottom, '-');
            $left += 1;
        }elsif( $element_left eq "(" and $element_right eq ')'){
            push(@top, ' ');
            push(@upper, $letter_left);
            push(@middle, '|');
            push(@lower, $letter_right);
            push(@bottom, ' ');
            $left += 1;
            $right -= 1;
            $parenthesis_counter += 1;
        }elsif( $element_left eq "(" and $element_right eq '.'){
            push(@top, '-');
            push(@upper, ' ');
            push(@middle, ' ');
            push(@lower, ' ');
            push(@bottom, $letter_right);
            $right -= 1;
        }elsif( $element_left eq "." and $element_right eq '.'){
            push(@top, $letter_left);
            push(@upper, ' ');
            push(@middle, ' ');
            push(@lower, ' ');
            push(@bottom, $letter_right);
            $left += 1;
            $right -= 1;
        }else{
            $stop = 1;
        }

        if ($left ge $parenthesis_number
            and length($sequence) - $right - 1 ge $parenthesis_number){
            $stop = 1;
        }

    }

    my $subsequence = substr($sequence, $left, $right - $left - 1);
    my @final_loop = make_loop($subsequence);
#    @{$hit}[1]
    my @top_end =  @{${final_loop[0]}};
    my @upper_end =  @{${final_loop[1]}};
    my @middle_end =  @{${final_loop[2]}};
    my @lower_end =  @{${final_loop[3]}};
    my @bottom_end =  @{${final_loop[4]}};

    my $result = <<END;
@top @top_end
@upper @upper_end
@middle @middle_end
@lower @lower_end
@bottom @bottom_end
END
    return $result;
}

1;

