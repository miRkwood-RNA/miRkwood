package PipelineMiRNA::Utils;

# ABSTRACT: Useful, reusable methods.

use strict;
use warnings;

use feature 'switch';

=method reverse_complement

Compute de reverse complement of the sequence

Input: 
	- ARN sequence string

Output: reverse complement of the ARN sequence

=cut

sub reverse_complement {
    my @args = @_;
    my $rna  = shift @args;

    # reverse the  sequence
    my $revcomp = reverse($rna);

    # complement the reversed sequence
    $revcomp =~ tr/ACGUTacgut/UGCAAugcaa/;

    return $revcomp;
}

=method cleanup_fasta_sequence

Cleans up the FASTA sequence:
- copy-paste symbols
- lowercase

=cut

sub cleanup_fasta_sequence {
    my @args     = @_;
    my $sequence = shift @args;
    $sequence =~ s/\r//g;
    $sequence = lc($sequence) . "\n";
}

=method get_position_from_opposite_strand


=cut

sub get_position_from_opposite_strand {
    my @args   = @_;
    my $start  = shift @args;
    my $end    = shift @args;
    my $length = shift @args;
    my @res = ($length - $end - 1, $length - $start - 1,);
    return \@res ;
}

=method parse_multi_fasta

Parse a multi FASTA

Input:
  - a file handle to the multi-fasta
  - a boolean on whether to upper case the sequence shortname
Output: An hash shortname=>sequence

=cut

sub parse_multi_fasta {
    my @args         = @_;
    my ($INPUT_FH)   = shift @args;
    my $to_uppercase = 0;
    if (@args) {
        $to_uppercase = shift @args;
    }
    my %tab   = ();
    my $EMPTY = q{};
    my $nameSeq;
    my @array       = ();
    my $current_seq = $EMPTY;
    while ( my $line = <$INPUT_FH> ) {
        if ( grep { /^>/smx } $line ) {
            if ($current_seq) {
                my @res = ( $nameSeq, $current_seq );
                push @array, \@res;
            }
            chomp $line;
            $nameSeq = $line;
            $nameSeq =~ s/\s/_/xmsg;
            $nameSeq =~ s/>//xmsg;

            $current_seq = $EMPTY;
        }
        else {
            chomp $line;
            $current_seq .= $line;
        }
    }
    if ($current_seq) {
        my @res = ( $nameSeq, $current_seq );
        push @array, \@res;
    }
    return @array;
}

=method sanitize_sequence_name

=cut

sub sanitize_sequence_name {
    my @args   = @_;
    my $header = shift @args;
    my $to_uppercase = 0;
    if (@args) {
        $to_uppercase = shift @args;
    }
    my $SPACE = q{ };
    my $nameSeq = substr( $header, 0, 30 );
    $nameSeq =~ s/^\s+|\s+$//xmsg;
    $nameSeq =~ s/\|/-/xmsg;
    $nameSeq =~ s/\s/_/xmsg;
    if ($to_uppercase) {
        $nameSeq = uc $nameSeq;
    }
    return $nameSeq;
}

=method rewrite_fasta_with_TU

Rewrite a multi FASTA by changing T with U (or the reverse).

Input:
 - the way to substitute: 'T' for changing U with T, 'U' for changing T with U
 - a file handle for the input multi-fasta
 - a file handle for the output multi-fasta
Output: -

=cut

sub rewrite_fasta_with_TU {
	my (@args)      = @_;
	my $nucleotide  = shift @args;
	my ($INPUT_FH)  = shift @args;
	my ($OUTPUT_FH) = shift @args;
	while ( my $line = <$INPUT_FH> ) {
		if ( grep { /^>/msx } $line ) {
			print {$OUTPUT_FH} $line;
		}
		else {
			if ( $nucleotide eq 'U' ) {
				$line =~ s/T/U/g;
				$line =~ s/t/u/g;
			}
			elsif ( $nucleotide eq 'T' ) {
				$line =~ s/U/T/g;
				$line =~ s/u/t/g;
			}
			print {$OUTPUT_FH} $line;
		}
	}
	return;
}

=method rewrite_fasta_with_TU_in_file

Wrapper method for rewrite_fasta_with_TU

Input:
 - the way to substitute: 'T' for changing U with T, 'U' for changing T with U
 - a file path to the input multi-fasta
 - a file path to the output multi-fasta
Output: -

=cut

sub rewrite_fasta_with_TU_in_file {
	my (@args)                   = @_;
	my $nucleotide               = shift @args;
	my $sequences_to_filter_file = shift @args;
	my $output_file              = shift @args;
	open my $ENTREE_FH, '<', $sequences_to_filter_file
	  or die "Error when opening file $sequences_to_filter_file: $!";
	open my $OUTPUT_FH, '>', $output_file
	  or die "Error when opening file $output_file: $!";
	rewrite_fasta_with_TU( $nucleotide, $ENTREE_FH, $OUTPUT_FH );
	close $ENTREE_FH or die "Unable to close: $!";
	close $OUTPUT_FH or die "Unable to close: $!";
	return;
}

=method filter_fasta

Filter a multi FASTA with a given list of sequence names

Input:
 - a file handle for the FASTA to filter
 - a file handle for the output
 - the sequences to remove, as an hash name => whatever
Output: -

=cut

sub filter_fasta {
	my $FASTA_FH           = shift;
	my $RES_FH             = shift;
	my $sequence_to_filter = shift;
	my %sequence_to_filter = %{$sequence_to_filter};

	my $lineSeq;
	while ( my $line = <$FASTA_FH> ) {
		if ( grep { /^>/msx } $line ) {
			$lineSeq = substr $line, 1, -1;
		}
		if ( !exists( $sequence_to_filter{$lineSeq} ) ) {
			printf {$RES_FH} $line;
		}
	}
	return;
}

=method is_fasta_header

Checks whether the given string is an accepted FASTA header

Input:
 - a sequence string
Output:
 - True or False

=cut

sub is_fasta_header {
    my (@args) = @_;
    my $line = shift @args;
    if ( $line =~ m{^\s*>}xms ) {
        return 1;
    }
    else {
        return 0;
    }
}

=method is_fasta_line

Checks whether the given string is an accepted FASTA sequence
(ie constituted of A, T, G, C or U)

Input:
 - a sequence string
Output:
 - True or False

=cut

sub is_fasta_line {
    my (@args) = @_;
    my $line = shift @args;
    if ( $line =~ m{^\s*[ATGCUatcgu]+\s*$}xms ) {
        return 1;
    }
    else {
        return 0;
    }
}

=method is_fasta

Checks whether the given string is a FASTA sequence

Input:
 - a sequence string
Output:
 - True or False

=cut

sub is_fasta {
    my (@args) = @_;
    my $sequence = shift @args;
    my @lines = split /\n/smx, $sequence;
    foreach my $line (@lines) {
        if ( !is_fasta_header($line) && !is_fasta_line($line) ) {
            warn "Problem with line $line";
            return 0;
        }
    }
    return 1;
}

=method mask_sequence_nucleotide

Mask a typical FASTA line with 'N' symbol

Input:
 - a sequence string
Output:
 - the masked sequence

=cut

sub mask_sequence_nucleotide {
	my (@args) = @_;
	my $sequence = shift @args;
	chomp $sequence;
	$sequence =~ s/[wrkmsyd]/n/g;
	return $sequence;
}

=method find_matching_count

Get the number of matchings in a stemloop.

=cut

sub find_matching_count {
	my $structure          = shift;
	my $stop               = 1;
	my $parenthesisCounter = 0;
	my $iteration          = -1;
	my $element;
	while ($stop) {
		$iteration += 1;
		$element = substr( $structure, $iteration, 1 );
		if ($element) {
			if ( $element eq '.' ) {
			}
			elsif ( $element eq '(' ) {
				$parenthesisCounter += 1;
			}
			elsif ( $element eq ')' ) {
				$stop = 0;
			}
		}
	}
	return $parenthesisCounter;
}

=method make_loop

Return an ASCII representation of the final loop.

=cut

sub make_loop {
	my $sequence = shift;
	my $SPACE    = q{ };
	my $EMPTY    = q{ };
	my ( @top, @upper, @middle, @lower, @bottom );
	my $len      = length($sequence);
	my $quotient = int( ( $len - 2 ) / 2 );
	my $modulo   = int( $len % 2 );
	push( @middle, ($SPACE) x $quotient );
	push( @upper,  ($SPACE) x $quotient );
	push( @lower,  ($SPACE) x $quotient );

	if ( $len < 4 ) {
		push( @middle, $SPACE );
	}

	if ( $modulo != 0 ) {
		push( @middle, substr( $sequence, $quotient + 1, 1 ) );
	}
	push( @top, split( $EMPTY, substr( $sequence, 0, $quotient ) ) );
	push( @upper, split( $EMPTY, substr( $sequence, $quotient, 1 ) ) );
	push( @lower,
		split( $EMPTY, substr( $sequence, $len - $quotient - 1, 1 ) ) );
	push(
		@bottom,
		split(
			$EMPTY, reverse substr( $sequence, $len - $quotient, $quotient )
		)
	);
	my @AOA;
	$AOA[0] = [@top];
	$AOA[1] = [@upper];
	$AOA[2] = [@middle];
	$AOA[3] = [@lower];
	$AOA[4] = [@bottom];
	return @AOA;
}

=method make_ASCII_viz

Return an ASCII representation of a stemloop.

=cut

sub make_ASCII_viz {
	my ( $sequence, $structure ) = @_;

	# Variables initialisation
	my $SPACE = q{ };
	my @top;
	my @upper;
	my @middle;
	my @lower;
	my @bottom;
	my $left                = 0;
	my $right               = length($sequence) - 1;
	my $parenthesis_counter = 0;
	my $stop                = 0;
	my $parenthesis_number  = find_matching_count($structure);

	while ( !$stop ) {
		my $element_left  = substr( $structure, $left,  1 );
		my $letter_left   = substr( $sequence,  $left,  1 );
		my $element_right = substr( $structure, $right, 1 );
		my $letter_right  = substr( $sequence,  $right, 1 );

		if ( $parenthesis_counter == $parenthesis_number ) {
			$stop = 1;
		}
		elsif ( $element_left eq '.' and $element_right eq ')' ) {
			push( @top,    $letter_left );
			push( @upper,  $SPACE );
			push( @middle, $SPACE );
			push( @lower,  $SPACE );
			push( @bottom, '-' );
			$left += 1;
		}
		elsif ( $element_left eq '(' and $element_right eq ')' ) {
			push( @top,    $SPACE );
			push( @upper,  $letter_left );
			push( @middle, '|' );
			push( @lower,  $letter_right );
			push( @bottom, $SPACE );
			$left += 1;
			$right -= 1;
			$parenthesis_counter += 1;
		}
		elsif ( $element_left eq '(' and $element_right eq '.' ) {
			push( @top,    '-' );
			push( @upper,  $SPACE );
			push( @middle, $SPACE );
			push( @lower,  $SPACE );
			push( @bottom, $letter_right );
			$right -= 1;
		}
		elsif ( $element_left eq '.' and $element_right eq '.' ) {
			push( @top,    $letter_left );
			push( @upper,  $SPACE );
			push( @middle, $SPACE );
			push( @lower,  $SPACE );
			push( @bottom, $letter_right );
			$left += 1;
			$right -= 1;
		}
		else {
			$stop = 1;
		}
	}

	my $subsequence = substr( $sequence, $left, $right - $left + 1 );

	my @final_loop = make_loop($subsequence);

	#    @{$hit}[1]
	my @top_end    = @{ ${ final_loop [0] } };
	my @upper_end  = @{ ${ final_loop [1] } };
	my @middle_end = @{ ${ final_loop [2] } };
	my @lower_end  = @{ ${ final_loop [3] } };
	my @bottom_end = @{ ${ final_loop [4] } };

	my $top    = join( '', @top ) . join( '',    @top_end );
	my $upper  = join( '', @upper ) . join( '',  @upper_end );
	my $middle = join( '', @middle ) . join( '', @middle_end );
	my $lower  = join( '', @lower ) . join( '',  @lower_end );
	my $bottom = join( '', @bottom ) . join( '', @bottom_end );
	my $result = <<"END";
$top
$upper
$middle
$lower
$bottom
END
	return $result;
}

=method compute_mature_boundaries

Given a harpin and the mature position,
computes the actual span of the mature
(because of gaps).

=cut

sub compute_mature_boundaries {
	my @args = @_;
	my ( $left, $size, $top ) = @args;

	my $real_left = 0;
	while ( $real_left < $left ) {
		my $element = substr( $top, $real_left, 1 );
		if ( $element eq '-' ) {
			$left += 1;
		}
		$real_left += 1;
	}

	my $sub_string = substr( $top, $real_left, $size );
	my $gap_count = 0;
	$gap_count++ while ( $sub_string =~ m/-/g );
	my $real_size = $size + $gap_count;
	return ( $real_left - 1, $real_size + 1 );
}

=method make_hairpin_with_mature

Given a harpin, the mature position and the length,
make an hairpin with the mature standing out in HTML.

=cut

sub make_hairpin_with_mature {
	my (@args) = @_;
	my ( $hairpin, $left, $right, $length, $mode ) = @args;
	my ( $top, $upper, $middle, $lower, $bottom )  = split( /\n/, $hairpin );
	my $hairpin_with_mature;
	my $pseudo_size = $right - $left;
	my $l           = length($top);
    if (!$mode){
        $mode = 'ascii';
    }
	if ( $right >= $l ) {

		#on the other side
		my $converted_left = $length - $right + 1;
		my ( $true_left, $size ) =
		  compute_mature_boundaries( $converted_left, $pseudo_size, $bottom );
		my $bottom_mature = substr( $bottom, $true_left, $size );
		my $lower_mature = substr( $lower, $true_left, $size );
        given ($mode) {
                when (/html/) {
                    $bottom_mature = '<span class="mature">' . $bottom_mature . '</span>';
                    $lower_mature  = '<span class="mature">' . $lower_mature . '</span>';
                }
                when (/ascii/) {
                    $bottom_mature = uc $bottom_mature;
                    $lower_mature  = uc $lower_mature;
                }
        }
		substr( $bottom, $true_left, $size ) = $bottom_mature;
		substr( $lower, $true_left, $size )  = $lower_mature;
	}
	else {
		my ( $true_left, $size ) =
		  compute_mature_boundaries( $left, $pseudo_size, $top );
	    my $top_mature = substr( $top, $true_left, $size );
		my $upper_mature = substr( $upper, $true_left, $size );
        given ($mode) {
                when (/html/) {
                    $top_mature   = '<span class="mature">' . $top_mature . '</span>';
                    $upper_mature = '<span class="mature">' . $upper_mature . '</span>';
                }
                when (/ascii/) {
                    $top_mature   = uc $top_mature;
                    $upper_mature = uc $upper_mature;
                }
        }
		substr( $top, $true_left, $size) = $top_mature;
		substr( $upper, $true_left, $size ) = $upper_mature;

	}
	$hairpin_with_mature = <<"END";
$top
$upper
$middle
$lower
$bottom
END
	return $hairpin_with_mature;
}

=method filter_mirbase_hairpins

Filter a multi-FASTA file to match the contents of another
multi-FASTA file, in the MiRbase hairpin/precursor way.

Usage: filter_mirbase_hairpins('data/MirbaseFile.txt',
                               'data/hairpin.fa',
                               'hairpin-filtered.fa');
=cut

sub filter_mirbase_hairpins {
	my (@args)                        = @_;
	my $sequences_to_filter_file      = shift @args;
	my $sequences_to_be_filtered_file = shift @args;
	my $output                        = shift @args;
	open my $ENTREE_FH, '<', $sequences_to_filter_file
	  or die "Error when opening sequences -$sequences_to_filter_file-: $!";
	my %sequences_to_filter =
	  PipelineMiRNA::Utils::parse_multi_fasta( $ENTREE_FH, 1 );
	close $ENTREE_FH or die "Unable to close: $!";

	open( my $FSeq, '<', $sequences_to_be_filtered_file )
	  or die "Error when opening file $sequences_to_be_filtered_file: $!";
	open( my $RES_FH, '>>', $output )
	  or die "Error when opening file $output: $!";
	my $lineSeq;
	my $nameSeq;
	while ( my $line = <$FSeq> ) {
		if ( grep { /^>/msx } $line ) {
			$nameSeq = ( split( q{ }, $line ) )[0];
			$nameSeq = uc $nameSeq;
		}
		if ( exists( $sequences_to_filter{ uc $nameSeq } ) ) {
			printf {$RES_FH} $line;
		}
	}
	close $FSeq   or die "Unable to close: $!";
	close $RES_FH or die "Unable to close: $!";
	return;
}

=method get_element_of_split

Get the n-th element of the split on a a given string using the given separator

Usage: get_element_of_split($b, '-', 0);
=cut

sub get_element_of_split {
	my @args  = @_;
	my $value = shift @args;
	my $sep   = shift @args;
	my $rank  = shift @args;
	my @split = split( /$sep/, $value );
	return $split[$rank];
}

1;

=method restrict_num_decimal_digits
 
restrict the number of digits after the decimal point

Usage:
    my $num = PipelineMiRNA::Utils::restrict_num_decimal_digits($number, 3);

=cut

sub restrict_num_decimal_digits {
	my $num         = shift;    #the number to work on
	my $digs_to_cut = shift;    # the number of digits after
	                            # the decimal point to cut
	                            #(eg: $digs_to_cut=3 will leave
	                            # two digits after the decimal point)

	if ( $num =~ /\d+\.(\d){$digs_to_cut,}/ ) {

		# there are $digs_to_cut or
		# more digits after the decimal point
		$num = sprintf( "%." . ( $digs_to_cut - 1 ) . "f", $num );
	}
	return $num;
}

=method compute_mfei_and_amfe

Compute the MFEI and the AMFE of a given sequence

=cut

sub compute_mfei_and_amfe {
    my @args     = @_;
    my $sequence = shift @args;
    my $energy   = shift @args;

    my $length   = length $sequence;
    my $gc_content = compute_gc_content($sequence);

    my $amfe = compute_amfe($sequence, $energy);
    my $mfei = 0;
    if ($gc_content) {
        $mfei = $amfe / $gc_content;
    }

    return ($mfei, $amfe);
}

=method compute_amfe

Compute the AMFE of a given sequence

=cut

sub compute_amfe {
    my @args     = @_;
    my $sequence = shift @args;
    my $energy   = shift @args;
    my $length   = length $sequence;
    my $amfe = ( $energy / $length ) * 100;
    return $amfe;
}

=method compute_gc_content

Compute the GC content of a given sequence
Returns the percentage (float between 0 and 100).
=cut

sub compute_gc_content {
    my @args     = @_;
    my $sequence = shift @args;

    my @dna      = split( //, $sequence );
    my $length   = scalar @dna;
    my $gc_count = 0;

    for my $nucl ( @dna ) {
        if (
            $nucl =~ m{
              [cg]     # G or C
            }smxi
          )
        {
            $gc_count++;
        }
    }
    return  $gc_count / $length * 100;
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

=method slurp_file

Return the contents of a file

=cut

sub slurp_file {
    my @args = @_;
    my $file = shift @args;
    open my $fh, '<', $file or die $!;
    my $contents = do { local $/; <$fh> };
    close $fh or die $!;
    return $contents;
}
