package miRkwood::Parsers;

# ABSTRACT: Parsing methods

use strict;
use warnings;

use miRkwood::Utils;

=method parse_pvalue

Parse the contents of the given pvalue file

=cut

sub parse_pvalue {
    my @args   = @_;
    my $pvalue = shift @args;
    my $result;
    open( my $FH, '<', $pvalue ) or die "Error when opening file $pvalue: $!";
    while ( my $line = <$FH> ) {

        if ( $line =~ m{
           ^                #Begin of line
           (.*?)            #A non-greedy sequence of characters, captured
           \t               #One tabulation
           (.*?)            #A non-greedy sequence of characters, captured
           \t               #One tabulation
           (.*?)            #A non-greedy sequence of characters, captured
           $                #End of line
           }smx ) {
           $result = $3;
       }

    }
    close $FH or die("Error when closing: $!");
    return $result;
}

=method parse_mfei

Parse the contents of the given MFEI file

=cut

sub parse_mfei {
    my @args     = @_;
    my $mfei_out = shift @args;
    my @res;
    open( my $FH, '<', $mfei_out )
      or die "Error when opening file: $!";
    while ( my $line = <$FH> ) {
        if ( $line =~ /(.*)\t(.*)\t(.*)\t(.*)/xms ) {
            push @res, $2;
            push @res, $3;
            push @res, $4;
        }
    }
    close $FH or die("Error when closing: $!");
    return @res;
}

=method parse_selfcontain

Parse the contents of the given SelfContain file

=cut

sub parse_selfcontain {
    my @args            = @_;
    my $selfcontain_out = shift @args;
    my $result;
    open( my $FH, '<', $selfcontain_out )
      or die "Error when opening file: $!";
    while ( my $line = <$FH> ) {
        if ( $line =~ m{
           ^                # Begin of line
           (.*?)            # A non-greedy sequence of characters, captured
           \s+?             # Some whitespace
           (.*?)            # A non-greedy sequence of characters, captured
           $                # End of line
           }xms ) {
           $result = $2;
       }
    }
    close $FH or die("Error when closing: $!");
    return $result;
}

=method parse_vienna

Parse the contents of the given Vienna file

=cut

sub parse_vienna {
    my @args       = @_;
    my $vienna_out = shift @args;
    my @res;
    open( my $FH, '<', $vienna_out )
      or die "Error when opening file: $!";
    while ( my $line = <$FH> ){
        if ( $line =~ m{
           ^\s*?            # Begin of line and possibly some whitespace
           (.*?)            # Whatever, non-greedy, captured
           \s+?             # Some whitespace
           ([aAcCgGtTuU]*)  # A sequence of nucleotides, captured
           \s+?             # Some whitespace
           ([\(\.\)]+)      # A sequence of either '.', '(' or ')'
           \s+?             # Some whitespace
           (.*?)            # Whatever, non-greedy, captured
           \s*$             # Some whitespace until the end
           }xms ) {
            $res[0] = $2;    # récupération sequence
            $res[1] = $3;    # récupération Vienna
       }
    }
    close $FH or die("Error when closing: $!");
    return @res;
}

=method parse_Vienna_line

Parse the given Vienna format line, return a couple (structure, energy)

=cut

sub parse_Vienna_line {
    my @args = @_;
    my $line = shift @args;
    my ( $structure, $energy ) =
                $line =~ m{
                       ^                #Begin of line
                       ([\.()]+?)       #A sequence of ( ) .
                       \s+?             #Some whitespace
                       \(\s*?           #Opening parenthesis and maybe whistespace
                           ([-.\d]*?)   #
                       \s*?\)           #Closing parenthesis and maybe whistespace
                   }smx;
    return ( $structure, $energy );
}

=method parse_RNAfold_output

Parse the output of RNAfold,
Return a string "<sequence name> <sequence> <vienna>"

=cut

sub parse_RNAfold_output {
    my (@args)      = @_;
    my ($file)  = shift @args;
    my ( $nameSeq, $dna, $structure, $energy );
    my $result = qw{};
    open( my $INPUT_FH, '<', $file ) or die "Error when opening file $file: $!";
    while ( my $line = <$INPUT_FH> ) {
        if ( $line =~ m{
           ^>               # Begin of line and a caret symbol
           (.*?)            # Whatever, non-greedy, captured
           \s*?$            # Possibly some whitespace, until the end
        }xms ){
            $nameSeq = $1;
        }
        elsif ( miRkwood::Utils::is_fasta_line_relaxed($line) ){
            $dna = $line;
            $dna =~ s/^\s+//;
            $dna =~ s/\s+$//;
        }
        elsif ( $line =~ m{^\s*?\(([-.\d ]*?)\)} ){

        }
        else{
            ( $structure, $energy ) = parse_Vienna_line($line);
        }
    } # while
    close $INPUT_FH or die "Error when closing file $file: $!";
    return ( $nameSeq, $dna, $structure, $energy );
}

=method parse_alternative_candidates_file

Parse the alternative candidates file

=cut

sub parse_alternative_candidates_file {
    my (@args)      = @_;
    my ($file)  = shift @args;
    my %results;
    open( my $INPUT_FH, '<', $file ) or die "Error when opening file $file: $!";
    while ( my $line = <$INPUT_FH> ) {
        chomp $line;
        my ($name, $sequence, $structure, $mfei) = split(/\t/, $line);
        my %result;
        $result{'sequence'}  = $sequence;
        $result{'structure'} = $structure;
        $result{'mfei'} = $mfei;
        $results{(substr $name, 1)} = \%result;
    } # while
    close $INPUT_FH or die "Error when closing file $file: $!";
    return %results;
}

=method index_blast_output

Looping in blast_output, indexing sequences found

=cut

sub index_blast_output {
    my @args         = @_;
    my $blast_output = shift @args;

    my %blast_seqs;
    open( my $FOut, '<', $blast_output )
      || die "Error opening $blast_output: $!";
    while ( my $line = <$FOut> ) {
        my @name = split( /\t/, $line );
        $blast_seqs{ $name[0] } = 1;
    }
    close $FOut
      || die "Error closing $blast_output: $!";

    return \%blast_seqs;
}

=method parse_blast_output

Looping in blast_output, indexing (sequences, positions) found

=cut

sub parse_blast_output {
    my @args = @_;
    my $file = shift @args;
    my %results;
    open( my $INPUT_FH, '<', $file ) or die "Error when opening file $file: $!";
    while ( my $line = <$INPUT_FH> ) {
        my @values = split( /\t/, $line );
        my $name  = $values[0];
        my $struct = {
            'start' => $values[6],
            'end'   => $values[7]
        };
        %results = push_to_array_in_hash(\%results, $name, $struct);
    }
    close $INPUT_FH or die "Error when closing file $file: $!";
    return %results;
}

=method parse_tRNAscanSE_output

Parse the contents of a tRNAscanSE output file

=cut

sub parse_tRNAscanSE_output {
    my @args         = @_;
    my ($file)  = shift @args;
    my %results;
    open( my $INPUT_FH, '<', $file ) or die "Error when opening file $file: $!";
    while ( my $line = <$INPUT_FH> ) {
        my @values = split( /[\t]/, $line );
        my $struct = {
            'start' => int(trim($values[2])),
            'end'   => int(trim($values[3]))
        };
        %results = push_to_array_in_hash(\%results, trim($values[0]), $struct);
    }
    close $INPUT_FH or die "Error when closing file $file: $!";
    return %results;
}

=method parse_rnammer_output

Parse the contents of a RNAmmer output file in GFF

=cut

sub parse_rnammer_output {
    my @args = @_;
    my ($file) = shift @args;
    my %results;
    open( my $INPUT_FH, '<', $file ) or die "Error when opening file $file: $!";
    while ( my $line = <$INPUT_FH> ) {
        if ( $line =~ /^#/ ) {
            next;
        }
        my @values = split( /[\t]/, $line );
        my $struct = {
            'start' => int( trim( $values[3] ) ),
            'end'   => int( trim( $values[4] ) )
        };
        %results =
          push_to_array_in_hash( \%results, trim( $values[0] ), $struct );
    }
    close $INPUT_FH or die "Error when closing file $file: $!";
    return %results;
}

=method trim

Trim a string

=cut

sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

=method push_to_array_in_hash

Helper method to push a value in a hash of arrays

=cut

sub push_to_array_in_hash {
    my @args = @_;
    my ($hash, $key, $value) = @args;
    my %hash = %{$hash};
    if (! $hash{$key}) {
        my @tab = ($value);
        $hash{$key} = \@tab;
    } else {
        my @values = @{$hash{$key}};
        push @values, $value;
        $hash{$key} = \@values;
    }
    return %hash;
}


=method parse_custom_exonerate_output

Parse our custom Exonerate output
(as defined in Programs::run_exonerate)
This method can also be used to parse
RNAcomp output.

Usage: parse_exonerate_alignment($alignment);

=cut

sub parse_custom_exonerate_output{
    my @args = @_;
    my $yaml_file = shift @args;

    my $yaml = miRkwood::get_yaml_file($yaml_file) or die("Error when parsing YAML file $yaml_file");
    my @contents = @{$yaml} or die("Error when parsing YAML file $yaml_file");

    my %results;
    foreach my $element (@contents){
        my ($begin_target, $end_target) = ($element->{'begin_target'} + 1, $element->{'end_target'});
        my ($begin_query, $end_query) = ($element->{'begin_query'}, $element->{'end_query'});
        if ($begin_query < $end_query){
            $begin_query += 1;
        }else{
            $end_query += 1;
        }
        my $key = "$begin_target-$end_target";
         my $value = {
             'name' => $element->{'name'},
             'seq' => $element->{'seq'},
             'score' => $element->{'score'},
             'begin_target' => $begin_target,
             'end_target' => $end_target,
             'strand_target' => $element->{'strand_target'},
             'begin_query' => $begin_query,
             'end_query' => $end_query,
             'def_query' => $element->{'def_query'},
             'seq' => $element->{'seq_query'},
             'score' => $element->{'score'},
             'alignment' => parse_exonerate_alignment($element->{'alignment'}),
         };
         if (! exists $results{$key}) {
             $results{$key} = ();
         }
         push @{$results{$key}}, $value;
    }
    return %results;
}

=method parse_exonerate_alignment

Parse the alignment given in the output of Exonerate

Usage: parse_exonerate_alignment($alignment);

=cut

sub parse_exonerate_alignment {
    my @args = @_;
    my $alignment = shift @args;

    my $SPACE = q{ };
    my @top;
    my @middle;
    my @bottom;

    for (split /\n/mxs, $alignment ) {
        $_ =~ s/T/U/g;
        $_ =~ s/t/u/g;
        my ($first, $second, $label) = split($SPACE, $_);
        if ($label ne 'none') {
            push( @top,    uc $second );
            push( @bottom, uc $first );
            if ( uc $first eq uc $second ) {
                push( @middle, '|' );
            }else{
                push( @middle, $SPACE );
            }
        }
    }
    my $top = join('', @top);
    my $middle = join('', @middle);
    my $bottom= join('', @bottom);

    my $result = <<"END";
$top
$middle
$bottom
END
    return $result;
}


1;
