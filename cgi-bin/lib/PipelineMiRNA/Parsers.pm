package PipelineMiRNA::Parsers;

# ABSTRACT: Parsing methods

use strict;
use warnings;


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
                       \s*$             #Whitespace until the end
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
        elsif ( $line =~ m{
           ^\s*?            # Begin of line and possibly some whitespace
           ([aAcCgGtTuU]*)  # A sequence of nucleotides, captured
           \s*?$            # Possibly some whitespace, until the end
        }xms ){
            # récupération de la sequence adn
            $dna = $1;
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

sub parse_tRNAscanSE_output {
    my @args         = @_;
    my ($file)  = shift @args;
    my @results;
    open( my $INPUT_FH, '<', $file ) or die "Error when opening file $file: $!";
    while ( my $line = <$INPUT_FH> ) {
        my @values = split( /[\t]/, $line );
        my @trna = ( trim($values[0]), int(trim($values[2])), int(trim($values[3])) );
        push @results, \@trna;
    }
    close $INPUT_FH or die "Error when closing file $file: $!";
    return @results;
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

1;
