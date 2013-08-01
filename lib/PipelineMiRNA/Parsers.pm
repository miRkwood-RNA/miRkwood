package PipelineMiRNA::Parsers;

use strict;
use warnings;


=method parse_pvalue

Parse the contents of the pvalue file the results structure of a given job identifier

=cut

sub parse_pvalue {
    my @args   = @_;
    my $pvalue = shift @args;
    my $result;
    open( my $FH, '<', $pvalue ) or die "Error when opening file $pvalue: $!";
    while ( my $line = <$FH> ) {
        if ( $line =~ /(.*)\t(.*)\t(.*)/xms ) {
            $result = $3;
        }
    }
    close $FH or die("Error when closing: $!");
    return $result;
}

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

sub parse_selfcontain {
    my @args            = @_;
    my $selfcontain_out = shift @args;
    my $result;
    open( my $FH, '<', $selfcontain_out )
      or die "Error when opening file: $!";
    while ( my $line = <$FH> ) {
        if ( $line =~ /(.*) (.*)/ ) {
            $result = $2;
        }
    }
    close $FH or die("Error when closing: $!");
    return $result;
}

sub parse_vienna {
    my @args       = @_;
    my $vienna_out = shift @args;
    my @res;
    open( my $FH, '<', $vienna_out )
      or die "Error when opening file: $!";
    while ( my $line = <$FH> ) {
        if ( $line =~ /(.*)\t(.*)\t(.*)/xms ) {
            $res[0] = $2;    # récupération sequence
            $res[1] = $3;    # récupération Vienna
        }
    }
    close $FH or die("Error when closing: $!");
    return @res;
}

sub parse_alignment {
    my @args            = @_;
    my $alignement_file = shift @args;
    open( my $FH, '<', $alignement_file )
      or die "Error when opening file: $!";
    my $align = 'none';
    while ( my $line = <$FH> ) {
        if ( $line =~ /^C4/xms ) {
            $align = $alignement_file;
            last;
        }
    }
    close $FH or die("Error when closing: $!");
    return $align;
}

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

1;
