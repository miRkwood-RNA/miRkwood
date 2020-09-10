#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Fasta;


########################################################################
# Script to build the miRNAs fasta file for use by miRkwood.
# Given a list of species (tabulated file with the species in the first column, 
# or simply a list of species, one per line) and fasta files (mature.fa 
# and hairpin.fa) found on miRBase web site, this script extracts sequences 
# from Viridiplantae clades and then merges identical miRNAs with a new 
# header containing all initial entries for each sequence.
# The list of species can be obtained here : http://www.mirbase.org/cgi-bin/browse.pl
# Expand all Viridiplantae clades and then copy/paste the corresponding lines.
#
# Usage : ./construct_mirbase_fasta_file.pl -species <mirbase_tab_file> -fasta <fasta_file> -output <output_file>
#
# The file used by miRkwood is created with mature.fa from miRBase.
########################################################################


########## Variables
my $fasta_file = '';
my $species_file = '';
my $output_file = '';
my @filtered_fasta;
my $help;
my $help_message = "construct_mirbase_fasta_file.pl
----------
Script to build the miRNAs fasta file for use by miRkwood.
Given a list of species (tabulated file with the species in the first column, 
or simply a list of species, one per line) and fasta files (mature.fa 
and hairpin.fa) found on miRBase web site, this script extracts sequences 
from Viridiplantae clades and then merges identical miRNAs with a new 
header containing all initial entries for each sequence.
The list of species can be obtained here : http://www.mirbase.org/cgi-bin/browse.pl
Expand all Viridiplantae clades and then copy/paste the corresponding lines.

Usage : ./construct_mirbase_fasta_file.pl -species <mirbase_tab_file> -fasta <fasta_file> -output <output_file>

The file used by miRkwood is created with mature.fa from miRBase.

\n";


########## Get options
GetOptions ('species=s' => \$species_file,
            'fasta=s'   => \$fasta_file,
            'output=s'  => \$output_file,
            'help'      => \$help);


########## Validate options
if ( $help ){
    print STDERR $help_message;
    exit;
}

if ( ! -r $species_file || ! -r $fasta_file ){
    print STDERR "Missing parameter !\n" . $help_message;
    exit;
}

if ( $output_file eq '' ){
    print STDERR "You must provide an output name!\n" . $help_message;
    exit;
}


########## Main

# Filter input fasta according to species file
my $filtered_file = extract_viridiplantae_sequences( $species_file, $fasta_file);

# Merge identical entries
merge_identical_entries( $filtered_file, $output_file);

# Delete temporary file (filtered, non merged file)
unlink $filtered_file;

########## Functions

sub extract_viridiplantae_sequences{
    my ( @args ) = @_;
    my $species = shift @args;
    my $fasta = shift @args;
    my $basename = '';
    my @tab;
    my @species_list;

    if ( $fasta =~ /(.*)[.](fa|fas|fasta)$/ ){
        $basename = $1;
    }
    else{
        print STDERR "$fasta is not a fasta file. Cannot work on it.\n";
        return;
    }

    open (my $SP, '<', $species) or die "ERROR while opening $species : $!";
    open (my $OUT, '>', "${basename}_viridiplantae.fa") or die "ERROR while creating ${basename}_viridiplantae.fa : $!";

    while ( <$SP> ){
        chomp;
        @tab = split ( /\t/ );
        $tab[0] =~ s/^\s+//;
        $tab[0] =~ s/\s+$//;
        push @species_list, $tab[0];
    }

    my $seq_db = Bio::DB::Fasta->new($fasta);
    my @ids = $seq_db->get_all_primary_ids;
    
    foreach my $id( sort(@ids) ) {
        foreach my $species (@species_list) {
            if ( $seq_db->header($id) =~ /$species/ ) {
                print $OUT '>'.$seq_db->header($id)."\n";
                print $OUT $seq_db->seq($id)."\n";
            }
        }
    }

    close $SP;
    close $OUT;

    return "${basename}_viridiplantae.fa";

}

sub merge_identical_entries {
    my ( @args ) = @_;
    my $fasta = shift @args;
    my $out_file = shift @args;
    my $basename = '';
    my $name = '';
    my %seq_noms;

    if ( $fasta =~ /(.*)[.](fa|fas|fasta)$/ ){
        $basename = $1;
    }
    else{
        print STDERR "$fasta is not a fasta file. Cannot work on it.\n";
        return;
    }

    open (my $FAS, '<', $fasta) or die "ERROR while opening $fasta : $!";

    while ( <$FAS> ){
        if ( /^>(.*)$/ ){
            $name = $1;
        }
        else{
            chomp;
            if ( ! defined( $seq_noms{$_} ) ){
                $seq_noms{$_} = '';
            }
            $seq_noms{$_} .= $name.' | ';
        }
    }

    close $FAS;

    open (my $OUT, '>', $out_file) or die "ERROR when creating $out_file : $!";

    foreach my $sequence (sort keys %seq_noms) {
        print $OUT '>' . $seq_noms{$sequence} . "\n";
        print $OUT $sequence . "\n";
    }

    close $OUT;

    return;

}
