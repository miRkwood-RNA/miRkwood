package PipelineMiRNA::Maskers;

# ABSTRACT: Masking sequences

use strict;
use warnings;

use File::Spec;
use Log::Message::Simple qw[msg error debug];

use PipelineMiRNA::Programs;
use PipelineMiRNA::Parsers;

=method get_coding_region_masking_information

Retrieve coding region masking information 
by running Blastx and parsing the output

=cut

sub get_coding_region_masking_information {
    my ( $dirData, $dirJob, $plant ) = @_;

    my $uploaded_sequences =
      File::Spec->catfile( $dirJob, 'input_sequences.fas' );
    my $input_sequences = File::Spec->catfile( $dirJob, 'Sequences.fas' );

    my $blast_database = File::Spec->catfile( $dirData, "$plant.fas" );
    my $blast_output   = File::Spec->catfile( $dirJob,  'outBlast.txt' );
    my $blastx_options = '-outfmt 6 -max_target_seqs 1 -evalue 1E-5';
    PipelineMiRNA::Programs::run_blast(
                                        $uploaded_sequences, $blast_database,
                                        $blastx_options,     $blast_output
    ) or die('Problem when running Blastx');

    my %blast_seqs = PipelineMiRNA::Parsers::parse_blast_output($blast_output);
    return %blast_seqs;
}

=method get_trna_masking_information

Retrieve tRNA masking information by running tRNAscan-SE and parsing the output

=cut

sub get_trna_masking_information {
    my ( $dirJob ) = @_;

    my $uploaded_sequences =
      File::Spec->catfile( $dirJob, 'input_sequences.fas' );
    my $output = File::Spec->catfile( $dirJob, 'trnascanse.out' );

    PipelineMiRNA::Programs::run_tRNAscanSE_on_file($uploaded_sequences, $output
    ) or die('Problem when running tRNAscanSE');

    my %trna_seqs = PipelineMiRNA::Parsers::parse_tRNAscanSE_output($output);
    return %trna_seqs;
}

=method get_rnammer_masking_information

Retrieve ribosomal RNA masking information by running RNAmmer and parsing the output

=cut

sub get_rnammer_masking_information {
    my ( $dirJob ) = @_;

    my $uploaded_sequences =
      File::Spec->catfile( $dirJob, 'input_sequences.fas' );
    my $output = File::Spec->catfile( $dirJob, 'rnammer.out' );
    my $kingdom = 'euk';
    PipelineMiRNA::Programs::run_rnammer_on_file( $uploaded_sequences, $kingdom, $output )
      or die('Problem when running RNAmmer');

    my %rnammer_seqs = PipelineMiRNA::Parsers::parse_rnammer_output($output);
    return %rnammer_seqs;
}

=method mask_sequences


 Usage : @sequences = mask_sequences(\%filter, @sequences_array);

=cut

sub mask_sequences {
    my @args = @_;
    my $filter = shift @args;
    my %filter = %{$filter};
    my @sequences_array = @args;
    my @sequences_results;
    foreach my $item (@sequences_array){
        my ($name, $sequence) = @{$item};
        my @split = split(/_/, $name);
        my $id = $split[0];
        if (exists $filter{$id}){
            foreach my $positions(@{$filter{$id}}){
                my %pos = %{$positions};
                debug( "Masking sequence $name between <$pos{'start'}> and <$pos{'end'}>", PipelineMiRNA->DEBUG() );
                $sequence = mask_sequence($sequence, $pos{'start'}, $pos{'end'});
            }
        }
        my @res = ( $name, $sequence );
        push @sequences_results, \@res;
    }
    return @sequences_results;
}

=method mask_sequence

 Usage: my $masked = mask_sequence($seq, 3, 6);

=cut

sub mask_sequence {
    my @args = @_;
    my ( $sequence, $start, $end ) = @args;
    if ($end < $start){
       ( $start, $end ) = ( $end, $start );
    }
    my $length = $end - $start + 1;
    my $mask = ( 'N' x $length );
    substr($sequence, $start, $length) = $mask;
    return $sequence;
}

1;