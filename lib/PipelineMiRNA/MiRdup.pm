package PipelineMiRNA::MiRdup;

use strict;
use warnings;

use File::Spec;
use File::Basename;
use PipelineMiRNA::Programs;

=method make_prediction_source_file

Given a set of sequences, write on disk a file appropriately
formatted for MiRdup prediction

=cut

sub make_prediction_source_file {
    my ($self, @args) = @_;
    my $output_file = shift @args;
    my $sequences = shift @args;
    my %sequences = %{$sequences};
    open( my $FOut, '>', $output_file )
      or die "Error when opening $output_file: $!";
    while ( my ($name, $seq) = each %sequences ){
          print {$FOut} "$name\t$seq\n";
    }
    close $FOut or die "Error when closing $output_file: $!";
    return (-e $output_file);
}

=method parse_predictions_file

Parse a miRdup prediction result file and return the
predictions as an hash <name> => <prediction>

Usage: 
    my %results = PipelineMiRNA::MiRdup->parse_predictions_file($res_file);

=cut

sub parse_predictions_file {
    my ($self, @args) = @_;
    my $prediction_file = shift @args;
    my %results = ();
    open( my $FH, '<', $prediction_file ) or die "Error when opening file $prediction_file: $!";
    while ( my $line = <$FH> ) {
        my ($name, $precursor, $structure, $predict5, $predict3) = split('\t', $line);
        $predict5 =~ m{
            ^([acgtuACGTU]*)
            (.*)$
            }xms;
        $results{$name} = $1;
    }
    close $FH or die("Error when closing: $!");
    return %results;
}

=method predict_with_mirdup

Given a set of sequences, predict their miRNAs
using miRdup in the given path.

Usage:
    my $res_file = PipelineMiRNA::MiRdup->predict_with_mirdup('prediction_results.txt', \%sequences);

=cut

sub predict_with_mirdup {
    my ($self, @args) = @_;
    my $prediction_file = shift @args;
    my %sequences = %{shift @args};
    my ($volume, $origin_directories, $origin_file) = File::Spec->splitpath($prediction_file);
    $self->make_prediction_source_file($prediction_file, \%sequences)
        or die "Error when making prediction_source_file: $!";
    my $result_file = PipelineMiRNA::Programs::run_mirdup_prediction_on_sequence_file($prediction_file);
    return $result_file;
}

1;