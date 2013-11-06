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

=method make_validation_source_file

Given a set of sequences, write on disk a file appropriately
formatted for MiRdup validation

=cut

sub make_validation_source_file {
    my ($self, @args) = @_;
    my $name_base = shift @args;
    my $mature_seq = shift @args;
    my $structure = shift @args;
    my @alignments = @args;

    my @results;
    foreach my $position (@alignments) {
        my ($left, $right) = split(/-/, $position);;
        my $name = $name_base . '__' . $position;
        my $precursor_seq = substr($mature_seq, $left, $right - $left);
        my $res = $name . "\t" .
                  $precursor_seq . "\t" .
                  $mature_seq . "\t" .
                  $structure;
        push @results, $res;
    }
    return join("\n", @results);
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
    return $self->parse_predictions_file($result_file);
}


=method predict_sequence_with_mirdup

Given a sequence, predict its miRNAs
using miRdup in the given path.

Usage:
    my $res_file = PipelineMiRNA::MiRdup->predict_sequence_with_mirdup('prediction_results.txt', $sequence);

=cut

sub predict_sequence_with_mirdup {
    my ($self, @args) = @_;
    my $prediction_path = shift @args;
    my $sequence = shift @args;
    my $result_file = PipelineMiRNA::Programs::run_mirdup_prediction_on_sequence($prediction_path, $sequence);
    return $self->parse_predictions_file($result_file);
}

=method validate_with_mirdup

Validate using miRdup

Usage:
    my $res_file = PipelineMiRNA::MiRdup->predict_with_mirdup('prediction_results.txt', \%sequences);

=cut

sub validate_with_mirdup {
    my ($self, @args) = @_;
    my $output_file = shift @args;
    my $name_base = shift @args;
    my $mature_seq = shift @args;
    my $structure = shift @args;
    my @alignments = @args;

    open( my $FOut, '>', $output_file )
      or die "Error when opening $output_file: $!";
    print {$FOut} $self->make_validation_source_file($name_base, $mature_seq,$structure, @alignments);
    close $FOut or die "Error when closing $output_file: $!";
    my $result_file = PipelineMiRNA::Programs::run_mirdup_validation_on_file($output_file);
    return $self->parse_validation_output($result_file);
}

=method parse_validation_output

Parse a miRdup validation result file and return the
validation

Usage:
    my %results = PipelineMiRNA::MiRdup->parse_validation_output($res_file);

=cut

sub parse_validation_output {
    my ($self, @args) = @_;
    my $validation_file = shift @args;
    my %results = ();
    my %boolean_mapping = ('true'=> 1, 'false' => 0);
    open( my $FH, '<', $validation_file ) or die "Error when opening file $validation_file: $!";
    while ( my $line = <$FH> ) {
        my ($name, $mature, $precursor, $structure, $validation, $confidence) = split('\t', $line);
        $validation = $boolean_mapping{$validation};
        $results{$name} = $validation;
    }
    close $FH or die("Error when closing: $!");
    return %results;
}

1;