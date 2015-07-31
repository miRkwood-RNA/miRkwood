package miRkwood::ResultsExporter::ResultsExporter;

# ABSTRACT: Abstract class for exporting results.

use strict;
use warnings;

use File::Spec;

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my $mirna_type = shift @args;
    my $self = {
        identifier => undef,
        results => undef,
        sequences_to_export => undef,
        mirna_type => $mirna_type,
    };
    bless $self, $class;
    return $self;
}

sub initialize {
    my ( $self, @args ) = @_;
    my ($identifier, $results, $sequences_to_export) = @args;
    $self->{'identifier'} = $identifier;
    $self->{'results'} = $results;
    $self->{'sequences_to_export'} = $sequences_to_export;
    return $self;
}

sub get_identifier {
    my ( $self, @args ) = @_;
    return $self->{'identifier'};
}

sub get_content_disposition {
    return 'attachment';
}

sub get_file_extension {
    die ('Unimplemented method get_file_extension');
}

sub get_type {
    my ( $self, @args ) = @_;
    if ( $self->{'mirna_type'} eq 'New' ) {
        return '_New';
    }
    elsif ( $self->{'mirna_type'} eq 'Known' ) {
        return '_Known';
    }
    else {
        return '';
    }
}

sub get_filename {
    my ( $self, @args ) = @_;
    if ( $self->get_identifier() ne '' ){
        return 'Results-' . $self->get_identifier() . $self->get_type() . '.' . $self->get_file_extension();
    }
    else {
        return 'Results' . $self->get_identifier() . $self->get_type() . '.' . $self->get_file_extension();
    }
}


sub get_header {
    return q{};
}

=method get_content_type

Abstract method

=cut

sub get_content_type {
    die ('Unimplemented method get_content_type');
}

=method export_candidate

Abstract method

=cut

sub export_candidate {
    die ('Unimplemented method export_candidate');
}

sub get_sequences_to_export {
    my ( $self, @args ) = @_;
    my $sequences_to_export_ref = $self->{'sequences_to_export'};
    my @sequences_to_export;
    if ( !eval { @sequences_to_export = @{$sequences_to_export_ref} } ) {
        @sequences_to_export = ();
    }
    return \@sequences_to_export;
}

sub is_no_sequence_selected {
    my ( $self, @args ) = @_;
    return !( scalar @{$self->get_sequences_to_export()} );
}

=method get_sorted_keys

Order the keys by name and starting position

=cut

sub get_sorted_keys {
    my ( $self, @args ) = @_;
    my %results = %{$self->{'results'}};
    my @keys = sort {
        ( $results{$a}->{'name'} cmp $results{$b}->{'name'} )
          || (
            $results{$a}->{'start_position'} <=> $results{$b}->{'start_position'} )
    } keys %results;
    return @keys;
}

=method is_sequence_to_export

Predicate on whether the given sequence has to
be exported, ie whether it figures in the
sequence to export array, or if no sequence
has been selected.

=cut

sub is_sequence_to_export {
    my ( $self, @args ) = @_;
    my $key = shift @args;
    my @sequences_to_export = @{$self->get_sequences_to_export()};
    return ( ( $key ~~@sequences_to_export )
            || $self->is_no_sequence_selected() );
}

sub perform_export{
    my ( $self, @args ) = @_;

    my %results = %{$self->{'results'}};
    my $output = '';
    $output .= $self->get_header();
    my @keys = $self->get_sorted_keys();
    foreach my $key (@keys) {
        if ( $self->is_sequence_to_export($key)){
            my $candidate = $results{$key};
            $output .= $self->export_candidate($candidate);
        }
    }
    return $output;
}

sub export_for_web {
    my ( $self, @args ) = @_;
    my $content_type = $self->get_content_type();
    my $content_disposition = $self->get_content_disposition();
    my $filename = $self->get_filename();
    my $contents = $self->perform_export();
    my $answer = <<"DATA"
Content-type: $content_type
Content-disposition: $content_disposition;filename=$filename

$contents
DATA
;
    return $answer;
}

sub export_on_disk {
    my ( $self, @args ) = @_;
    my $directory = shift @args;
    my $file = File::Spec->catfile($directory, $self->get_filename());
    open( my $FILE, '>', $file )
      or die("Cannot open $file: $!");
    print {$FILE} $self->perform_export()
      or die("Cannot write in file $file: $!");
    close($FILE)
      or die("Cannot close file $file: $!");
    return $file;
}

1;
