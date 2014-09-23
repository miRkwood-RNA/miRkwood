package miRkwood::Candidate;

# ABSTRACT: Code directly tied to the candidate data structure

use strict;
use warnings;

use File::Spec;
use YAML::XS;

use miRkwood;
use miRkwood::Paths;
use miRkwood::Parsers;
use miRkwood::Utils;
use miRkwood::WebPaths;

=method new

Constructor
    my $attributes = @args;
=cut

sub new {
    my ( $class, @args ) = @_;
    my %attributes = ();
    if(@args){
    my $attributes = shift @args;
        %attributes = %{$attributes};
    }
    my $self = bless \%attributes, $class;
    return $self;
}

=method new_from_serialized

Constructor and deserializer

=cut

sub new_from_serialized {
    my ( $class, @args ) = @_;
    my $serialization_file = shift @args;
    (-e $serialization_file && -f $serialization_file)
        or die("File $serialization_file does not exist");
    my %attributes = YAML::XS::LoadFile($serialization_file);
    (%attributes) or die("Desarialization of $serialization_file failed");
    return $class->new(\%attributes);
}

=method get_identifier

Accessor for the identifier

=cut

sub get_identifier{
    my ( $self, @args ) = @_;
    return $self->{'identifier'};
}

=method has_mirdup_validation

Return whether the given candidate has at least
one alignment which has been validated by MirDup.

=cut

sub has_mirdup_validation{
    my ($self, @args) = @_;
    if ($self->{'mirdup_validation'}){
        my %mirdup_results = %{$self->{'mirdup_validation'}};
        if (scalar (grep { /^1$/ } values %mirdup_results ) >= 1){
            return 1;
        }
    }else{
        return 0;
    }
}

=method compute_quality

Compute a general quality score

=cut

sub compute_quality {
    my ( $self, @args ) = @_;
    my $quality = 0;
    if ( $self->{'mfei'} ) {
        if ( $self->{'mfei'} < -0.8 ){
            $quality += 1;
        }
    }
    $quality += $self->{'alignment'};
    $self->{'quality'} =  $quality;
    return;
}

=method compute_alignment_quality

Compute the alignment quality score

=cut

sub compute_alignment_quality {
    my ( $self, @args ) = @_;
    my $alignment_existence = 0;
    if ($self->{'alignment_existence'}){
        $alignment_existence = $self->{'alignment_existence'};
    }else{
        $alignment_existence = 0;
    }
    my $has_mirdup_validation = $self->has_mirdup_validation();
    $self->{'alignment'} = $alignment_existence + $has_mirdup_validation;
    return;
}

=method get_absolute_image

Return the path to the image

=cut

sub get_absolute_image {
    my ( $self, @args ) = @_;
    my $image = $self->{'image'};
    return $image;
}

=method get_relative_image

Return the path to the image

=cut

sub get_relative_image {
    my ( $self, @args ) = @_;
    my $image = $self->get_absolute_image();
    return miRkwood::WebPaths->filesystem_to_relative_path($image);
}

=method get_name

Return the name of the given candidate,
constructed by concatenating sequence name and position.

=cut

sub get_name {
    my ( $self, @args ) = @_;
    return $self->{'name'}.'__'.$self->{'position'};
}

=method get_shortened_sequence_name

Return the name of the sequence for the given candidate,
shortened and sanitized

=cut

sub get_shortened_sequence_name {
    my ( $self, @args ) = @_;
    my $name = $self->{'name'};
    return miRkwood::Utils::sanitize_sequence_name($name);
}

=method get_shortened_name

Return the shortened name of the candidate

=cut

sub get_shortened_name {
    my ( $self, @args ) = @_;
    my $name = $self->get_shortened_sequence_name();
    return $name . '__' . $self->{'position'};
}


=method candidateAsVienna

Convert a given candidate to Vienna dot-bracket format

=cut

sub candidateAsVienna {
    my ( $self, @args ) = @_;
    my $optimal = shift @args;
    my $output = "";
    my $candidate_name = $self->get_name();
    my $header = ">$candidate_name";
    my $structure;
    if ($optimal){
        $header .= ", MFE structure";
        $structure = $self->{'structure_optimal'};
    }else{
        $structure = $self->{'structure_stemloop'};
        $header .= ", stemloop structure";
    }
    $output .= $header . "\n" . $self->{'sequence'} . "\n" . "$structure" . "\n";
    return $output;
}

=method candidateAsFasta

Convert a given candidate to FASTA format

=cut

sub candidateAsFasta {
    my ( $self, @args ) = @_;
    my $output = "";
    my $candidate_name = $self->get_name();
    $output .= '>'.$candidate_name . "\n" . $self->{'sequence'} . "\n";
    return $output;
}

=method candidate_as_gff

Convert a given candidate to a GFF line.

Usage:
my $gff_line = $candidate->candidate_as_gff();

=cut

sub candidate_as_gff {
    my ( $self, @args ) = @_;
    my $candidate_name = 'preMir_' . $self->get_shortened_name();
    my $text .= q{} .                       # BEGIN
      $self->{'name'} . "\t" .              # seqid
      'miRkwood' . "\t" .                   # source
      'miRNA_primary_transcript' . "\t" .   # type
      $self->{'start_position'} . "\t" .    # start
      $self->{'end_position'} . "\t" .      # end
      '.' . "\t" .                          # score
      $self->{'strand'} . "\t" .            # strand
      '.' . "\t" .                          # phase
      "Name=$candidate_name" . "\t" .       # attributes
      "\n";
    return $text;
}


=method alternativeCandidatesAsVienna

Return alternative candidates as Vienna dot-bracket

=cut

sub alternativeCandidatesAsVienna {
    my ( $self, @args ) = @_;
    my $alternatives = $self->{'alternatives'};
    my $output = "";
    if ($alternatives) {
        my %alternatives = %{$alternatives};
        foreach my $name (sort keys %alternatives) {
            my %alternative = %{$alternatives{$name}};
            $output .= '>'.$name . ' (MFEI: ' . $alternative{'mfei'} . ')'. "\n" .
                       $alternative{'sequence'} . "\n" .
                       $alternative{'structure'} . "\n";
        }
    }
    return $output;
}

=method make_alignments_HTML


=cut

sub make_alignments_HTML {
    my ($self, @args) = @_;

    # Alignments

    my %alignments = %{$self->{'alignments'}};
    my %mirdup_results = %{$self->{'mirdup_validation'}};

    my $contents = "";
    my @TOC;
    my $predictionCounter = 0;

    # Sorting by position
    my @keys = sort { ( miRkwood::Utils::get_element_of_split($a, '-', 0)  <=>
                        miRkwood::Utils::get_element_of_split($b, '-', 0)
                      ) ||
                      ( miRkwood::Utils::get_element_of_split($a, '-', 1)  <=>
                        miRkwood::Utils::get_element_of_split($b, '-', 1))
                    } keys %alignments;

    foreach my $position (@keys) {
        my ($left, $right) = split(/-/, $position);

        # MiRdup
        my $mirdup_key = $self->{'name'} . '__' . $position;
        my $mirdup_prediction;
        if ( $mirdup_results{$mirdup_key} ){
            $mirdup_prediction = 'This prediction is validated by miRdup.';
        } else {
            $mirdup_prediction = 'This prediction is not validated by miRdup.';
        }

        # Hairpin
        my $hairpin_with_mature =
            miRkwood::Utils::make_hairpin_with_mature($self->{'hairpin'},
                                                           $left, $right,
                                                           length $self->{'sequence'},
                                                           'html');
        $predictionCounter += 1;

        # Sorting the hit list by descending value of the 'score' element
        my @hits = sort { $b->{'score'} <=> $a->{'score'} } @{$alignments{$position}};
        my $title = "Prediction $predictionCounter: $position";
        $contents .= "<h3 id='$position'>$title</h3>
        <pre style='height: 80px;'>$hairpin_with_mature</pre>
        $mirdup_prediction
        <h4>Alignments</h4>
        ";

        my $toc_element = "<a href='#$position'>$position</a>";
        push @TOC, $toc_element;
        foreach my $hit (@hits){
            my $alignment = $hit->{'alignment'};
            my $names = $hit->{'name'} . q{ } . $hit->{'def_query'};
            my $additional_content = "";
            my $name;
            my $html_name;
            my @splitted = split('\|', $names);

            my $spacing = 15;
            my ($top, $middle, $bottom) = split(/\n/, $alignment);
            $top    = sprintf "%-${spacing}s %3s %s %s", 'query', $hit->{'begin_target'}, $top,   $hit->{'end_target'};
            $middle = sprintf "%-${spacing}s %3s %s %s", '',      '',                     $middle, '';

            my $title = '';
            if( (scalar @splitted) > 1 ) {
                $title = 'miRBase sequences: ';
            }else{
                $title = 'miRBase sequence: ';
            }

            $name = "miRBase";
            my @sequences;
            foreach my $seq (@splitted){
                $seq =~ s/^\s+//;
                $seq =~ s/\s+$//;
                if ($seq =~ 'revcomp'){
                } else {
                    my @splitted_one = split(/ /, $seq);
                    my $name = $splitted_one[0];
                    my $mirbase_id = $splitted_one[1];
                    my $mirbase_link = miRkwood::Utils::make_mirbase_link($mirbase_id);
                    my $html_name = "<a href='$mirbase_link'>$name</a>";
                    push @sequences, $html_name;
                }
            }
            $additional_content = "<span class='others'>$title" . join(', ', @sequences) . "</span>";


            $bottom = sprintf "%-${spacing}s %3s %s %s", $name,   $hit->{'begin_query'},  $bottom, $hit->{'end_query'};
            my $additional_space = "";
            my $sub_string = substr($bottom, 0, $spacing);
            $additional_space .= ' ' while ($sub_string =~ m/ /g);
            substr($bottom, 0, $spacing) = $name . $additional_space;
            $contents .= <<"INNER";
<pre class='alignment'>
$top
$middle
$bottom
</pre>
$additional_content
INNER
        }

    }
    my $toc = "<span class='toc'>Putative miRNA locus: " . join(", ", @TOC) . '</span>';
    return $toc . "\n" . $contents;

}

=method get_optional_candidate_fields

Return the optional fields based on the current configuration

=cut

sub get_optional_candidate_fields {
    my ( $self, @args ) = @_;
    my @fields = ();
    my $cfg    = miRkwood->CONFIG();
    push @fields, ( 'mfe', 'mfei', 'amfe' );
    if ( $cfg->param('options.randfold') ) {
        push @fields, ('shuffles');
    }
    if ( $cfg->param('options.align') ) {
        push @fields, ('alignment');
    }
    return @fields;
}

=method candidate_as_pseudoXML

Convert a given candidate to a pseudo XML

Usage:
my $xml_element = $candidate->candidate_as_pseudoXML();

=cut

sub candidate_as_pseudoXML {
    my ( $self, @args ) = @_;

    my $name = $self->get_shortened_sequence_name();

    my @fields_to_truncate = ( 'mfe', 'mfei', 'amfe' );

    my @optional_fields = $self->get_optional_candidate_fields();
    my @headers1        =
      ( 'position', 'length', 'strand', 'quality', @optional_fields );
    my @headers2 = ( 'structure_stemloop', 'sequence', 'identifier' );

    my $result = "<Sequence";

    $result .= " name='$name'";
    for my $header (@headers1) {
        my $contents = $self->{$header};
        if ( $header ~~ @fields_to_truncate){
            $contents = miRkwood::Utils::restrict_num_decimal_digits($contents, 3);
        }
        if ( !defined $contents ) {
            $contents = q{};
        }
        if ( $header eq 'shuffles' && $contents == 1){
            $contents = q{};
        }
        $result .= " $header='$contents'";
    }
    my $img = $self->get_relative_image();
    $result .= " image='$img'";
    for my $header (@headers2) {
        $result .= " $header='$self->{$header}'";
    }
    $result .= "></Sequence>";

    return $result;
}

1;
