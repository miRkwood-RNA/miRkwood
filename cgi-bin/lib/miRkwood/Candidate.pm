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
use miRkwood::Components;
use miRkwood::WebPaths;

my $candidate_base_filename = 'candidate.yml';


=method retrieve_candidate_information

Check correctness and get the result for a given candidate

Arguments:
- $job - the job directory
- $id - the candidate identifier
=cut

sub retrieve_candidate_information {
    my ( $self, @args ) = @_;
    my $job = shift @args;
    my $id = shift @args;
    my $candidate_file = $self->get_candidate_filename($job, $id);
    if ( ! -e $candidate_file ){
        die("Unvalid candidate information");

    }else{
        return $self->deserialize_candidate($candidate_file);
    }
}

=method get_candidate_filename

Get the candidate filename given its identifier and the job directory

=cut

sub get_candidate_filename {
    my ( $self, @args ) = @_;
    my $job = shift @args;
    my $id = shift @args;
    return File::Spec->catfile($job, 'candidates', $self->make_candidate_filename($id));
}

=method make_candidate_filename

Return the candidate filename based on the identifier.

=cut

sub make_candidate_filename {
    my ( $self, @args ) = @_;
    my $identifier = shift @args;
    return $identifier . '.yml';
}

=method serialize_candidate_information


=cut

sub serialize_candidate_information {
    my ( $self, @args ) = @_;
    my $job_dir = shift @args;
    my $seq_dir = shift @args;
    my $can_dir = shift @args;
    my $serialization_dir = shift @args;

    my $cfg    = miRkwood->CONFIG();

    my $full_candidate_dir = miRkwood::Paths->get_candidate_paths($job_dir,  $seq_dir, $can_dir);
    my $candidate_file = File::Spec->catfile($full_candidate_dir, $candidate_base_filename);
    my %candidate = $self->parse_candidate_information($full_candidate_dir);
    $candidate{'identifier'} = "$seq_dir-$can_dir";
#    $candidate{'name'} = $seq_dir;    #récupération nom séquence

    if ( $cfg->param('options.varna') ) {
        $candidate{'image'} = File::Spec->catfile($full_candidate_dir, 'image.png');
    } else {
        $candidate{'image'} = '';
    }

    $candidate{'position'} = "$candidate{'start_position'}-$candidate{'end_position'}";
    $candidate{'length'} = $candidate{'end_position'} - $candidate{'start_position'} +1;
    $candidate{'%GC'} = miRkwood::Utils::restrict_num_decimal_digits(
                            miRkwood::Utils::compute_gc_content($candidate{'DNASequence'}),
                            3);

    my $alternative_candidates_file = File::Spec->catfile($full_candidate_dir, 'alternativeCandidates.txt');
    if (-e $alternative_candidates_file){
        my %alternatives = miRkwood::Parsers::parse_alternative_candidates_file($alternative_candidates_file);
        $candidate{'alternatives'} = \%alternatives;
    }

    my $hairpin = miRkwood::Utils::make_ASCII_viz($candidate{'DNASequence'}, $candidate{'Vienna'});
    $candidate{'hairpin'} = $hairpin;
    my %sequence;
#    $sequence{$candidate{'name'}} = $candidate{'DNASequence'};
#    my $tmp_file = File::Spec->catfile($full_candidate_dir, "mirdup_prediction.txt");
#    $candidate{'mirdup_prediction'} = \%{miRkwood::MiRdup->predict_with_mirdup($tmp_file, \%sequence)};

    return $self->serialize_candidate( \%candidate, $serialization_dir );
}

=method parse_candidate_information

Get the results for a given candidate

Arguments:
- $full_candidate_dir - the prefixed path to the candidate results

=cut

sub parse_candidate_information {
    my ( $self, @args ) = @_;
    my $full_candidate_dir = shift @args;
    my %result = ();

    my $seq_info_file =
      File::Spec->catfile( $full_candidate_dir, 'sequence_information.txt' );
    if ( -e $seq_info_file )    # si fichier existe
    {
        my @res = miRkwood::Components::get_sequence_information($seq_info_file);
        ($result{'strand'}, $result{'start_position'}, $result{'end_position'}) = @res;
    }

    my $randfold_output =
      File::Spec->catfile( $full_candidate_dir, 'randfold.out' );
    if ( -e $randfold_output )    # si fichier existe
    {
        $result{'shuffles'} = miRkwood::Parsers::parse_pvalue($randfold_output);
    }

    #Récupération valeur MFEI
    my $mfei_out =
      File::Spec->catfile( $full_candidate_dir, 'outMFEI.txt' );
    if ( -e $mfei_out )                 # si fichier existe
    {
        my @mfeis = miRkwood::Parsers::parse_mfei($mfei_out);
        $result{'mfei'} = $mfeis[0];
        $result{'mfe'} = $mfeis[1];
        $result{'amfe'} = $mfeis[2];
    }

    #Récupération séquence et format Vienna
    my $rnafold_stemloop_out = File::Spec->catfile( $full_candidate_dir,
                                       'outRNAFold_stemloop.txt' );
    if ( -e $rnafold_stemloop_out )                  # si fichier existe
    {
        my @res = miRkwood::Components::get_data_from_rnafold_out($rnafold_stemloop_out);
        my $devnull;
        ($result{'name'}, $devnull, $result{'DNASequence'}, $result{'Vienna'}) = @res;
    }

    #Récupération séquence et format Vienna
    my $rnafold_optimal_out = File::Spec->catfile( $full_candidate_dir,
                                                   'outRNAFold_optimal.txt' );
    if ( -e $rnafold_optimal_out )                  # si fichier existe
    {
        my @vienna_res = miRkwood::Parsers::parse_RNAfold_output($rnafold_optimal_out);

        $result{'Vienna_optimal'} = $vienna_res[2];
    }

    #Récupération alignement avec mirBase
    my $alignments_results_file = File::Spec->catfile($full_candidate_dir, 'merged_alignments.yml');
    my $mirdup_results_file = File::Spec->catfile($full_candidate_dir, 'mirdup_results.yml');
    $result{'alignment_existence'} = ( -e $alignments_results_file && ! -z $alignments_results_file );
    if ($result{'alignment_existence'}){
        my %mirdup_results = YAML::XS::LoadFile($mirdup_results_file) or die("Error when parsing YAML file $mirdup_results_file");
        my %alignments = YAML::XS::LoadFile($alignments_results_file);
        $result{'alignments'} = \%alignments;
        $result{'mirdup_validation'} = \%mirdup_results;
    }
    # Computing general quality
    $result{'alignment'} = $self->compute_alignment_quality(\%result);
    $result{'quality'} = $self->compute_quality(\%result);

    return %result;
}

=method serialize_candidate

Serialize the given candidate on disk

Arguments:
- $serialization_path - the filepath to serialize to
- %candidate - the candidate

=cut

sub serialize_candidate{
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my $serialization_path = shift @args;
    my $candidate_base_filename = $self->make_candidate_filename($candidate{'identifier'});
    my $serialization_file = File::Spec->catfile($serialization_path, $candidate_base_filename);
    return YAML::XS::DumpFile($serialization_file, %candidate);
}

=method deserialize_candidate

Deerialize the given candidate on disk

Arguments:
- $serialization_file - the filepath to serialize to

=cut

sub deserialize_candidate{
    my ( $self, @args ) = @_;
    my $serialization_file = shift @args;
    (-e $serialization_file)
        or die("File $serialization_file does not exists");
    return YAML::XS::LoadFile($serialization_file);
}

=method has_mirdup_validation

Return whether the given candidate has at least
one alignment which has been validated by MirDup.

=cut

sub has_mirdup_validation{
    my ($self, @args) = @_;
    my %candidate = %{shift @args};
    if ($candidate{'mirdup_validation'}){
        my %mirdup_results = %{$candidate{'mirdup_validation'}};
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
    my %candidate = %{shift @args};
    my $quality = 0;
    if ( $candidate{'mfei'} ) {
        if ( $candidate{'mfei'} < -0.8 ){
            $quality += 1;
        }
    }
    $quality += $self->compute_alignment_quality(\%candidate);
    return $quality;
}

=method compute_alignment_quality

Compute the alignment quality score

=cut

sub compute_alignment_quality {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my $alignment_existence = 0;
    if ($candidate{'alignment_existence'}){
        $alignment_existence = $candidate{'alignment_existence'};
    }else{
        $alignment_existence = 0;
    }
    my $has_mirdup_validation = $self->has_mirdup_validation(\%candidate);
    return $alignment_existence + $has_mirdup_validation;
}

=method get_absolute_image

Return the path to the image

=cut

sub get_absolute_image {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my $image = $candidate{'image'};
    return $image;
}

=method get_relative_image

Return the path to the image

=cut

sub get_relative_image {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my $image = $candidate{'image'};
    return miRkwood::WebPaths->filesystem_to_relative_path($image);
}

=method get_name

Return the name of the given candidate,
constructed by concatenating sequence name and position.

=cut

sub get_name {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    return $candidate{'name'}.'__'.$candidate{'position'};
}

=method get_shortened_sequence_name

Return the name of the sequence for the given candidate,
shortened and sanitized

=cut

sub get_shortened_sequence_name {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    return miRkwood::Utils::sanitize_sequence_name($candidate{'name'});
}

=method get_shortened_name

Return the shortened name of the candidate

=cut

sub get_shortened_name {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my $name = $self->get_shortened_sequence_name(\%candidate);
    return $name . '__' . $candidate{'position'};
}


=method candidateAsVienna

Convert a given candidate to Vienna dot-bracket format

=cut

sub candidateAsVienna {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my $optimal = shift @args;
    my $output = "";
    my $candidate_name = $self->get_name(\%candidate);
    my $header = ">$candidate_name";
    my $structure;
    if ($optimal){
        $header .= ", MFE structure";
        $structure = $candidate{'Vienna_optimal'};
    }else{
        $structure = $candidate{'Vienna'};
        $header .= ", stemloop structure";
    }
    $output .= $header . "\n" . $candidate{'DNASequence'} . "\n" . "$structure" . "\n";
    return $output;
}

=method candidateAsFasta

Convert a given candidate to FASTA format

=cut

sub candidateAsFasta {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my $output = "";
    my $candidate_name = $self->get_name(\%candidate);
    $output .= '>'.$candidate_name . "\n" . $candidate{'DNASequence'} . "\n";
    return $output;
}

=method candidate_as_gff

Convert a given candidate to a GFF line.

Usage:
my $gff_line = miRkwood::Candidate->candidate_as_gff($value);

=cut

sub candidate_as_gff {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my $candidate_name = 'preMir_' . $self->get_shortened_name(\%candidate);
    my $text .= q{} .                       # BEGIN
      $candidate{'name'} . "\t" .           # seqid
      'miRkwood' . "\t" .                   # source
      'miRNA_primary_transcript' . "\t" .   # type
      $candidate{'start_position'} . "\t" . # start
      $candidate{'end_position'} . "\t" .   # end
      '.' . "\t" .                          # score
      $candidate{'strand'} . "\t" .         # strand
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
    my %candidate = %{shift @args};
    my $alternatives = $candidate{'alternatives'};
    my $output = "";
    if ($alternatives) {
        my %alternatives = %{$alternatives};
        foreach my $name (keys %alternatives) {
            my %alternative = %{$alternatives{$name}};
            $output .= '>'.$name . ' (MFEI: ' . $alternative{'mfei'} . ')'. "\n" .
                       $alternative{'sequence'} . "\n" .
                       $alternative{'structure'} . "\n";
        }
    }
    return $output;
}

=method make_Vienna_viz

Make a nicer Vienna display by cutting too long lines.

Usage:
my $string = make_Vienna_viz($Vienna, $DNASequence)

=cut

sub make_Vienna_viz {
    my ($self, @args) = @_;
    my $Vienna = shift @args;
    my $DNASequence = shift @args;

    my $viennaString   = q{};
    my $sequenceString = q{};
    my $string         = q{};
    for ( 1 .. length($Vienna) ) {

        $viennaString   .= substr $Vienna,      $_ - 1, 1;
        $sequenceString .= substr $DNASequence, $_ - 1, 1;
        if ( $_ % 60 == 0 ) {

            $string .= $sequenceString . "\n" . $viennaString . "\n\n";
            $viennaString   = q{};
            $sequenceString = q{};
        }
        if ( ( $viennaString ne q{} ) && ( $_ == length($Vienna) ) ) {
            $string .= $sequenceString . "\n" . $viennaString . "\n\n";
        }
    }
    return $string
}


=method make_alignments_HTML


=cut

sub make_alignments_HTML {
    my ($self, @args) = @_;
    my %candidate = %{shift @args};

    # Alignments

    my %alignments = %{$candidate{'alignments'}};
    my %mirdup_results = %{$candidate{'mirdup_validation'}};

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
        my $mirdup_key = $candidate{'name'} . '__' . $position;
        my $mirdup_prediction;
        if ( $mirdup_results{$mirdup_key} ){
            $mirdup_prediction = 'This prediction is validated by miRdup.';
        } else {
            $mirdup_prediction = 'This prediction is not validated by miRdup.';
        }

        # Hairpin
        my $hairpin_with_mature =
            miRkwood::Utils::make_hairpin_with_mature($candidate{'hairpin'},
                                                           $left, $right,
                                                           length $candidate{'DNASequence'},
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
my $xml_element = miRkwood::Candidate->candidate_as_pseudoXML($value);

=cut

sub candidate_as_pseudoXML {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};

    my $name = $self->get_shortened_sequence_name(\%candidate);

    my @fields_to_truncate = ( 'mfe', 'mfei', 'amfe' );

    my @optional_fields = $self->get_optional_candidate_fields();
    my @headers1        =
      ( 'position', 'length', 'strand', 'quality', @optional_fields );
    my @headers2 = ( 'Vienna', 'DNASequence', 'identifier' );

    my $result = "<Sequence";

    $result .= " name='$name'";
    for my $header (@headers1) {
        my $contents = $candidate{$header};
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
    my $img = $self->get_relative_image(\%candidate);
    $result .= " image='$img'";
    for my $header (@headers2) {
        $result .= " $header='$candidate{$header}'";
    }
    $result .= "></Sequence>";

    return $result;
}

1;
