package PipelineMiRNA::Candidate;

# ABSTRACT: Code directly tied to the candidate data structure

use strict;
use warnings;

use Data::Dumper;
use File::Spec;
use YAML::XS;

use PipelineMiRNA;
use PipelineMiRNA::Paths;
use PipelineMiRNA::MiRdup;
use PipelineMiRNA::Parsers;
use PipelineMiRNA::Utils;
use PipelineMiRNA::WebTemplate;
use PipelineMiRNA::Components;


my $candidate_base_filename = 'candidate.yml';


=method retrieve_candidate_information

Check correctness and get the result for a given candidate

Arguments:
- $job - the job identifier
- $id - the candidate identifier
=cut

sub retrieve_candidate_information {
    my ( $self, @args ) = @_;
    my $job = shift @args;
    my $id = shift @args;
    my $candidate_file = File::Spec->catfile($job, 'candidates', $self->make_candidate_filename($id));
    if ( ! -e $candidate_file ){
        die("Unvalid candidate information");

    }else{
        return $self->deserialize_candidate($candidate_file);
    }
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

    my $full_candidate_dir = PipelineMiRNA::Paths->get_candidate_paths($job_dir,  $seq_dir, $can_dir);
    my $candidate_file = File::Spec->catfile($full_candidate_dir, $candidate_base_filename);
    my %candidate = $self->parse_candidate_information($full_candidate_dir);
    $candidate{'identifier'} = "$seq_dir-$can_dir";
#    $candidate{'name'} = $seq_dir;    #récupération nom séquence
    $candidate{'image'} = File::Spec->catfile($full_candidate_dir, 'image.png');

    $candidate{'length'} = PipelineMiRNA::Utils::get_element_of_split($candidate{'position'},'-', 1) - PipelineMiRNA::Utils::get_element_of_split($candidate{'position'},'-', 0) +1;


    my $alternative_candidates_file = File::Spec->catfile($full_candidate_dir, 'alternativeCandidates.txt');
    if (-e $alternative_candidates_file){
        my %alternatives = PipelineMiRNA::Parsers::parse_alternative_candidates_file($alternative_candidates_file);
        $candidate{'alternatives'} = \%alternatives;
    }

    my $hairpin = PipelineMiRNA::Utils::make_ASCII_viz($candidate{'DNASequence'}, $candidate{'Vienna'});
    $candidate{'hairpin'} = $hairpin;
    my %sequence;
#    $sequence{$candidate{'name'}} = $candidate{'DNASequence'};
#    my $tmp_file = File::Spec->catfile($full_candidate_dir, "mirdup_prediction.txt");
#    $candidate{'mirdup_prediction'} = \%{PipelineMiRNA::MiRdup->predict_with_mirdup($tmp_file, \%sequence)};

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
    my $pvalue =
      File::Spec->catfile( $full_candidate_dir, 'pvalue.txt' );
    if ( -e $pvalue )    # si fichier existe
    {
        $result{'p_value'} = PipelineMiRNA::Parsers::parse_pvalue($pvalue);
    }

    #Récupération valeur MFEI
    my $mfei_out =
      File::Spec->catfile( $full_candidate_dir, 'outMFEI.txt' );
    if ( -e $mfei_out )                 # si fichier existe
    {
        my @mfeis = PipelineMiRNA::Parsers::parse_mfei($mfei_out);
        $result{'mfei'} = $mfeis[0];
        $result{'mfe'} = $mfeis[1];
        $result{'amfe'} = $mfeis[2];
    }

    #Récupération séquence et format Vienna
    my $rnafold_stemloop_out = File::Spec->catfile( $full_candidate_dir,
                                       'outRNAFold_stemloop.txt' );
    if ( -e $rnafold_stemloop_out )                  # si fichier existe
    {
        my @res = PipelineMiRNA::Components::get_data_from_rnafold_out($rnafold_stemloop_out);
        ($result{'name'}, $result{'position'}, $result{'DNASequence'}, $result{'Vienna'}) = @res;

        my @vienna_res = PipelineMiRNA::Parsers::parse_RNAfold_output($rnafold_stemloop_out);

    }

    #Récupération séquence et format Vienna
    my $rnafold_optimal_out = File::Spec->catfile( $full_candidate_dir,
                                                   'outRNAFold_optimal.txt' );
    if ( -e $rnafold_optimal_out )                  # si fichier existe
    {
        my @vienna_res = PipelineMiRNA::Parsers::parse_RNAfold_output($rnafold_optimal_out);

        $result{'Vienna_optimal'} = $vienna_res[2];
    }

    #Récupération alignement avec mirBase
    my $file_alignement = File::Spec->catfile($full_candidate_dir, 'alignement.txt');
    $result{'alignment_existence'} = ( -e $file_alignement && ! -z $file_alignement );
    my %alignments;
    if (! eval {%alignments = PipelineMiRNA::Components::parse_custom_exonerate_output($file_alignement);}) {
        # Catching exception
    } else {
        %alignments = $self->merge_alignments(\%alignments);
        my $tmp_file = File::Spec->catfile($full_candidate_dir, "mirdup_validation.txt");
        my %mirdup_results = PipelineMiRNA::MiRdup->validate_with_mirdup($tmp_file, $result{'name'},
                                                                         $result{'DNASequence'}, $result{'Vienna'},
                                                                         keys %alignments);
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
    my %mirdup_results = %{$candidate{'mirdup_validation'}};
    if (scalar (grep { /^1$/ } values %mirdup_results ) >= 1){
        return 1;
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
        if ( $candidate{'mfei'} < -0.7 ){
            $quality += 1;
        }
    }
    my $length = length ($candidate{'DNASequence'});

    if ( $length > 80 && $length < 200 ){
        $quality += 1;
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
    return PipelineMiRNA::Paths->filesystem_to_relative_path($image);
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
        if ( $_ % 50 == 0 ) {

            $string .= $viennaString . "\n" . $sequenceString . "\n\n";
            $viennaString   = q{};
            $sequenceString = q{};
        }
        if ( ( $viennaString ne q{} ) && ( $_ == length($Vienna) ) ) {
            $string .= $viennaString . "\n" . $sequenceString . "\n\n";
        }
    }
    return $string
}

=method merge_alignments

Merge overlapping alignments.
Given ordered positions, merge in [a..b] all [c..d] if a<=c and d<=b+2

=cut

sub merge_alignments {
    my ($self, @args) = @_;
    my $alignments = shift @args;
    my %alignments = %{$alignments};

    my %merged_alignments;
    my ($stocked_left, $stocked_right) = (-10, -10);

    my @keys = sort { ( PipelineMiRNA::Utils::get_element_of_split($a, '-', 0) <=>
                        PipelineMiRNA::Utils::get_element_of_split($b, '-', 0) )
                      ||
                      ( PipelineMiRNA::Utils::get_element_of_split($a, '-', 1)  <=>
                        PipelineMiRNA::Utils::get_element_of_split($b, '-', 1) )
                    } keys %alignments;
    my @stocked_hits;
    my $final_key;
    my $final_hit_count = -1;
    foreach my $current_key (@keys) {
        my ($current_left, $current_right) = split(/-/, $current_key);
        my $current_hit_count = scalar @{$alignments{$current_key}};

        if ( $current_right > $stocked_right + 3 ) {
            # No merge ; drop the list of current hits in the hash (only if there are some)
            if (@stocked_hits){
                push @{$merged_alignments{$final_key}}, @stocked_hits;
            }

            # Reinitialise
            $final_hit_count = -1;
            @stocked_hits = ();
            ($stocked_left, $stocked_right) = ($current_left, $current_right);
        }
        if ($current_hit_count > $final_hit_count){
            # This position holds more hits than the previous, so it will be our new key
            $final_hit_count = $current_hit_count;
            $final_key = $current_key;
        }
        # Stock the current hits
        push @stocked_hits, @{$alignments{$current_key}};
    }
    # Drop the remaining hits in the hash
    push @{$merged_alignments{$final_key}}, @stocked_hits;
    return %merged_alignments;
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
    my @keys = sort { ( PipelineMiRNA::Utils::get_element_of_split($a, 0)  <=>
                        PipelineMiRNA::Utils::get_element_of_split($b, 0)
                      ) ||
                      ( PipelineMiRNA::Utils::get_element_of_split($a, 1)  <=>
                        PipelineMiRNA::Utils::get_element_of_split($b, 1))
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
            PipelineMiRNA::Utils::make_hairpin_with_mature($candidate{'hairpin'},
                                                           $left, $right,
                                                           length $candidate{'DNASequence'});
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
                    my $mirbase_link = PipelineMiRNA::WebTemplate::make_mirbase_link($mirbase_id);
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

1;
