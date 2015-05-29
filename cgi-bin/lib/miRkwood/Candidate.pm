package miRkwood::Candidate;

# ABSTRACT: Code directly tied to the candidate data structure

use strict;
use warnings;

use File::Spec;
use YAML::XS;

use Log::Message::Simple qw[msg error debug];

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
    my $cfg = miRkwood->CONFIG();
    my $mode = $cfg->param('job.mode');

    if ( $self->{'mfei'} ) {
        if ( $self->{'mfei'} < -0.8 ){
            $quality += 1;
        }
    }

    if ( $mode eq 'WebBAM' ){
        my ($end_arm_1, $start_arm_2) = $self->determine_precursor_arms( );
        $quality += $self->compute_quality_from_reads( $end_arm_1, $start_arm_2 );
        $quality += $self->has_mirdup_validation();
    }
    else{
        $quality += $self->{'alignment'};
    }

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
    my $output = '';
    my $candidate_name = $self->get_name();
    my $header = ">$candidate_name";
    my $structure;
    if ($optimal){
        $header .= ', MFE structure';
        $structure = $self->{'structure_optimal'};
    }else{
        $structure = $self->{'structure_stemloop'};
        $header .= ', stemloop structure';
    }
    $output .= $header . "\n" . $self->{'sequence'} . "\n" . "$structure" . "\n";
    return $output;
}

=method candidateAsFasta

Convert a given candidate to FASTA format

=cut

sub candidateAsFasta {
    my ( $self, @args ) = @_;
    my $output = '';
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
    my $text = q{} .                        # BEGIN
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
    my $output = '';
    if ($alternatives) {
        my %alternatives = %{$alternatives};
        foreach my $name (sort keys %alternatives) {
            my %alternative = %{$alternatives{$name}};
            $output .= '>'.$name . ' (MFEI: ' . $alternative{'mfei'} . ')'. "\n" .
                       $alternative{'sequence'} . "\n" .
                       $alternative{'structure_stemloop'} . "\n";
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

    my $contents = '';
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
            my $additional_content = '';
            my $name;
            my $html_name;
            my @splitted = split('\|', $names);

            my $spacing = 15;
            my ($top, $middle, $bottom) = split(/\n/, $alignment);
            $top    = sprintf "%-${spacing}s %3s %s %s", 'query', $hit->{'begin_target'}, $top,   $hit->{'end_target'};
            $middle = sprintf "%-${spacing}s %3s %s %s", '',      '',                     $middle, '';

            my $mirbase_title = '';
            if( (scalar @splitted) > 1 ) {
                $mirbase_title = 'miRBase sequences: ';
            }else{
                $mirbase_title = 'miRBase sequence: ';
            }

            $name = 'miRBase';
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
            $additional_content = "<span class='others'>$mirbase_title" . join(', ', @sequences) . "</span>";


            $bottom = sprintf "%-${spacing}s %3s %s %s", $name,   $hit->{'begin_query'},  $bottom, $hit->{'end_query'};
            my $additional_space = '';
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
    my $toc = "<span class='toc'>Putative miRNA locus: " . join(', ', @TOC) . '</span>';
    return $toc . "\n" . $contents;

}

=method get_optional_candidate_fields

Return the optional fields based on the current configuration

=cut

sub get_optional_candidate_fields {
    my ( $self, @args ) = @_;
    my @fields = ();
    my $cfg    = miRkwood->CONFIG();
    push @fields, qw{mfe mfei amfe};
    if ( $cfg->param('options.randfold') ) {
        push @fields, ('shuffles');
    }
    if ( $cfg->param('options.align') ) {
        push @fields, ('alignment');
    }
    if ( defined($cfg->param('job.mode')) and ($cfg->param('job.mode') eq 'bam' or $cfg->param('job.mode') eq 'WebBAM') ){
        push @fields, ('reads');
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

    my @fields_to_truncate = qw{mfe mfei amfe};

    my @optional_fields = $self->get_optional_candidate_fields();
    my @headers1        =
      ( 'position', 'length', 'strand', 'quality', @optional_fields );
    my @headers2 = qw{structure_stemloop sequence identifier};
    my $result = "<Sequence";

    $result .= " name='$name'";
    for my $header (@headers1) {
        my $contents = $self->{$header};
        if (grep { $header eq $_ } @fields_to_truncate){
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

sub get_basic_informations {

	my ( $self, @args ) = @_;

	my @optional_fields = $self->get_optional_candidate_fields();
    my @headers;
    if ( defined( $self->{'precursor_name'} ) ){
        push @headers, 'precursor_name';
    }
    push @headers, ( 'identifier', 'position', 'start_position', 'length', 'strand', 'quality', @optional_fields );
	my $result = {};

	foreach (@headers){
		$result->{$_} = $self->{$_};
	}
	$result->{'name'} = $self->get_shortened_sequence_name();
	$result->{'image'} = $self->get_relative_image();
	return $result;

}

=method turn_relative_positions_into_absolute_positions

=cut
sub turn_relative_positions_into_absolute_positions {
    my ($self) = @_;

    my $locus = $self->{'name'};

    my $chromosome = '';
    my $start_cluster = -1;
    my $end_cluster = -1;
    if ( $locus =~ /([^_]+)__(\d+)-(\d+)/ ){
        $chromosome = $1;
        $start_cluster = $2;
        $end_cluster = $3;
    }

    $self->{'start_position'} = $start_cluster + $self->{'start_position'} -1;
    $self->{'end_position'}   = $start_cluster + $self->{'end_position'};
    $self->{'position'}       = $self->{'start_position'} . '-' . $self->{'end_position'};
    $self->{'name'}           = $chromosome;

    return $self;

}

=method get_reads_from_bam_file

Method to get all reads corresponding to a candidate sequence
Parameter : bam file full path
Modifies the object candidate to add a variable (hash) containing
all reads with start, stop and depth.

=cut
sub get_reads_from_bam_file {

    my ( $self, @args ) = @_;
    my $bam_file = shift @args;

    my $candidate_position = $self->{'name'} . ':' . $self->{'start_position'} . '-' . $self->{'end_position'};
    #~ $candidate_position =~ s/__/:/;

    my $reads = {};

    my $samtools_cmd = "samtools view $bam_file $candidate_position |";

    open(my $SAMVIEW, $samtools_cmd);
    while ( <$SAMVIEW> ){
        my @line = split(/\t/);
        my $start = $line[3];
        my $end = length($line[9]) + $start - 1;

        # Count the depth of each read
        if ( !exists( $reads->{"$start-$end"} ) ){
            $reads->{"$start-$end"} = 0;
        }
        $reads->{"$start-$end"}++;

    }
    close $SAMVIEW;

    $self->{'reads'} = $reads;

    return $self;

}

=method get_reads_from_bed_file

Method to get all reads corresponding to a candidate sequence
Parameter : bed file full path
Modifies the object candidate to add a variable (hash) containing
all reads with start, stop and depth.

=cut
sub get_reads_from_bed_file {
    my ( $self, @args ) = @_;
    my $bed_file = shift @args;
    my $reads = {};

    open (my $BED, '<', $bed_file) or die "ERROR while opening $bed_file : $!";
    while ( <$BED> ){
        chomp;
        my @fields = split( /\t/ );
        if ( $fields[0] eq $self->{'name'} && $fields[5] eq $self->{'strand'} ){
            my $start_read = $fields[1] + 1;
            my $end_read   = $fields[2];
            if ( $start_read >= $self->{'start_position'} && $end_read <= $self->{'end_position'} ){
                $reads->{"$start_read-$end_read"} = $fields[4];
            }
        }
    }

    $self->{'reads'} = $reads;

    return $self;

}

=method determine_precursor_arms

  Method to determine the positions of the precursor
  based on positions of last '(' and first ')' in the
  structure

=cut
sub determine_precursor_arms {
    my ( $self, @args ) = @_;

    my $start_arm_2 = 0;
    my $end_arm_1 = 0;
    my @structure_array = split('', $self->{'structure_stemloop'});
    my $stop = 0;

    my $index = 0;

    while ( $index < scalar(@structure_array) and ! $stop ){
        if ( $structure_array[$index] eq '(' ){
            $end_arm_1 = $index;
        }
        elsif ( $structure_array[$index] eq ')' ){
            $start_arm_2 = $index;
            $stop = 1;
        }
        $index++;
    }

    $start_arm_2 += $self->{'start_position'};
    $end_arm_1   += $self->{'start_position'};

    return ( $end_arm_1, $start_arm_2 );
}

=method compute_quality_from_reads

=cut
sub compute_quality_from_reads {
    my ( $self, @args ) = @_;
    my $end_arm_1 = shift @args;
    my $start_arm_2 = shift @args;

    my $start_arm_1 = $self->{'start_position'};
    my $end_arm_2   = $self->{'end_position'};
    my $reads_arm_1 = {};
    my $reads_arm_2 = {};
    my $count_arm_1 = 0;
    my $count_arm_2 = 0;
    my $total_reads = 0;
    my $criteria_nb_reads = 0;
    my $criteria_duplex = 0;
    my $criteria_precision_arm_1 = 0;
    my $criteria_precision_arm_2 = 0;
    my $max_count_arm_1 = 0;
    my $max_count_arm_2 = 0;
    my $max_read_arm_1 = '';
    my $max_read_arm_2 = '';

    if ( scalar(keys %{$self->{'reads'}}) > 0 ){
        foreach my $read_position (keys %{$self->{'reads'}}){
            my ($start_read, $end_read) = split(/-/, $read_position);
            if ( $start_read >= $start_arm_1 and $end_read <= $end_arm_1 ){
                $reads_arm_1->{ $start_read } += $self->{'reads'}{ $read_position };
                $count_arm_1 += $self->{'reads'}{ $read_position };
                if ( $self->{'reads'}{ $read_position } > $max_count_arm_1 ){
                    $max_read_arm_1  = $read_position;
                    $max_count_arm_1 = $self->{'reads'}{ $read_position };
                }
            }
            if ( $start_read >= $start_arm_2 and $end_read <= $end_arm_2 ){
                $reads_arm_2->{ $start_read } += $self->{'reads'}{ $read_position };
                $count_arm_2 += $self->{'reads'}{ $read_position };
                if ( $self->{'reads'}{ $read_position } > $max_count_arm_2 ){
                    $max_read_arm_2  = $read_position;
                    $max_count_arm_2 = $self->{'reads'}{ $read_position };
                }
            }
            $total_reads += $self->{'reads'}{ $read_position };
        }
    }

    # Criteria nb of reads
    if ( $count_arm_1 >= 10 and $count_arm_2 >= 10 ){
        $criteria_nb_reads = 1;
    }
    if ( $count_arm_1 >= 5 and $count_arm_2 >= 5 and $total_reads >= 100 ){
        $criteria_nb_reads = 1;
    }

    # Criteria duplex
    if ( $max_read_arm_1 ne '' && $max_read_arm_2 ne '' ){

        my ($start_max_1, $end_max_1) = split(/-/, $max_read_arm_1);
        my ($start_max_2, $end_max_2) = split(/-/, $max_read_arm_2);
        my $relative_start_max_1 = $start_max_1 - $start_arm_1;
        my $relative_end_max_2   = $end_max_2 - $start_arm_1;

        if ( $self->{'CT'}{($relative_start_max_1)} ){
            if ( $self->{'CT'}{$relative_start_max_1} - $relative_end_max_2 >= 0
             and $self->{'CT'}{$relative_start_max_1} - $relative_end_max_2 <= 4 ){
                $criteria_duplex = 1;
            }
        }
        elsif ( $self->{'CT'}{($relative_start_max_1 +1 )} ) {
            if ( $self->{'CT'}{$relative_start_max_1 +1} - ($relative_end_max_2 -1) >= 0
             and $self->{'CT'}{$relative_start_max_1 +1} - ($relative_end_max_2 -1) <= 4 ){
                $criteria_duplex = 1;
            }
        }
        elsif( $self->{'CT'}{($relative_start_max_1 -1 )} ) {
            if ( $self->{'CT'}{$relative_start_max_1 -1} - ($relative_end_max_2 +1) >= 0
             and $self->{'CT'}{$relative_start_max_1 -1} - ($relative_end_max_2 +1) <= 4 ){
                $criteria_duplex = 1;
            }
        }
    }

    # Criteria precision
    foreach my $start_read ( keys %{$reads_arm_1} ){
        if ( $reads_arm_1->{$start_read} >= ( $count_arm_1 / 2 ) ){
            $criteria_precision_arm_1 = 1;
        }
    }
    foreach my $start_read ( keys %{$reads_arm_2} ){
        if ( $reads_arm_2->{$start_read} >= ( $count_arm_2 / 2 ) ){
            $criteria_precision_arm_2 = 1;
        }
    }

    return $criteria_nb_reads + $criteria_duplex + ( $criteria_precision_arm_1 * $criteria_precision_arm_2 );
}


sub store_attribute_ct {
    my ( @args ) = @_;
    my $candidate = shift @args;
    my $directory = shift @args;

    open (my $CT, '<', "$directory/outB2ct_stemloop.ct")
        or die "Error while opening $directory/outB2ct_stemloop.ct : $!";
    my @ct = <$CT>;
    close $CT;
    foreach ( @ct ){
        if ( /(\d+)\s+[a-zA-Z]\s+\d+\s+\d+\s+(\d+)\s+\d+/ ){
            $candidate->{'CT'}{ $1 } = $2;
        }
    }

    return $candidate;

}

=method create_reads_length_diagramm

  Method to draw a diagramm of reads length in raw text.
  
=cut
sub create_reads_length_diagramm {
    my ($self, @args) = @_;
    my $max_width = 100;

    my $diagramm = '<pre>';
    my %reads_length;

    foreach ( keys%{$self->{'reads'}} ){
        my @tab = split ( /-/ );
        my $length = $tab[1] - $tab[0] + 1;
        if ( !defined( $reads_length{ $length } ) ){
            $reads_length{ $length } = 0;
        }
        $reads_length{ $length } += $self->{'reads'}{$_};        
    }

    my $max = 0;
    foreach my $key ( keys%reads_length ){
        if ( $reads_length{$key} > $max ){
            $max = $reads_length{$key};
        }
    }

    foreach my $key ( sort( keys%reads_length ) ){
        my $width = int( $reads_length{$key} * $max_width / $max + 0.5 );
        $diagramm .= "$key nt ($reads_length{$key})\t| ";
        my $i = 0;
        while ( $i < $width ){
            $diagramm .= '*';
            $i++;
        }
        $diagramm .= "\n";
    }
    $diagramm .= '</pre>';

    return $diagramm;
}


1;
