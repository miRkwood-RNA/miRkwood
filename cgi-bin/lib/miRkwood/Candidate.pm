package miRkwood::Candidate;

# ABSTRACT: Code directly tied to the candidate data structure

use strict;
use warnings;

use File::Spec;
use YAML::XS;

use Log::Message::Simple qw[msg error debug];

use miRkwood;
use miRkwood::CLI;
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
    }
    return 0;
}

=method compute_quality

Compute a general quality score

=cut

sub compute_quality {
    my ( $self, @args ) = @_;
    my $quality = 0;
    my $cfg = miRkwood->CONFIG();
    my $mode = $cfg->param('job.pipeline');

    if ( $mode eq 'smallRNAseq' ){
        $self->compute_quality_from_reads();
    }
    else{
        if ( $self->{'mfei'} ) {
            if ( $self->{'mfei'} < -0.8 ){
                $quality += 1;
            }
        }
        $quality += $self->{'alignment'};
        $self->{'quality'} =  $quality;
    }

    return;
}

=method compute_quality_for_known_miRNAs

Compute a general quality score for known miRNAs

=cut
sub compute_quality_for_known_miRNAs {
    my ( $self, @args ) = @_;
    my $quality = 0;

    my $precursor_reads = 0;
    my $mature_reads = 0;

    ##### Count number of reads
    foreach (keys %{$self->{'reads'}}){
        $precursor_reads += $self->{'reads'}{$_};
    }
    foreach my $mature_id ( keys %{$self->{'matures'}} ){
        foreach my $read ( keys %{$self->{'matures'}{$mature_id}{'mature_reads'}} ){
            $mature_reads += $self->{'matures'}{$mature_id}{'mature_reads'}{$read};
        }
    }

    $self->{'criteria_nb_reads'} = 0;
    $self->{'criteria_reads_mirna'} = 0;

    ##### Calculate score
    $self->{'quality'} = 0;
    if ( $precursor_reads >= 10 ){
        $self->{'quality'}++;
        $self->{'criteria_nb_reads'} = 1;
    }
    if ( $mature_reads >= ( $precursor_reads / 2 ) ){
        $self->{'quality'}++;
        $self->{'criteria_reads_mirna'} = 1;
    }

    return $self;
}

=method compute_alignment_quality_for_abinitio

Compute the alignment quality score for abinitio pipeline.
To calculate this score :
* 1 point if there is at least one alignment with miRBase
* + 1 point if at least one alignment is validated by miRdup 


=cut

sub compute_alignment_quality_for_abinitio {
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

=method compute_alignment_quality_for_smallRNAseq

Compute the alignment quality score for smallRNAseq pipeline.
To calculate this score :
* 1 point if there is at least one alignment with miRBase
* + 1 point if at least one alignment intersects with at least
40 % of the reads.

=cut

sub compute_alignment_quality_for_smallRNAseq {
    my ( $self, @args ) = @_;
    my $threshold = 0.4;
    my $alignment = 0;
    my $nb_matching_reads = 0;

    my $alignment_existence = 0;
    my $validation = 0;

    if ($self->{'alignment_existence'}){
        $alignment_existence = $self->{'alignment_existence'};
    }else{
        $alignment_existence = 0;
    }

    foreach my $alignment_position (keys%{ $self->{'alignments'} }){
        my ($alignment_start, $alignment_end) = split ( /-/, $alignment_position);
        my $abs_alignment_start = 0;
        my $abs_alignment_end = 0;
        $abs_alignment_start = $self->{'start_position'} + $alignment_start - 1;
        $abs_alignment_end = $self->{'start_position'} + $alignment_end - 1;

        $nb_matching_reads = 0;
        foreach my $read_position (keys%{ $self->{'reads'} }) {
            if ( miRkwood::Utils::is_read_overlapping( "$abs_alignment_start-$abs_alignment_end", $read_position ) ){
                $nb_matching_reads += $self->{'reads'}{$read_position};
            }
        }
        my $ratio = $nb_matching_reads / $self->{'nb_reads'};
        if ( $ratio >= $threshold ){
            $validation = 1;
        }
    }

    $self->{'alignment'} = $alignment_existence + $validation;
    return $self;
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
    if ( exists($self->{'mirna_sequence'}) && $self->{'mirna_sequence'} ne '' ){
        my ($mirna_start, $mirna_end) = split(/-/, $self->{'mirna_position'});
        my $mir_name = 'Mir_' . $self->{'name'} . '__' . $self->{'mirna_position'};
        $text .= $self->{'name'} . "\t" .        # seqid
                    'miRkwood' . "\t" .          # source
                    'miRNA' . "\t" .             # type
                    $mirna_start . "\t" .        # start
                    $mirna_end . "\t" .          # end
                    '.' . "\t" .                 # score
                    $self->{'strand'} . "\t" .   # strand
                    '.' . "\t" .                 # phase
                    "Name=$mir_name" . "\t" .    # attributes
                    "\n";
    }
    return $text;
}

=method candidateAsOrg

Convert a candidate to ORG-mode format

=cut

sub candidateAsOrg {
    my ( $self, @args ) = @_;
    my $cfg = miRkwood->CONFIG();
    my $boolean = { 1 => 'yes', 0 => 'no' };
    my $known_miRNA = 0;
    my $mirna_type = 'novel_miRNA';
    my $output = '';
    my $position = miRkwood::Utils::make_numbers_more_readable( $self->{'position'} );
    $output .= "* Results for $self->{'name'}: $position ($self->{'strand'})\n";
    if ( defined( $self->{'mirbase_id'} ) ){
        $mirna_type = 'known_miRNA';
        $known_miRNA = 1;
        my $mirbase_link = miRkwood::Utils::make_mirbase_link( $self->{'mirbase_id'} );
        $output .= "miRBase name: [[$mirbase_link][$self->{'identifier'}]]\n";
    }
    $output .= "*Chromosome:* $self->{'name'}\n";
    $output .= "*Position:* $position ($self->{'length'} nt)\n";
    $output .= "*Strand:* $self->{'strand'}\n";
    $output .= "*G+C content:* $self->{'%GC'} %\n";
    $output .= "*miRNA sequence:* $self->{'mirna_sequence'}\n";
    if ( ! $known_miRNA ){
        $output .= "*miRNA depth:* $self->{'mirna_depth'} (weigth: $self->{'weight'})\n";
        if ( defined( $self->{'list_id_with_same_mirna'} ) && scalar( @{ $self->{'list_id_with_same_mirna'} } ) ){
            my $list_mirna = '';
            foreach ( @{ $self->{'list_id_with_same_mirna'} } ){
                $list_mirna .= "$_ ";
            }
            $output .= "*Candidates with the same miRNA:* $list_mirna\n";
        }
        else{
            $output .= "*Candidates with the same miRNA:* none\n";
        }
    }
    $output .= "*Stability of the secondary structure of the precursor:* /MFE/ $self->{'mfe'} kcal/mol | /AMFE/ $self->{'amfe'} | /MFEI/ $self->{'mfei'}\n";
    if ( ! $known_miRNA ){
        $output .= "*Stability of the miRNA duplex (mirdup):* $boolean->{ $self->{'criteria_mirdup'} }\n";
    }
    $output .= "*Total number of reads mapped to the precursor:* $self->{'nb_reads'}\n";
    if ( ! $known_miRNA ){
        if ( $self->{'criteria_reads_mirna'} == 1 && $self->{'criteria_star'} == 1 ){
            $output .= "*Distribution of reads:* two islands\n";
        }
        if ( ($self->{'criteria_reads_mirna'} + $self->{'criteria_star'}) == 1 ){
            $output .= "*Distribution of reads:* one island\n";
        }
        if ( $self->{'criteria_reads_mirna'} == 0 && $self->{'criteria_star'} == 0 ){
            $output .= "*Distribution of reads:* random\n";
        } 
    }

    # Read cloud
    my $absolute_read_cloud_path = File::Spec->catfile(
                miRkwood::Paths::get_dir_reads_path_from_job_dir( $cfg->param('job.directory') ),
                $mirna_type,
                $self->{'identifier'}.'.txt');
    $output .= "=\n" . miRkwood::CLI::include_read_cloud_in_html( $absolute_read_cloud_path, $self->{'length'}, $self->{'nb_reads'} ) . "=\n";

    # Alignments
    if ( ! $known_miRNA && $cfg->param('options.align') ){
        if ( $self->{'alignment'} == 0 ){
            $output .= "*miRBase alignment:* none\n";
        }
        else {
            if ( $self->{'alignment'} == 2 ){
                $output .= "*miRBase alignment:* presence of alignments that cover the miRNA locus (see reads cloud above)\n";
            }
            elsif ( $self->{'alignment'} == 1 ){
                $output .= "*miRBase alignment:* presence of alignments, which do not overlap the miRNA locus (see reads cloud above)\n";
            }
            else {
                $output .= 'none\n';
            }
            $output .= "\n" . $self->include_alignments_in_html( 'org' );
        }
    }
    $output .= "\n";
    return $output;
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
            $output .= '>'.$name . "\n" .
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
    my %mirdup_results;
    if ( defined( $self->{'mirdup_validation'} ) ){
        %mirdup_results = %{$self->{'mirdup_validation'}};
    }

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
        my $mirdup_prediction = '';

        # MiRdup
        if ( defined( $self->{'mirdup_validation'} ) ){
            my $mirdup_key = $self->{'name'} . '__' . $position;
            if ( $mirdup_results{$mirdup_key} ){
                $mirdup_prediction = 'This prediction is validated by miRdup.';
            } else {
                $mirdup_prediction = 'This prediction is not validated by miRdup.';
            }
        }

        # Hairpin
        my $hairpin_with_mature = '';
        if ( !eval {
            $hairpin_with_mature =
              miRkwood::Utils::make_hairpin_with_mature($self->{'hairpin'},
                                                        $left, $right,
                                                        length $self->{'sequence'},
                                                        'html');
            }
        ){
            $hairpin_with_mature = $self->{'hairpin'};
        }

        $predictionCounter += 1;

        # Sorting the hit list by descending value of the 'score' element
        my @hits = sort { $b->{'score'} <=> $a->{'score'} } @{$alignments{$position}};
        my $title = "Prediction $predictionCounter: $position";
        $contents .= <<"END_TXT";
        <h3 id='$position'>$title</h3>
        <pre style='height: 80px;'>$hairpin_with_mature</pre>
        $mirdup_prediction
        <h4>Alignments</h4>
END_TXT

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
            $additional_content = "<span class='others'>$mirbase_title" . join(', ', @sequences) . '</span>';


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

    if ( $cfg->param('options.randfold') ) {
        push @fields, ('shuffles');
    }
    if ( $cfg->param('options.align') ) {
        push @fields, ('alignment');
    }

    return @fields;
}

=method get_basic_informations

Return a hash with only the features from the candidate that are 
needed to create or display the results.
/!\ If we add new features to the candidate YML, they will be
automatically added to the basic_candidates YML file unless they
are also added to @elements_to_be_deleted.

=cut
sub get_basic_informations {

	my ( $self, @args ) = @_;

	my @optional_fields = $self->get_optional_candidate_fields();
	my $result = {};

    my @headers = keys %{$self};
    my @elements_to_be_deleted = qw(alignments alignment_existence CT energy_optimal name image mirdup_validation reads cluster chromosome matures);
    foreach (@elements_to_be_deleted){
        @headers = miRkwood::Utils::delete_element_in_array( $_, \@headers);
    }

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
    my ( @args ) = @_;
    my $start_position = shift @args;
    my $structure = shift @args;
    my $start_arm_2 = 0;
    my $end_arm_1 = 0;
    my @structure_array = split(//, $structure);
    my $stop = 0;

    my $index = 0;

    while ( $index < scalar(@structure_array) && ! $stop ){
        if ( $structure_array[$index] eq '(' ){
            $end_arm_1 = $index;
        }
        elsif ( $structure_array[$index] eq ')' ){
            $start_arm_2 = $index;
            $stop = 1;
        }
        $index++;
    }

    $start_arm_2 += $start_position;
    $end_arm_1   += $start_position;

    return ( $end_arm_1, $start_arm_2 );
}

=method compute_quality_from_reads

=cut
sub compute_quality_from_reads {
    my ( $self, @args ) = @_;
    my ($end_arm_1, $start_arm_2) = determine_precursor_arms( $self->{'start_position'}, $self->{'structure_stemloop'});
    my $precursor_length = $self->{'end_position'} - $self->{'start_position'} + 1;
    my $start_arm_1 = $self->{'start_position'};
    my $end_arm_2   = $self->{'end_position'};
    my $reads_arm_1 = {};
    my $reads_arm_2 = {};
    my $count_arm_1 = 0;
    my $count_arm_2 = 0;
    my $criteria_nb_reads = 0;
    my $criteria_star = 0;
    my $criteria_reads_mirna = 0;
    my $criteria_duplex = 0;
    my $criteria_precision_arm_1 = 0;
    my $criteria_precision_arm_2 = 0;
    my $reads_around_mirna = 0;
    my $max_count_arm_1 = 0;
    my $max_count_arm_2 = 0;
    my $max_read_arm_1 = '';
    my $max_read_arm_2 = '';

    # Count the number of reads on each arm.
    # Don't count reads which match (even partially) the loop.
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
        }
    }

    # Criteria nb of reads
    if ( $count_arm_1 >= 10 and $count_arm_2 >= 10 ){
        $criteria_nb_reads = 1;
    }
    if ( $self->{'nb_reads'} >= 100 ){
        $criteria_nb_reads = 1;
    }

    # The two other criteria need the existence of a miRNA
    if ( $self->{'mirna_sequence'} eq '' ){
        $criteria_star = 0;
        $criteria_reads_mirna = 0;
        #~ debug( "             Candidate $self->{'identifier'} : no mirna", miRkwood->DEBUG() );
    }
    else {
        my ($start_mirna, $end_mirna) = split(/-/, $self->{'mirna_position'});
        my $relative_pairing_start_mirna = 0;
        my $relative_pairing_end_mirna = 0;
        my $relative_start_mirna = 0;
        my $relative_end_mirna = 0;
        my $pairing_start_mirna = 0;
        my $pairing_end_mirna = 0;

        $relative_start_mirna = $start_mirna - $self->{'start_position'} + 1;
        if ( $self->{'strand'} eq '-' ){
            $relative_start_mirna = $self->{'length'} + $self->{'start_position'} - $end_mirna;  # vive les maths
        }
        $relative_end_mirna = $relative_start_mirna + $end_mirna - $start_mirna;

        $relative_pairing_start_mirna = $self->find_pairing_position( $relative_start_mirna );
        $relative_pairing_end_mirna = $self->find_pairing_position( $relative_end_mirna );

        if ( $relative_pairing_start_mirna != -1 && $relative_pairing_end_mirna != -1 ){
            if ( $self->{'strand'} eq '+' ){
                $pairing_start_mirna = $relative_pairing_start_mirna + $self->{'start_position'} - 1;
                $pairing_end_mirna = $relative_pairing_end_mirna + $self->{'start_position'} - 1;
            }
            else {
                $pairing_end_mirna = $self->{'length'} + $self->{'start_position'} - $relative_pairing_start_mirna;
                $pairing_start_mirna = $relative_pairing_start_mirna + $pairing_end_mirna - $relative_pairing_end_mirna;
            }
        }
        foreach my $read_position (sort(keys %{$self->{'reads'}})){
            my ($start_read, $end_read) = split(/-/, $read_position);
            # Criteria nb of reads starting in a window [-3; +3] around the
            # mirna start or ending in a window [-5; +5] around the "star" end
            if ( ($start_read >= ($start_mirna - 3) and  $start_read <= ($start_mirna + 3)) ){
                $reads_around_mirna += $self->{'reads'}{$read_position};
            }
            elsif ( $end_read >= ($pairing_start_mirna - 5) and $end_read <= ($pairing_start_mirna + 5) ) {
                $reads_around_mirna += $self->{'reads'}{$read_position};
            }

            # Criteria presence of a star
            if ( ($end_mirna <= $end_arm_1) || ($start_mirna >= $start_arm_2) ){  
                my $overlap = miRkwood::Utils::size_overlap( $read_position, "$pairing_end_mirna-$pairing_start_mirna");
                if ( $overlap > 12 && $self->{'reads'}{$read_position} > 1 ){
                    $criteria_star = 1;
                }
            }

        }
        #~ debug("             $reads_around_mirna reads around the miRNA on a total of $self->{'nb_reads'} (" . ( 100 * $reads_around_mirna / $self->{'nb_reads'}) . ' %)', 1);
        if ( $reads_around_mirna / $self->{'nb_reads'} >= 0.75 ){
            $criteria_reads_mirna = 1;
        }

    }

    $self->{'criteria_nb_reads'} = $criteria_nb_reads;
    $self->{'criteria_star'} = $criteria_star;
    $self->{'criteria_reads_mirna'} = $criteria_reads_mirna;

    # reads distribution takes into account the existence of a star, the
    # precision of the processing (% of reads around the miRNA and/or its
    # star) and the miRdup validation
    if ( ! defined( $self->{'criteria_mirdup'} ) ){
        $self->{'criteria_mirdup'} = 0;
    }
    $self->{'reads_distribution'} = $criteria_star + $criteria_reads_mirna + $self->{'criteria_mirdup'};

    # quality takes into account the reads distribution quality, the
    # number of reads, the existence of a miRNA and the MFEI value
    $self->{'quality'} = $self->{'reads_distribution'};
    $self->{'quality'} += $criteria_nb_reads;
    if ( $self->{'mfei'} < -0.8 ){
        $self->{'quality'} += 1;
    }
    if ( $self->{'mirna_sequence'} ne '' ){
        $self->{'quality'} += 1;
    }

    #~ debug( "             Candidate $self->{'identifier'} ($self->{'quality'}): criteria nb reads : $criteria_nb_reads;  criteria mirdup : $self->{'criteria_mirdup'}; criteria reads around mirna : $criteria_reads_mirna; criteria star : $criteria_star", miRkwood->DEBUG() );

    #~ return $criteria_nb_reads + $criteria_star + $criteria_reads_mirna;
    return $self;
}


sub store_attribute_ct {
    my ( @args ) = @_;
    my $candidate = shift @args;
    my $directory = shift @args;

    if ( -e "$directory/outB2ct_stemloop.ct" ){
        open (my $CT, '<', "$directory/outB2ct_stemloop.ct")
            or warn "Error while opening $directory/outB2ct_stemloop.ct : $!\n";
        my @ct = <$CT>;
        close $CT;
        foreach ( @ct ){
            if ( /(\d+)\s+[a-zA-Z]\s+\d+\s+\d+\s+(\d+)\s+\d+/ ){
                $candidate->{'CT'}{ $1 } = $2;
            }
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

=method find_mirna

  Method to find the majority read, which is considered
  as the more likely miRNA.
  
=cut
sub find_mirna {
    my (@args) = @_;
    my $self = shift @args;
    my $genome_db = shift @args;
    my $max = 0;
    my $position_max = 0;
    my $same_start_reads = 0;
    my $mirna_start = 0;
    my $mirna_end = 0;
    my $threshold = 0.5;
    my $arm = 'both';
    my $count_arm = {};

    # determine the positions of the arms
    my ($end_arm_1, $start_arm_2) = determine_precursor_arms( $self->{'start_position'}, $self->{'structure_stemloop'});
    my $start_arm_1 = $self->{'start_position'};
    my $end_arm_2   = $self->{'end_position'};

    # first browse to find the most abundant read
    # and count the number of reads on each arm
    foreach ( keys%{$self->{'reads'}} ){
        my ($read_start, $read_end) = split( /-/, $_ );
        if ( $self->{'reads'}{$_} > $max ){
            $max = $self->{'reads'}{$_};
            $position_max = $_;
        }
        if ( $read_start >= $start_arm_1 && $read_end <= $end_arm_1 ){
            $count_arm->{ 'first' } += $self->{'reads'}{$_};
        }
        if ( $read_start >= $start_arm_2 && $read_end <= $end_arm_2 ){
            $count_arm->{ 'second' } += $self->{'reads'}{$_};
        }
    }
    $count_arm->{ 'both' } = $self->{'nb_reads'};
    ($mirna_start, $mirna_end) = split( /-/, $position_max );
    if ( $mirna_end <= $end_arm_1 ){
        $arm = 'first';
    }
    elsif ( $mirna_start >= $start_arm_2 ){
        $arm = 'second';
    }

    # second browse to count how many reads start at the same
    # position than the most abundant read
    foreach ( keys%{$self->{'reads'}} ){
        my ($read_start, $read_end) = split( /-/, $_ );
        if ( $read_start eq $mirna_start ){
            $same_start_reads += $self->{'reads'}{$_};
        }
    }

    my $percentage_majority_read = $same_start_reads / $count_arm->{ $arm };

    #~ debug( "             $percentage_majority_read % reads start at the same position than the major read.", 1);
    if ( $percentage_majority_read >= $threshold ){
        my $chromosome = $self->{'name'};
        if ( $self->{'name'} =~ /(.*)__.*/ ){
            $chromosome = $1;
        }
        $self->{'mirna_position'} = "$mirna_start-$mirna_end";
        $self->{'mirna_sequence'} = $genome_db->seq( $chromosome, $mirna_start => $mirna_end );
        $self->{'mirna_sequence'} =~ s/T/U/g;
        if ( $self->{'strand'} eq '-' ){
            $self->{'mirna_sequence'} = miRkwood::Utils::reverse_complement( $self->{'mirna_sequence'} );
        }
        $self->{'mirna_length'}   = $mirna_end - $mirna_start + 1;
        $self->{'mirna_depth'}    = $self->{'reads'}{ $self->{'mirna_position'} };
    }
    else {
        $self->{'mirna_position'} = '';
        $self->{'mirna_sequence'} = '';
        $self->{'mirna_length'}   = '';
        $self->{'mirna_depth'}    = '';
    }

    return $self;
}


=method find_mirna_for_known_candidate
  
=cut

sub find_mirna_for_known_candidate {
    my ($self, @args) = @_;
    my $genome_db = shift @args;
    my $mirna_start = 0;
    my $mirna_end = 0;

    $self->{'mirna_position'} = '';
    $self->{'mirna_sequence'} = '';
    $self->{'mirna_length'}   = '';

    if ( defined( $self->{'matures'} ) and ( scalar(keys%{$self->{'matures'}}) > 0 ) ){
        my @matures_id = keys%{$self->{'matures'}};
        if ( scalar( @matures_id ) == 1 ){
            $mirna_start = $self->{'matures'}{$matures_id[0]}{'mature_start'};
            $mirna_end = $self->{'matures'}{$matures_id[0]}{'mature_end'};
        }
        else {
            my $max = 0;
            my $mirna_max = '';
            foreach my $mature_id ( @matures_id ){
                my $nb_reads = 0;
                foreach ( keys%{ $self->{'matures'}{$mature_id}{'mature_reads'} } ){
                    $nb_reads += $self->{'matures'}{$mature_id}{'mature_reads'}{$_};
                }
                if ( $nb_reads > $max ){
                    $max = $nb_reads;
                    $mirna_max = $mature_id;
                }
            }
            $mirna_start = $self->{'matures'}{$mirna_max}{'mature_start'};
            $mirna_end = $self->{'matures'}{$mirna_max}{'mature_end'};
        }
        $self->{'mirna_position'} = "$mirna_start-$mirna_end";
        $self->{'mirna_sequence'} = $genome_db->seq( $self->{'name'}, $mirna_start => $mirna_end );
        if ( $self->{'strand'} eq '-' ){
            $self->{'mirna_sequence'} = miRkwood::Utils::reverse_complement( $self->{'mirna_sequence'} );
        }
        $self->{'mirna_sequence'} =~ s/T/U/g;
        $self->{'mirna_length'}   = $mirna_end - $mirna_start + 1;
    }

    return $self;

}

sub find_pairing_position {
    my ($self, @args) = @_;
    my $relative_position = shift @args;
    my $pairing_position = 0;

    my $precursor_length = $self->{'end_position'} - $self->{'start_position'} + 1;

    if ( ! defined( $self->{'CT'} ) ){
        warn "No CT available.\n";
        return -1;
    }
    if ( $self->{'CT'}{$relative_position} ne '0' ){
        $pairing_position = $self->{'CT'}{$relative_position};
    }
    else{
        # if the position is not paired (CT value is 0), we look at the closest
        # paired positions surrounding the position, we look for the corresponding
        # paired values and we take as the paired position a value in-between
        # according to a pro-rata.
        my $before = 0;
        my $after = 0;
        while ( $before < $relative_position && $self->{'CT'}{$relative_position - $before} eq '0' ){
            $before++;
        }
        while ( $after < $precursor_length && $self->{'CT'}{$relative_position + $after} eq '0' ){
            $after++;
        }

        $pairing_position = int( $self->{'CT'}{$relative_position - $before} + ( $before / ( $before + $after ) ) * ( $self->{'CT'}{$relative_position + $after} - $self->{'CT'}{$relative_position - $before} ) );
    }

    return $pairing_position;
}



sub count_total_nb_of_reads_for_candidate {
    my (@args) = @_;
    my $self = shift @args;
    my $total_reads = 0;
    if ( scalar(keys %{$self->{'reads'}}) > 0 ){
        foreach my $key (keys( %{$self->{'reads'}} )){
            $total_reads += $self->{'reads'}{$key};
        }
    }
    $self->{'nb_reads'} = $total_reads;
    return $self;
}

=method tag_orphan_hairpins

  Method to set an attribute
  'orphan_hairpin' to 1 if the candidate has a global score
  of 0 and no conservation.

=cut
sub tag_orphan_hairpins{
    my (@args) = @_;
    my $self = shift @args;
    #~ debug("quality: $self->{'quality'}, alignment: $self->{'alignment'}", 1);
    if ( $self->{'quality'} eq '0' && $self->{'alignment'} eq '0'){
        $self->{'orphan_hairpin'} = 1;
    }
    else{
        $self->{'orphan_hairpin'} = 0;
    }
    return $self;
}

=method include_alignments_in_html

=cut
sub include_alignments_in_html {
    my ($self, @args) = @_;
    my $format = shift @args;   # should be 'html' or 'org'
    my $mode_for_hairpin = $format;
    if ( $format eq 'org' ){
        $mode_for_hairpin = 'ascii';
    }
    my $cfg = miRkwood->CONFIG();
    my $alignment_file = File::Spec->catfile(
                miRkwood::Paths::get_dir_alignments_path_from_job_dir( $cfg->param('job.directory') ),
                $self->{'identifier'}.'_aln.txt');
    my $result = '';
    open(my $IN, '<', $alignment_file) or return '';
    my $line;
    while ( $line = <$IN> ){
        chomp $line;
        if ( $line =~ /Prediction : ([\d]+)-([\d]+)/ ){
            if ( $format eq 'org' ){
                $line = "*$line*\n\n";
            }
            else {
                $line = "<b style='font-size:1.2em;'>$line</b>\n\n";
            }
            my $hairpin_with_mature = '';
            if ( !eval {
                $hairpin_with_mature =
                  miRkwood::Utils::make_hairpin_with_mature($self->{'hairpin'},
                                                            $1, $2,
                                                            length $self->{'sequence'},
                                                            $mode_for_hairpin);
                }
            ){
                $hairpin_with_mature = $self->{'hairpin'};
            }
            $line .= $hairpin_with_mature;
            if ( $format eq 'org' ){
                $line .= "\n*Alignments*\n";
            }
            else {
                $line .= "\n<b>Alignments</b>\n";
            }
        }
        if ( $line =~ /miRBase (\d+):(.*)/ ){
            my @list_results = split( '\|', $2);
            my $mirbase_link = "miRBase $1 : ";
            foreach my $result ( @list_results ){
                if ( $result =~ /^ ([^ ]*) / ){
                    if ( $format eq 'org' ){
                        $mirbase_link .= "[[miRkwood::Utils::make_mirbase_link( $1 )][$1]] | ";
                    }
                    else {
                        $mirbase_link .= '<a href="' . miRkwood::Utils::make_mirbase_link( $1 ) . '">' . $1 . '</a> | ';
                    }
                }
            }
            $line = $mirbase_link;
        }
        $result .= "$line\n";
    }
    close $IN;
    return $result;
}


1;
