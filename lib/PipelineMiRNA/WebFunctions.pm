package PipelineMiRNA::WebFunctions;

# ABSTRACT: Code directly tied to the web interface

use strict;
use warnings;

use Data::Dumper;
use File::Spec;
use Time::gmtime;
use Switch;

use PipelineMiRNA;
use PipelineMiRNA::Paths;
use PipelineMiRNA::Parsers;
use PipelineMiRNA::WebTemplate;
use PipelineMiRNA::Components;


=method get_optional_candidate_fields

Return the optional fields based on the current configuration

=cut

sub get_optional_candidate_fields {
    my ( $self, @args ) = @_;
    my @fields = ();
    my $cfg = PipelineMiRNA->CONFIG();

    if ($cfg->param('options.mfe')){
        push @fields, ('mfe', 'mfei', 'amfe');
    }
    if ($cfg->param('options.randfold')){
        push @fields, ('p_value');
    }
    if ($cfg->param('options.align')){
        push @fields, ('alignment')
    }
    return @fields;
}

=method make_job_id

Return a jobId (based on the current time)

=cut

sub make_job_id {
    my ( $self, @args ) = @_;
    my $now = gmctime();
    $now =~ s/[: ]//g;
    $now = substr( $now, 3 );
    return $now;
}

=method jobId_to_jobPath

Get the job path from a job identifier

=cut

sub jobId_to_jobPath {
    my ( $self, @args ) = @_;
    my $id_job      = shift @args;
    my $dirJob_name = 'job' . $id_job;
    my $results_dir = PipelineMiRNA::Paths->get_results_dir_name();
    my $jobPath = File::Spec->catdir( $results_dir, $dirJob_name);
    return $jobPath;
}

=method is_valid_jobID

Test whether a jobID is valid - ie if there are results for it.

=cut

sub is_valid_jobID {
    my ( $self, @args ) = @_;
    my $id_job          = shift @args;
    my $jobPath = $self->jobId_to_jobPath($id_job);
    my $full_path = PipelineMiRNA::Paths->get_absolute_path($jobPath);
    return (-e $full_path);
}

=method get_structure_for_jobID

Get the results structure of a given job identifier

Usage:
my %results = PipelineMiRNA::WebFunctions->get_structure_for_jobID($jobId);

=cut

sub get_structure_for_jobID {
    my ( $self, @args ) = @_;
    my $jobId   = shift @args;
    my $job = $self->jobId_to_jobPath($jobId);
    my $job_dir = PipelineMiRNA::Paths->get_absolute_path($job);
    PipelineMiRNA->CONFIG_FILE(PipelineMiRNA::Paths->get_job_config_path($job_dir));
    my %myResults = ();

    opendir DIR, $job_dir;    #ouverture répertoire job
    my @dirs;
    @dirs = readdir DIR;
    closedir DIR;
    foreach my $dir (@dirs)    # parcours du contenu
    {
        my $full_dir = File::Spec->catdir( $job_dir, $dir );
        if (    $dir ne "."
             && $dir ne ".."
             && -d $full_dir )    #si fichier est un répertoire
        {
            opendir DIR, $full_dir;    # ouverture du sous répertoire
            my @files;
            @files = readdir DIR;
            closedir DIR;
            foreach my $subDir (@files) {
                my $subDir_full = File::Spec->catdir( $job_dir, $dir, $subDir );
                if (    ( $subDir ne "." )
                     && ( $subDir ne ".." )
                     && -d $subDir_full ) # si le fichier est de type repertoire
                {
                    my %candidate = $self->retrieve_candidate_information($job, $dir, $subDir);
                    $myResults{$subDir} = \%candidate;
                }
            }
        }
    }
    return %myResults;
}

=method retrieve_candidate_information

Check correctness and get the result for a given candidate

Arguments:
- $job - the job identifier
- $dir - the sequence name
- $subDir - the candidate name

=cut

sub retrieve_candidate_information {
    my ( $self, @args ) = @_;
    my $job = shift @args;
    my $dir = shift @args;
    my $subDir = shift @args;

    my ($candidate_dir, $full_candidate_dir) = PipelineMiRNA::Paths->get_candidate_paths($job,  $dir, $subDir);

    if ( ! -e $full_candidate_dir ){
        die('Unvalid candidate information');

    }else{
        my %result = $self->actual_retrieve_candidate_information($candidate_dir, $full_candidate_dir);
        $result{'name'} = $dir;    #récupération nom séquence
        my @position = split( /__/, $subDir );
        $result{'position'} = $position[1]; # récupération position
        return %result;
    }
}

=method actual_retrieve_candidate_information

Get the results for a given candidate

Arguments:
- $candidate_dir - the unprefixed path to the candidate results
- $full_candidate_dir - the prefixed path to the candidate results

=cut

sub actual_retrieve_candidate_information {
    my ( $self, @args ) = @_;
    my $candidate_dir = shift @args;
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

    #Récupération valeur self contain
    my $selfcontain_out =
      File::Spec->catfile( $full_candidate_dir, 'selfContain.txt' );
    if ( -e $selfcontain_out )
    {Dumper
        $result{'self_contain'} = PipelineMiRNA::Parsers::parse_selfcontain($selfcontain_out);
    }

    #Récupération séquence et format Vienna
    my $rnafold_stemloop_out = File::Spec->catfile( $full_candidate_dir,
                                       'outRNAFold_stemloop.txt' );
    if ( -e $rnafold_stemloop_out )                  # si fichier existe
    {
        my @vienna_res = PipelineMiRNA::Parsers::parse_RNAfold_output($rnafold_stemloop_out);


        $result{'DNASequence'} = $vienna_res[1];
        $result{'Vienna'} = $vienna_res[2];
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
    $result{'alignment'} = ( -e $file_alignement && ! -z $file_alignement );

    my $image_path = File::Spec->catfile($candidate_dir, 'image.png');
    $result{'image'} = $image_path;

    # Computing general quality
    $result{'quality'} = $self->compute_quality(\%result);

    return %result;
}

=method compute_quality

Compute a general quality score

=cut

sub compute_quality(){
    my ( $self, @args ) = @_;
    my %result = %{shift @args};
    my $quality = 0;
    if ( $result{'mfei'} < -0.5 ){
        $quality += 1;
    }
    my $length = length ($result{'DNASequence'});

    if ( $length > 80 && $length < 200 ){
        $quality += 1;
    }
    return $quality;
}

=method candidateAsVienna

Convert a given candidate to Vienna dot-bracket format

=cut

sub candidateAsVienna {
    my ( $self, @args ) = @_;
    my %candidate = %{shift @args};
    my $optimal = shift @args;
    my $output = "";
    my $candidate_name = $candidate{'name'}.'__'.$candidate{'position'};
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
    $output .= '>'.$candidate{'name'} . '__' . $candidate{'position'} . "\n" . $candidate{'DNASequence'} . "\n";
    return $output;
}


=method export

Convert the results

=cut

sub export {
    my ( $self, @args ) = @_;
    my $export_type = shift @args;
    my $results = shift @args;
    my @sequences_to_export = shift @args;
    my %results = %{$results};
    my $output = "";

    my @keys = sort keys %results;
    foreach my $key(@keys)
    {
        if (  $key ~~ \@sequences_to_export )
        {
            my $value = $results{$key};
            switch ($export_type) {
                case 'fas'        { $output .= $self->candidateAsFasta($value); }
                case 'dot'        { $output .= $self->candidateAsVienna($value); }
            }
        }
    }
    return $output;
}

=method resultstruct2csv

Convert the results structure to CSV

=cut

sub resultstruct2csv {
    my ( $self, @args ) = @_;
    my $results = shift @args;
	my @tab = shift @args;
	my %results = %{$results};
	my @optional_fields = $self->get_optional_candidate_fields();
    my @csv_headers = ('name', 'position', @optional_fields, 'Vienna', 'DNASequence');
	my $result = join( ',', @csv_headers ) . "\n";

    my @keys = sort keys %results;
    foreach my $key(@keys)
    {
    	if (  $key ~~ \@tab ) 
    	{
    	    my $value = $results{$key};
	        for my $header (@csv_headers)
 	        {
	            $result .= "${$value}{$header},";
	        }
	        $result .= "\n";
    	}
    }
    return $result;
}

=method resultstruct2pseudoXML

Convert the results structure to to pseudo XML format

=cut

sub resultstruct2pseudoXML {
	
    my ( $self, @args ) = @_;
    my $results = shift @args;
    my %results = %{$results};

    my @optional_fields = $self->get_optional_candidate_fields();
    my @headers1 = ('name', 'position', 'quality', @optional_fields);
    my @headers2 = ('Vienna', 'DNASequence');

    my $result = "<results id='all'>\n";
	my @keys = sort keys %results;

    foreach my $key (@keys) {
   		my $value = $results{$key};
       	$result .= "<Sequence";
        for my $header (@headers1){
            $result .= " $header='${$value}{$header}'";
        }
        my $img = PipelineMiRNA::Paths->get_server_path(${$value}{'image'});
        $result .= " image='$img'";
        for my $header (@headers2){
            $result .= " $header='${$value}{$header}'";
        }
        $result .= "></Sequence>\n";
    }
    $result .= "</results>";
    $result .= "<results id='all2'>\n";
   	@keys = sort keys %results;
    @keys = sort {$results{$b}{'quality'} <=> $results{$a}{'quality'}||$results{$a}{'position'} <=> $results{$b}{'position'} } keys %results; 
  	foreach my $key (@keys) {
   		my $value = $results{$key};
       	$result .= "<Sequence";
        for my $header (@headers1){
            $result .= " $header='${$value}{$header}'";
        }
        my $img = PipelineMiRNA::Paths->get_server_path(${$value}{'image'});
        $result .= " image='$img'";
        for my $header (@headers2){
            $result .= " $header='${$value}{$header}'";
        }
        $result .= "></Sequence>\n";
    }
    $result .= "</results>";
    return $result;
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

sub get_element_of_split {
    my @args = @_;
    my $value = shift @args;
    my $rank  = shift @args;
    my @split = split(/-/, $value);
    return $split[$rank];
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

    my @keys = sort { get_element_of_split($a, 0)  <=> get_element_of_split($b, 0) || get_element_of_split($a, 1)  <=> get_element_of_split($b, 1)} keys %alignments;
    my @stocked_hits;
    my $final_key;
    my $final_hit_count = -1;
    foreach my $current_key (@keys) {
        my ($current_left, $current_right) = split(/-/, $current_key);
        my $current_hit_count = scalar @{$alignments{$current_key}};

        if ( ! ($current_right <= $stocked_right + 2 )) {

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
    my $job = shift @args;
    my $dir = shift @args;
    my $subDir = shift @args;
    my $hairpin = shift @args;
    my ($candidate_dir, $full_candidate_dir) = PipelineMiRNA::Paths->get_candidate_paths($job,  $dir, $subDir);
    my $file_alignement = File::Spec->catfile($full_candidate_dir, 'alignement.txt');

    my %alignments;
    %alignments = PipelineMiRNA::Components::parse_custom_exonerate_output($file_alignement);
    if (! eval {%alignments = PipelineMiRNA::Components::parse_custom_exonerate_output($file_alignement);}) {
        # Catching exception
        my $error = "Error with alignments $file_alignement.";
        print PipelineMiRNA::WebTemplate::get_error_page($error);
        die($error);
    }

    my $contents = "";
    my @TOC;
    my $predictionCounter = 0;

    %alignments = $self->merge_alignments(\%alignments);

    my @keys = sort { get_element_of_split($a, 0)  <=> get_element_of_split($b, 0) || get_element_of_split($a, 1)  <=> get_element_of_split($b, 1)} keys %alignments;

    foreach my $position (@keys) {
        my ($left, $right) = split(/-/, $position);
        my ($top, $upper, $middle, $lower, $bottom) = split(/\n/, $hairpin);

        my $hairpin_with_mature;

        if ($left > length $top)
        {
            #on the other side
            $hairpin_with_mature = $hairpin;
        } else {
            my $size = PipelineMiRNA::Utils::compute_mature_boundaries($left, $right, $top);
            substr($top, $left, $size)   = '<span class="mature">' . substr($top, $left, $size) . '</span>';
            substr($upper, $left, $size) = '<span class="mature">' . substr($upper, $left, $size) . '</span>';
            $hairpin_with_mature = <<"END";
$top
$upper
$middle
$lower
$bottom
END
        }
        $predictionCounter += 1;
        # Sorting the hit list by descending value of the 'score' element
        my @hits = sort { $b->{'score'} <=> $a->{'score'} } @{$alignments{$position}};
        my $title = "Prediction $predictionCounter: $position";
        $contents .= "<h3 id='$position'>$title</h3>
        <pre style='height: 80px;'>$hairpin_with_mature</pre>
        <ul>
            <li>Evaluation score of MIRdup: TODO</li>
        </ul>
        <h4>Alignments</h4>
        ";

        my $toc_element = "<a href='#$position'>$position</a>";
        push @TOC, $toc_element;
        foreach my $hit (@hits){
            my $alignment = $hit->{'alignment'};
            my $name = $hit->{'name'};
            my @splitted = split(/ /, $hit->{'def_query'});
            my $mirbase_id = $splitted[0];
            my $mirbase_link = PipelineMiRNA::WebTemplate::make_mirbase_link($mirbase_id);
            my $html_name = "<a href='$mirbase_link'>$name</a>";
            my $spacing = 15;
            my ($top, $middle, $bottom) = split(/\n/, $alignment);
            $top    = sprintf "%-${spacing}s %3s %s %s", 'query', $hit->{'begin_target'}, $top,   $hit->{'end_target'};
            $middle = sprintf "%-${spacing}s %3s %s %s", '',      '',                     $middle, '';
            $bottom = sprintf "%-${spacing}s %3s %s %s", $name,   $hit->{'begin_query'},  $bottom, $hit->{'end_query'};
            my $additional_space = "";
            my $sub_string = substr($bottom, 0, $spacing);
            $additional_space .= ' ' while ($sub_string =~ m/ /g);
            substr($bottom, 0, $spacing) = $html_name . $additional_space;
            $contents .= <<"INNER";
<pre>
$top
$middle
$bottom
</pre>
INNER
        }

    }
    my $toc = "<span class='toc'>" . join(" - ", @TOC) . '</span>';
    return $toc . "\n" . $contents;

}

1;
