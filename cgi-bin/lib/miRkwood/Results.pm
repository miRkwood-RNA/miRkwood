package miRkwood::Results;

# ABSTRACT: Code directly tied to the results data structure

use strict;
use warnings;

use feature 'switch';
use Time::gmtime;
use File::Spec;
use YAML::XS;

use miRkwood;
use miRkwood::Paths;
use miRkwood::Candidate;
use miRkwood::CandidateHandler;;
use miRkwood::Utils;

=method make_job_id

Return a jobId (based on the current time)

=cut

sub make_job_id {
	my ( $self, @args ) = @_;
	my $type = shift @args;     # type should be 'BAM' or 'Fasta'
	my $now = gmctime();
	$now =~ s/[: ]//g;
	$now = $type . substr( $now, 3 );
	return $now;
}

=method jobId_to_jobPath

Get the job path from a job identifier

=cut

sub jobId_to_jobPath {
	my ( $self, @args ) = @_;
	my $id_job      = shift @args;
	my $dirJob_name = 'job' . $id_job;
	my $results_dir = miRkwood::Paths->get_results_filesystem_path();
	my $jobPath     = File::Spec->catdir( $results_dir, $dirJob_name );
	return $jobPath;
}

=method is_job_finished

Return whether a job is finished or not

=cut

sub is_job_finished {
    my ( $self, @args ) = @_;
    my $id_job      = shift @args;
    my $job_dir     = $self->jobId_to_jobPath($id_job);
    my $is_finished_file = File::Spec->catfile( $job_dir, 'finished' );
    return (-e $is_finished_file);
}

=method is_valid_jobID

Test whether a jobID is valid - ie if there are results for it.

=cut

sub is_valid_jobID {
	my ( $self, @args ) = @_;
	my $id_job    = shift @args;
	my $full_path = $self->jobId_to_jobPath($id_job);
	return ( -e $full_path );
}

=method get_structure_for_jobID

Get the results structure of a given job identifier

Usage:
my %results = miRkwood::Results->get_structure_for_jobID($jobId, $type);
$type should be 'known_miRNA', 'novel_miRNA' (for BEDPipeline) or '' (for other pipelines)

=cut

sub get_structure_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
	my $mirna_type = shift @args;   # should be 'known_miRNA' or anything else
	my $job_dir = $self->jobId_to_jobPath($jobId);
	miRkwood->CONFIG_FILE(
		miRkwood::Paths->get_job_config_path($job_dir) );
	my $candidates_file = '';
    if ( $mirna_type eq 'known_miRNA' ){
        $candidates_file = File::Spec->catfile( $job_dir, 'basic_known_candidates.yml');
    }
    else{
        $candidates_file = File::Spec->catfile( $job_dir, 'basic_candidates.yml');
    }
	return $self->deserialize_results($candidates_file);
}

=method get_basic_structure_for_jobID

=cut

sub get_basic_structure_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
    my $type    = shift @args;  # $type should be 'novel_miRNA' or 'known_miRNA'
	my $job_dir = $self->jobId_to_jobPath($jobId);
    my $yml_file = '';
    if ( $type eq 'known_miRNA' ){
        $yml_file = 'basic_known_candidates.yml';
    }
    else{
        $yml_file = 'basic_candidates.yml';
    }
	miRkwood->CONFIG_FILE(
		miRkwood::Paths->get_job_config_path($job_dir) );
	my $candidates_file = File::Spec->catfile( $job_dir, $yml_file);
	return miRkwood::get_yaml_file( $candidates_file );
}

sub get_basic_pseudoXML_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
    my $pipeline_type = shift @args;    # should be 'abinitio' or 'smallRNAseq'
    my $type    = shift @args;  # $type should be 'novel_miRNA' or 'known_miRNA'

	my $results = $self->get_basic_structure_for_jobID($jobId, $type);

	my $output = '';

    $output .= "<results id='all'>\n";
    my @candidates = sort {
        ( $a->{'name'} cmp  $b->{'name'} )
          || (
            $a->{'start_position'} <=> $b->{'start_position'} )
    } @{$results};

    foreach my $candidate (@candidates) {
        $output .= $self->convert_basic_to_pseudoXML($candidate, $pipeline_type, $type) . "\n";
    }
    $output .= "</results>\n";

    $output .= "<results id='all2'>\n";
    @candidates = sort {
        ( $b->{'quality'} cmp $a->{'quality'} )
          ||
        ( $a->{'name'} cmp  $b->{'name'} )
          || (
            $a->{'start_position'} <=> $b->{'start_position'} )
    } @candidates;
    foreach my $candidate (@candidates) {
        $output .= $self->convert_basic_to_pseudoXML($candidate, $pipeline_type, $type) . "\n";
    }
    $output .= '</results>';

    return $output;
}

sub convert_basic_to_pseudoXML {
	my ( $self, @args ) = @_;
	my $candidate = shift @args;
	my %candidate = %{$candidate};
    my $pipeline_type = shift @args;    # should be 'abinitio' or 'smallRNAseq'
    my $type      = shift @args;  # $type should be 'novel_miRNA' or 'known_miRNA'

	my @headers;
    my @fields_to_truncate = qw{mfe mfei amfe};
    my @optional_candidate_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @optional_mirna_fields = miRkwood::Candidate->get_optional_mirna_fields();

    if ( $pipeline_type eq 'smallRNAseq' ){
        if ( $type eq 'known_miRNA' ){  # known miRNAs for pipeline smallRNAseq
            @headers = qw{name precursor_name position strand mirna_sequence mirna_length quality nb_reads identifier criteria_nb_reads};
        }
        else {  # new miRNAs for pipeline smallRNAseq
            push @headers, ( 'name', 'position', 'strand', 'mirna_sequence', 'mirna_length', 'mirna_depth', 'weight', 'reads_distribution', 'mfei', 'nb_reads', @optional_candidate_fields, @optional_mirna_fields, 'identifier', 'criteria_nb_reads' );
        }
    }
    else {  # pipeline ab initio
       push @headers, ( 'name', 'position', 'length', 'strand', 'quality', 'mfe', 'mfei', 'amfe', @optional_candidate_fields, @optional_mirna_fields, 'image', 'identifier' );
    }

    my $result = '<Sequence';
    for my $header (@headers) {
        my $contents = $candidate->{$header};
        if (grep { $header eq $_ } @fields_to_truncate){
            $contents = miRkwood::Utils::restrict_num_decimal_digits($contents, 3);
        }
        if ( !defined $contents ) {
            $contents = q{};
        }
        if ( $header eq 'shuffles' && $contents == 1){
            $contents = q{};
        }
        if ( $header eq 'position' ){
            my ($start, $end) = split( /-/, $contents);
            $start = miRkwood::Utils::make_numbers_more_readable( $start );
            $end = miRkwood::Utils::make_numbers_more_readable( $end );
            $contents = "$start-$end";
        }
        $result .= " $header='$contents'";
    }
    $result .= '></Sequence>';

    return $result;
}

=method has_candidates

Parse and serialize the results structure of $job_dir

Usage:
miRkwood::Results->has_candidates( \%myResults );

=cut

sub has_candidates {
	my ( $self, @args ) = @_;
	my $results = shift @args;
	my %results = %{$results};
	return ( keys %results > 0 );
}

=method deserialize_results

Retrieve the results in the given basic yml file

Usage:
my %results = miRkwood::Results->deserialize_results( $basic_candidates_file );

=cut

sub deserialize_results {
	my ( $self, @args ) = @_;
	my $basic_candidates_file = shift @args;
    my %myResults = ();
    my $results = miRkwood::get_yaml_file( $basic_candidates_file );
    foreach my $element ( @{$results} ){
        $myResults{ $element->{'identifier'} } = miRkwood::Candidate->new( $element );
    }
    return %myResults;
}

=method number_of_results

return total number of candidates 

=cut

sub number_of_results {
	my ( $self, @args ) = @_;
	my $results = shift @args;
	my %results = %{$results};
	my $size    = scalar keys %results;
	return $size;
}

sub number_of_results_bis {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
    my $type    = shift @args;  # $type should be 'novel_miRNA' or 'known_miRNA'
	my $results = $self->get_basic_structure_for_jobID($jobId, $type);
	my $size    = scalar @{$results};
	return $size;
}

sub count_reads_in_basic_yaml_file {
    my ( @args ) = @_;
    my $yaml = shift @args;
    my $nb_reads = 0;

    my @attributes = YAML::XS::LoadFile($yaml);

    foreach my $data ( @{$attributes[0]} ){
        $nb_reads += $data->{'nb_reads'};
    }

    return $nb_reads;
}

=method count_candidates_per_quality

Method to count the number of novel candidates with quality 0, 1, ...
Input : the basic_candidates.yml file
Output: a hash with key=quality and value=nb of novel candidates

=cut
sub count_candidates_per_quality {
    my ( @args ) = @_;
    my $yaml = shift @args;
    my $quality_count = {};

    my @attributes = YAML::XS::LoadFile($yaml);

    foreach my $data ( @{$attributes[0]} ){
        if ( ! defined( $quality_count->{ $data->{'quality'} } ) ){
            $quality_count->{ $data->{'quality'} } = 0;
        }
        $quality_count->{ $data->{'quality'} } += 1;
    }

    return $quality_count;
}

=method count_candidates_per_reads_distribution

Method to count the number of novel candidates with reads_distribution 0, 1, ...
Input : the basic_candidates.yml file
Output: a hash with key=reads_distribution and value=nb of novel candidates

=cut
sub count_candidates_per_reads_distribution {
    my ( @args ) = @_;
    my $yaml = shift @args;
    my $reads_distrib_count = {};

    my @attributes = YAML::XS::LoadFile($yaml);

    foreach my $data ( @{$attributes[0]} ){
        if ( ! defined( $reads_distrib_count->{ $data->{'reads_distribution'} } ) ){
            $reads_distrib_count->{ $data->{'reads_distribution'} } = 0;
        }
        $reads_distrib_count->{ $data->{'reads_distribution'} } += 1;
    }

    return $reads_distrib_count;
}

sub make_reads_barchart {
    my ( $self,
         $total_witdh,
         $percentage_CDS_reads,
         $percentage_other_reads,
         $percentage_multi_reads,
         $percentage_known_miRNAs_reads,
         $percentage_new_miRNAs_reads,
         $percentage_orphan_clusters_reads,
         $percentage_orphan_hairpins_reads ) = @_;

    my $width;
    $width->{'CDS'}             = int($percentage_CDS_reads * $total_witdh / 100 + 0.5 );
    $width->{'other'}           = int($percentage_other_reads * $total_witdh / 100 + 0.5 );
    $width->{'multi'}           = int($percentage_multi_reads * $total_witdh / 100 + 0.5 );
    $width->{'orphan_clusters'} = int($percentage_orphan_clusters_reads * $total_witdh / 100 + 0.5 );
    $width->{'orphan_hairpins'} = int($percentage_orphan_hairpins_reads * $total_witdh / 100 + 0.5 );
    $width->{'known_miRNAs'}    = int($percentage_known_miRNAs_reads * $total_witdh / 100 + 0.5 );
    $width->{'new_miRNAs'}      = int($percentage_new_miRNAs_reads * $total_witdh / 100 + 0.5 );
    $width->{'unclassified_reads'}         = $total_witdh
                                        - $width->{'CDS'}
                                        - $width->{'other'}
                                        - $width->{'multi'} 
                                        - $width->{'known_miRNAs'}
                                        - $width->{'new_miRNAs'}
                                        - $width->{'orphan_clusters'}
                                        - $width->{'orphan_hairpins'};

    my $barchart = "<div style='width:${total_witdh}px'>\n";
    $barchart .= "<table id='barchart_table'>\n";
    $barchart .= "<tr>\n";
    my @categories = qw{CDS other multimapped orphan_clusters orphan_hairpins unclassified_reads known_miRNAs new_miRNAs};
    foreach my $category ( @categories ){
        if ( $width->{$category} > 0){
            $barchart .= "<td id='$category' style='width:$width->{$category}px'></td>\n";
        }
    }
    $barchart .= "</tr>\n";
    $barchart .= "</table>\n";
    $barchart .= "</div>\n";

    return $barchart;

}

sub create_reads_archive {
 	my ( $self, @args ) = @_;
	my $job_id   = shift @args;
    my $mirna_type = shift @args;
    my $id_results = shift @args;
    my $job_dir = miRkwood::Results->jobId_to_jobPath($job_id);
    my $reads_path = miRkwood::Paths::get_dir_reads_path_from_job_dir_and_mirna_type( $job_dir, $mirna_type );
    my $list_reads_files = '';
    foreach my $id ( @{$id_results} ){
        $list_reads_files .= $id.'.txt ';
    }
    my $archive_path = "$job_dir/reads.tar.gz";
    my $cmd = "tar zcf $archive_path -C $reads_path $list_reads_files";
    system($cmd);
    return $archive_path;
}

=method

  Create a page with a summary of options
  and results (number of reads...)

=cut
sub create_summary_page {
    my (@args) = @_;
    my $absolute_job_dir = shift @args;
    my $output_file = shift @args;
    my $run_options_file = miRkwood::Paths->get_job_config_path($absolute_job_dir);
    miRkwood->CONFIG_FILE($run_options_file);
    my $cfg = miRkwood->CONFIG();
    my $annotation_gff = $cfg->param( 'options.annotation_gff' );
    my @annotation_gff = split( /\&/, $annotation_gff );

    my %boolean_mapping = ( 0 => 'No', 1 => 'Yes', '' => 'No' );
    my $bed_sizes;
    my $basename_bed = $cfg->param('job.bed');
    my $basic_known_yaml = File::Spec->catfile( $absolute_job_dir, 'basic_known_candidates.yml');
    my $basic_yaml = File::Spec->catfile( $absolute_job_dir, 'basic_candidates.yml');
    my $nb_reads_known_miRNAs = miRkwood::Results::count_reads_in_basic_yaml_file( $basic_known_yaml );
    my $nb_reads_new_miRNAs = miRkwood::Results::count_reads_in_basic_yaml_file( $basic_yaml );
    #~ my $nb_candidates_per_quality = miRkwood::Results::count_candidates_per_quality( $basic_yaml );
    #~ my $nb_candidates_per_reads_distrib = miRkwood::Results::count_candidates_per_reads_distribution( $basic_yaml );

    my %known_results = miRkwood::Results->deserialize_results($basic_known_yaml);
    my %novel_results = miRkwood::Results->deserialize_results($basic_yaml);

    my $nb_results_known_miRNAs = miRkwood::Results->number_of_results( \%known_results );
    my $nb_results_new_miRNAs = miRkwood::Results->number_of_results( \%novel_results );

    my $bed_sizes_file = File::Spec->catfile( $absolute_job_dir, miRkwood::Paths::get_bed_size_file_name() );
    open (my $FH, '<', $bed_sizes_file) or die "ERROR while opening $bed_sizes_file : $!";
    while ( <$FH> ){
        if ( $_ !~ /^#/ ){
            chomp;
            my @line = split(/\t/);
            my $name = '';
            if ( $line[0] eq "$basename_bed.bed" ){
                $name = $basename_bed;
            }
            elsif ( $line[0] =~ /${basename_bed}_(.*)[.]bed/ ){
                $name = $1;
            }
            $bed_sizes->{$name}{'reads'} = $line[1];
            $bed_sizes->{$name}{'unique_reads'} = $line[2];
        }
    }
    close $FH;

    my $nb_orphan_reads = $bed_sizes->{$basename_bed}{'reads'} - $nb_reads_known_miRNAs - $nb_reads_new_miRNAs;
    foreach my $category ( keys%{$bed_sizes} ){
        if ( $category ne $basename_bed && $category ne 'miRNAs' ){
            $nb_orphan_reads -= $bed_sizes->{ $category }{'reads'};
        }
    }

    my $multimapped_interval = $cfg->param('options.multimapped_interval');
    my $min_nb_positions = 0;
    my $max_nb_positions = 5;
    if ( $multimapped_interval =~ /\[(.*);(.*)\]/ ){
        $min_nb_positions = $1;
        $max_nb_positions = $2;
    }

    my $output_txt = '';
    $output_txt .= "OPTIONS SUMMARY\n";
    $output_txt .= "===============\n\n";
    $output_txt .= 'BED file: ' . $basename_bed . ".bed\n";
    $output_txt .= 'Reference species: ' . $cfg->param('job.plant') . "\n";
    $output_txt .= 'Flag conserved mature miRNAs: ' . $boolean_mapping{ $cfg->param('options.align') } . "\n";
    $output_txt .= 'Select only sequences with MFEI < -0.6: ' . $boolean_mapping{ $cfg->param('options.mfei') } . "\n";
    $output_txt .= 'Compute thermodynamic stability: ' . $boolean_mapping{ $cfg->param('options.randfold') } . "\n";
    foreach my $gff ( @annotation_gff ){
        if ( $gff =~ /[\/\\]([^\/\\]+[.](dat|gtf|gff3?))/ ){
            $output_txt .= "Filter out features given in $1\n";
        }
    }
    if ( $cfg->param('options.mirbase_gff') =~ /[\/\\]([^\/\\]+[.](dat|gtf|gff3?))/ ){
        $output_txt .= "Filter out known miRNAs present in $1\n";
    }
    if ( $min_nb_positions == 0 && $max_nb_positions == 0){
        $output_txt .= "Filter multiple mapped reads: No\n";
    }
    else{
        $output_txt .= "Filter multiple mapped reads: keep reads mapping at $min_nb_positions to $max_nb_positions positions\n";
    }
    $output_txt .= 'Filter low quality hairpins: ' . $boolean_mapping{ $cfg->param('options.filter_bad_hairpins') } . "\n";

    $output_txt .= "\nRESULTS SUMMARY\n";
    $output_txt .= "===============\n\n";
    $output_txt .= "Total number of reads: $bed_sizes->{$basename_bed}{'reads'} ($bed_sizes->{$basename_bed}{'unique_reads'} unique reads)\n";
    foreach my $gff ( @annotation_gff ){
        if ( $gff =~ /[\/\\]([^\/\\]+)_(.*)[.](dat|gtf|gff3?)/ ){
            $output_txt .= "${1}_$2: $bed_sizes->{ $2 }{'reads'} reads\n";
        }
    }
    $output_txt .= "Multiple mapped reads: $bed_sizes->{'multimapped'}{'reads'} reads\n";
    $output_txt .= "Orphan clusters of reads: $bed_sizes->{'orphan_clusters'}{'reads'} reads\n";
    if ( $cfg->param('options.filter_bad_hairpins') ){
        $output_txt .= "Orphan hairpins: $bed_sizes->{'orphan_hairpins'}{'reads'} reads\n";
    }
    $output_txt .= "Unclassified reads: $nb_orphan_reads reads\n";
    $output_txt .= "Known miRNAs: $nb_results_known_miRNAs sequence(s) - $nb_reads_known_miRNAs reads\n";
    $output_txt .= "Novel miRNAs: $nb_results_new_miRNAs sequence(s) - $nb_reads_new_miRNAs reads\n";
    #~ $output_txt .= "\nDistribution of novel miRNAs according to quality:\n";
    #~ foreach my $qual (reverse(sort(keys%{$nb_candidates_per_quality}))){
        #~ $output_txt .= "    quality $qual: $nb_candidates_per_quality->{$qual} miRNAs\n";
    #~ }
    #~ $output_txt .= "\nDistribution of novel miRNAs according to reads distribution:\n";
    #~ foreach my $qual (reverse(sort(keys%{$nb_candidates_per_reads_distrib}))){
        #~ $output_txt .= "    reads distribution quality $qual: $nb_candidates_per_reads_distrib->{$qual} miRNAs\n";
    #~ }

    open (my $OUT, '>', $output_file) or die "ERROR while creating $output_file : $!";
    print $OUT $output_txt;
    close $OUT;
    return;
}

sub clean_job_dir_for_cli_pipeline {
    my (@args) = @_;
    my $absolute_job_dir = shift @args;

    opendir (my $DIR, $absolute_job_dir) or die "ERROR : cannot open directory $absolute_job_dir : $!";
    while ( readdir $DIR ){
        if ( $_ eq 'basic_candidates.yml'
          || $_ eq 'basic_known_candidates.yml'
          || $_ eq 'finished' ){
            unlink "$absolute_job_dir/$_";
        }
        elsif ( $_ =~ /.*_filtered.bed/ ){
            miRkwood::BEDHandler::zipBEDfile( "$absolute_job_dir/$_" );
        }

    }
    closedir $DIR;

}

1;
