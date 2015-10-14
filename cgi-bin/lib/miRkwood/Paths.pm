package miRkwood::Paths;

# ABSTRACT: Managing paths and path construction

use strict;
use warnings;
use Cwd;
use File::Spec;
use File::Basename;

use miRkwood;
use miRkwood::Results;

=method get_config

Get the configuration file contents.

=cut

sub get_config {
    my ($self, @args) = @_;
    return miRkwood->PIPELINE_CONFIG();
}

=method get_job_config_path

Given a job directory, return the path to the job configuration file

=cut

sub get_job_config_path {
    my ($self, @args) = @_;
    my ($dir_job) = shift @args;
    my $job_config_path = File::Spec->catfile( $dir_job, 'run_options.cfg' );
    return $job_config_path;
}

=method get_data_path

Return the project data directory

=cut

sub get_data_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(miRkwood->MIRKWOOD_PATH(), 'data');
}

=method get_local_programs_path

Return the project local programs directory

=cut

sub get_local_programs_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(miRkwood->MIRKWOOD_PATH(), 'programs');
}

=method get_scripts_path

Return the project static directory

=cut

sub get_scripts_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(miRkwood->MIRKWOOD_PATH(), 'scripts');
}

=method get_lib_path

Return the project static directory

=cut

sub get_lib_path {
    my ($self, @args) = @_;
    return File::Spec->catdir(miRkwood->MIRKWOOD_PATH(), 'lib');
}

=method get_results_filesystem_path

Return the path to the results, as seen from the get_results_filesystem_path

=cut

sub get_results_filesystem_path {
    my ($self, @args) = @_;
    my %config = $self->get_config();
    return $config{'filesystem_results'};
}

=method get_workspace_path

Return the path to the job workspace

=cut

sub get_workspace_path {
    my ($self, @args) = @_;
    my ($job_dir) = @args;
    return File::Spec->catdir($job_dir, 'workspace');
}

=method get_candidate_paths

Return the path for a given candidate
in the workspace directory

=cut

sub get_candidate_paths {
    my ($self, @args) = @_;
    my ($job_dir,  $dir, $subDir) = @args;
    my $workspace = $self->get_workspace_path($job_dir);
    my $candidate_dir = File::Spec->catdir($workspace, $dir, $subDir);
    return $candidate_dir;
}

=method get_candidate_dir_name

Return the basename for candidates directory
(containing YAML files)

=cut

sub get_candidate_dir_name {
    my (@args) = @_;
    return 'YML';
}

=method get_dir_candidates_path

Return the path to the candidates directory
Parameter : job id

=cut

sub get_dir_candidates_path {
    my (@args) = @_;
    my $job_id = shift @args;
    my $job_dir = miRkwood::Results->jobId_to_jobPath( $job_id );
    return File::Spec->catdir( $job_dir, get_candidate_dir_name() );
}

=method get_dir_candidates_path_from_job_dir

Return the path to the candidates directory
Parameter : job directory

=cut
sub get_dir_candidates_path_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir( $job_dir, get_candidate_dir_name() );
}

=method get_basename_for_novel_miRNA

Return the basename for new candidates directory
(used for candidates directory and reads directory)

=cut

sub get_basename_for_novel_miRNA {
    my (@args) = @_;
    return 'novel_miRNA';
}

=method get_new_candidates_dir

Return the path to the new candidates directory
Parameter : job id

=cut
sub get_new_candidates_dir {
    my (@args) = @_;
    my $job_id = shift @args;
    return File::Spec->catdir( get_dir_candidates_path($job_id), get_basename_for_novel_miRNA() );
}

=method get_new_candidates_dir_from_job_dir

Return the path to the new candidates directory
Parameter : job directory

=cut
sub get_new_candidates_dir_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir( get_dir_candidates_path_from_job_dir($job_dir), get_basename_for_novel_miRNA() );
}

=method get_basename_for_known_miRNA

Return the basename for known candidates directory
(used for candidates directory and reads directory)

=cut

sub get_basename_for_known_miRNA {
    my (@args) = @_;
    return 'known_miRNA';
}

=method get_known_candidates_dir

Return the path to the new candidates directory
Parameter : job id

=cut
sub get_known_candidates_dir {
    my (@args) = @_;
    my $job_id = shift @args;
    return File::Spec->catdir( get_dir_candidates_path($job_id), get_basename_for_known_miRNA() );
}

=method get_known_candidates_dir_from_job_dir

Return the path to the new candidates directory
Parameter : job directory

=cut
sub get_known_candidates_dir_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir( get_dir_candidates_path_from_job_dir($job_dir), get_basename_for_known_miRNA() );
}

=method get_reads_dir_name

Return the basename for reads clouds directory

=cut

sub get_reads_dir_name {
    my (@args) = @_;
    return 'read_clouds';
}

=method get_dir_reads_path

Return the path to the reads directory
Parameter : job id

=cut

sub get_dir_reads_path {
    my (@args) = @_;
    my $job_id = shift @args;
    my $job_dir = miRkwood::Results->jobId_to_jobPath( $job_id );
    return File::Spec->catdir($job_dir, get_reads_dir_name() );
}

=method get_dir_reads_path_from_job_dir

Return the path to the reads directory
Parameter : job directory

=cut

sub get_dir_reads_path_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir($job_dir, get_reads_dir_name() );
}

=method get_new_reads_dir

Return the path to the novel miRNAs reads directory
Parameter : job id

=cut
sub get_new_reads_dir {
    my (@args) = @_;
    my $job_id = shift @args;
    return File::Spec->catdir( get_dir_reads_path($job_id), get_basename_for_novel_miRNA() );
}

=method get_new_reads_dir_from_job_dir

Return the path to the novel miRNAs reads directory
Parameter : job directory

=cut
sub get_new_reads_dir_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir( get_dir_reads_path_from_job_dir($job_dir), get_basename_for_novel_miRNA() );
}

=method get_known_reads_dir

Return the path to the known miRNAs reads directory
Parameter : job id

=cut
sub get_known_reads_dir {
    my (@args) = @_;
    my $job_id = shift @args;
    return File::Spec->catdir( get_dir_reads_path($job_id), get_basename_for_known_miRNA() );
}

=method get_known_reads_dir_from_job_dir

Return the path to the known miRNAs reads directory
Parameter : job id

=cut
sub get_known_reads_dir_from_job_dir {
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir( get_dir_reads_path_from_job_dir($job_dir), get_basename_for_known_miRNA() );
}

=method get_dir_reads_path_from_job_dir_and_mirna_type

Return the path to the reads directory for novel or known
candidates.
Parameters : - job dir
             - mirna type ('novel_miRNA' or 'known_miRNA')

=cut
sub get_dir_reads_path_from_job_dir_and_mirna_type{
    my (@args) = @_;
    my $job_dir = shift @args;
    my $mirna_type = shift @args;
    if ( $mirna_type eq 'known_miRNA' ){
        return get_known_reads_dir_from_job_dir( $job_dir );
    }
    else{
        return get_new_reads_dir_from_job_dir( $job_dir );
    }
}


=method get_images_dir_name

  Return the basename for images directory
  (images created by VARNA)

=cut
sub get_images_dir_name{
    my (@args) = @_;
    return 'images';
}

=method get_dir_images_path_from_job_dir

Return the path to the images directory.
Parameters :  job dir

=cut
sub get_dir_images_path_from_job_dir{
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir($job_dir, get_images_dir_name() );
}


=method get_alignments_dir_name

  Return the basename for alignments directory

=cut
sub get_alignments_dir_name{
    my (@args) = @_;
    return 'alignments';
}

=method get_dir_alignments_path_from_job_dir

Return the path to the alignments directory.
Parameters :  job dir

=cut
sub get_dir_alignments_path_from_job_dir{
    my (@args) = @_;
    my $job_dir = shift @args;
    return File::Spec->catdir($job_dir, get_alignments_dir_name() );
}


=method get_results_folder_basename_for_CLI

  Return the basename for results directory
  (CLI pipeline)

=cut
sub get_results_folder_basename_for_CLI {
    my (@args) = @_;
    return 'results';
}

=method get_results_folder_for_CLI_from_job_dir

  Return the directory in which the results will
  be stored (CLI pipeline), either :
  - jobDir/results
  - jobDir/results/known_miRNA
  - jobDir/results/novel_miRNA

=cut
sub get_results_folder_for_CLI_from_job_dir {
    my (@args) = @_;
    my $job_dir       = shift@args;
    my $pipeline_type = shift @args;
    my $mirna_type    = shift @args;

    my $results_folder = create_folder( File::Spec->catdir( $job_dir, get_results_folder_basename_for_CLI() ) );

    if ( $pipeline_type eq 'smallRNAseq' ){
        if ( $mirna_type eq 'known_miRNA' ){
            return create_folder( File::Spec->catdir( $results_folder, get_basename_for_known_miRNA() ) );
        }
        else{
            return create_folder( File::Spec->catdir( $results_folder, get_basename_for_novel_miRNA() ) );
        }
    }
    else{
        return $results_folder;
    }
}

=method get_pieces_folder_basename_for_CLI

  Return the basename for results directory
  (CLI pipeline)

=cut
sub get_pieces_folder_basename_for_CLI {
    my (@args) = @_;
    return 'pieces';
}

=method get_bed_file

Return the path to the BED file corresponding to given type
'', '_filtered', '_CDS, '_miRNAs', '_otherRNA', '_multimapped'
for the given job ID.
Parameter : job id and bed type

=cut
sub get_bed_file {
    my (@args) = @_;
    my $absolute_job_dir = shift @args;
    my $bed_type = shift @args; # should be '', '_filtered', '_CDS, '_miRNAs', '_otherRNA', '_multimapped'
    my $extension = shift @args;
    my $run_options_file = miRkwood::Paths->get_job_config_path($absolute_job_dir);
    miRkwood->CONFIG_FILE($run_options_file);
    my $cfg      = miRkwood->CONFIG();
    my $bed_name = $cfg->param('job.bed');

    return File::Spec->catdir( $absolute_job_dir, $bed_name. $bed_type . ".$extension");
}

=method get_workspace_chromosome_dir

  Return the path to a given chromosome directory
  in the workspace directory
  Parameters : - workspace directory
               - chromosome id

=cut
sub get_workspace_chromosome_dir {
    my (@args) = @_;
    my $workspace_dir = shift @args;
    my $chromosome = shift @args;
    return File::Spec->catdir( $workspace_dir, $chromosome);
}

=method get_workspace_candidate_dir

  Return the path to a Candidate directory
  in the workspace directory
  Parameters : - workspace directory
               - chromosome id
               - cluster positions
               - strand  

=cut
sub get_workspace_candidate_dir {
    my (@args) = @_;
    my $workspace_dir = shift @args;
    my $chromosome = shift @args;
    my $cluster = shift @args;
    my $strand = shift @args;

    my $chromosome_dir = get_workspace_chromosome_dir( $workspace_dir, $chromosome );
    return File::Spec->catdir( $chromosome_dir, $cluster . $strand );
}

=method get_bed_size_file_name

  Return the name of the text file
  which contains the number of reads
  in each BED file.

=cut
sub get_bed_size_file_name {
    my (@args) = @_;
    return 'bed_sizes.txt';
}

=method get_orphan_hairpin_file_name

  Return the name of the bed file which
  contains data about the orphan hairpins,
  that is the candidates with a quality score of 0
  and no conservation.

=cut
sub get_orphan_hairpin_file_name {
    my (@args) = @_;
    my $basename_bed = shift @args;
    return $basename_bed.'_orphan_hairpins.bed';
}

=method create_folder

  Method to create a folder if it not already exists
  Parameter : the folder name

=cut
sub create_folder {
    my (@args) = @_;
    my $folder = shift @args;
    if ( !-e $folder ) {
        mkdir $folder or die "ERROR when creating $folder : $!";
    }
    return $folder;
}

1;
