package PipelineMiRNA::Results;

# ABSTRACT: Code directly tied to the results data structure

use strict;
use warnings;

use feature 'switch';
use Time::gmtime;

use PipelineMiRNA::Candidate;
use PipelineMiRNA::Utils;

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
    my $results_dir = PipelineMiRNA::Paths->get_results_filesystem_path();
    my $jobPath = File::Spec->catdir( $results_dir, $dirJob_name);
    return $jobPath;
}

=method get_candidates_dir


=cut

sub get_candidates_dir {
    my ( $self, @args ) = @_;
    my $id_job      = shift @args;
    my $results_dir = $self->jobId_to_jobPath($id_job);
    my $candidates_dir = File::Spec->catdir( $results_dir, 'candidates');
    return $candidates_dir;
}

=method is_valid_jobID

Test whether a jobID is valid - ie if there are results for it.

=cut

sub is_valid_jobID {
    my ( $self, @args ) = @_;
    my $id_job          = shift @args;
    my $full_path = $self->jobId_to_jobPath($id_job);
    return (-e $full_path);
}

=method get_structure_for_jobID

Get the results structure of a given job identifier

Usage:
my %results = PipelineMiRNA::Results->get_structure_for_jobID($jobId);

=cut

sub get_structure_for_jobID {
    my ( $self, @args ) = @_;
    my $jobId   = shift @args;
    my $job_dir = $self->jobId_to_jobPath($jobId);
    PipelineMiRNA->CONFIG_FILE(PipelineMiRNA::Paths->get_job_config_path($job_dir));
    my $candidates_dir = $self->get_candidates_dir($jobId);
    return $self->deserialize_results($candidates_dir);
}

=method has_candidates

Parse and serialize the results structure of $job_dir

Usage:
PipelineMiRNA::Results->has_candidates( \%myResults );

=cut

sub has_candidates {
    my ( $self, @args ) = @_;
    my $results = shift @args;
    my %results = %{$results};
    return ( keys %results > 0);
}

=method deserialize_results

Retrieve the results in the given directory

Usage:
my %results = PipelineMiRNA::Results->deserialize_results( $candidates_dir );

=cut

sub deserialize_results {
    my ( $self, @args ) = @_;
    my $candidates_dir   = shift @args;
    my %myResults = ();
    opendir DIR, $candidates_dir;    #ouverture répertoire job
    my @files;
    @files = readdir DIR;
    closedir DIR;
    foreach my $file (@files)    # parcours du contenu
    {
        my $full_file = File::Spec->catfile( $candidates_dir, $file );
        if (    $file ne "."
             && $file ne "..")
        {
            my %candidate;
            %candidate = PipelineMiRNA::Candidate->deserialize_candidate($full_file);
            if (! eval { %candidate = PipelineMiRNA::Candidate->deserialize_candidate($full_file) } ){
                # Catching, do nothing
            }else{
                my $identifier = $candidate{'identifier'};
                $myResults{$identifier} = \%candidate;
            }
        }
    }
    return %myResults;
}


=method serialize_results

Parse and serialize the results structure of $job_dir

Usage:
PipelineMiRNA::Results->serialize_results($job_dir);

=cut

sub serialize_results {
    my ( $self, @args ) = @_;
    my $relative_job_dir   = shift @args;
    my $job_dir = PipelineMiRNA::Paths->get_absolute_path($relative_job_dir);
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
                    if (! eval { PipelineMiRNA::Candidate->serialize_candidate_information($relative_job_dir, $dir, $subDir)} ) {
                        # Catching
                    }else{
                        # All is well
                    }
                }
            }
        }
    }
    return;
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
            given ($export_type) {
                when (/fas/) { $output .= PipelineMiRNA::Candidate->candidateAsFasta($value); }
                when (/dot/) { $output .= PipelineMiRNA::Candidate->candidateAsVienna($value); }
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
    my @headers1 = ('name', 'position','length','strand','quality', @optional_fields);
    my @headers2 = ('Vienna', 'DNASequence', 'identifier');

    my $result = "<results id='all'>\n";
    my @keys = sort { PipelineMiRNA::Utils::get_element_of_split($results{$a}{'position'}, '-', 0) <=>
                      PipelineMiRNA::Utils::get_element_of_split($results{$b}{'position'}, '-', 0)
                    } keys %results;

    foreach my $key (@keys) {
        my $value = $results{$key};
        $result .= "<Sequence";
        for my $header (@headers1){
            my $contents = ${$value}{$header};
            if (!$contents){
                $contents = q{};
            }
            $result .= " $header='$contents'";
        }
        my $img = PipelineMiRNA::Candidate->get_relative_image($value);
        $result .= " image='$img'";
        for my $header (@headers2){
            $result .= " $header='${$value}{$header}'";
        }
        $result .= "></Sequence>\n";
    }
    $result .= "</results>\n";
    $result .= "<results id='all2'>\n";
    @keys = sort keys %results;
    @keys = sort { ( $results{$b}{'quality'} cmp
                     $results{$a}{'quality'} )
                   ||
                   ( PipelineMiRNA::Utils::get_element_of_split($results{$a}{'position'}, '-', 0) <=>
                     PipelineMiRNA::Utils::get_element_of_split($results{$b}{'position'}, '-', 0) )
                 } keys %results;
    foreach my $key (@keys) {
        my $value = $results{$key};
        $result .= "<Sequence";
        for my $header (@headers1){
            my $contents = ${$value}{$header};
            if (!$contents){
                $contents = q{};
            }
            $result .= " $header='$contents'";
        }
        my $img = PipelineMiRNA::Candidate->get_relative_image($value);
        $result .= " image='$img'";
        for my $header (@headers2){
            $result .= " $header='${$value}{$header}'";
        }
        $result .= "></Sequence>\n";
    }
    $result .= "</results>";
    return $result;
}

1;
