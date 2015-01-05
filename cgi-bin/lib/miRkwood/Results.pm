package miRkwood::Results;

# ABSTRACT: Code directly tied to the results data structure

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use feature 'switch';
use Time::gmtime;
use File::Spec;

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


=method get_candidates_dir


=cut

sub get_candidates_dir {
	my ( $self, @args ) = @_;
	my $id_job         = shift @args;
	my $results_dir    = $self->jobId_to_jobPath($id_job);
	my $candidates_dir = File::Spec->catdir( $results_dir, 'candidates' );
	return $candidates_dir;
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
my %results = miRkwood::Results->get_structure_for_jobID($jobId);

=cut

sub get_structure_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
	my $job_dir = $self->jobId_to_jobPath($jobId);
	miRkwood->CONFIG_FILE(
		miRkwood::Paths->get_job_config_path($job_dir) );
	my $candidates_dir = $self->get_candidates_dir($jobId);
	return $self->deserialize_results($candidates_dir);
}

=method get_basic_structure_for_jobID

=cut

sub get_basic_structure_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
	my $job_dir = $self->jobId_to_jobPath($jobId);
	miRkwood->CONFIG_FILE(
		miRkwood::Paths->get_job_config_path($job_dir) );
	my $candidates_file = File::Spec->catfile( $job_dir, 'basic_candidates.yml');
	return miRkwood::get_yaml_file( $candidates_file );
}

sub get_basic_pseudoXML_for_jobID {
	my ( $self, @args ) = @_;
	my $jobId   = shift @args;
	my $results = $self->get_basic_structure_for_jobID($jobId);

	my $output = "";

    $output .= "<results id='all'>\n";
    my @candidates = sort {
        ( $a->{'name'} cmp  $b->{'name'} )
          || (
            $a->{'start_position'} <=> $b->{'start_position'} )
    } @{$results};
    
    foreach my $candidate (@candidates) {
        $output .= $self->convert_basic_to_pseudoXML($candidate);
    }
    $output .= "</results>\n";
    
    $output .= "<results id='all2'>\n";
    @candidates = sort {
        ( $b->{'quality'} cmp $a->{'quality'} )
          || (
            $a->{'start_position'} <=> $b->{'start_position'} )
    } @candidates;
    foreach my $candidate (@candidates) {
        $output .= $self->convert_basic_to_pseudoXML($candidate) . "\n";
    }
    $output .= "</results>";

    return $output;
}

sub convert_basic_to_pseudoXML {
	my ( $self, @args ) = @_;
	my $candidate = shift @args;
	my %candidate = %{$candidate};
	
	my @fields_to_truncate = ( 'mfe', 'mfei', 'amfe' );
    my $result = "<Sequence";
	my @optional_fields = miRkwood::Candidate->get_optional_candidate_fields();
    my @headers = ( 'name', 'position', 'length', 'strand', 'quality', @optional_fields, 'image', 'identifier' );

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
        $result .= " $header='$contents'";
    }
    $result .= "></Sequence>";
}

=method known_mirnas_for_jobID

Print an html table for known miRNAs, with IDs, locus, strand, and
number of reads matching to the mature and the precursor.

Input : a BED file
Output : a string containing an HTML table

=cut
sub known_mirnas_for_jobID {
	my ( $self, @args ) = @_;
    my $jobId = shift @args;
    my $bed_file = shift @args;
    my $gff_file = shift @args;
    my $genome_file = shift @args;
    my $reads_dir = miRkwood::Paths::get_known_reads_dir( $jobId );
    my $candidates_dir = miRkwood::Paths::get_known_candidates_dir( $jobId );

    my $html = '';
    my @field;
    my ($id, $name, $chromosome, $strand, $score);
    my ($precursor_name, $precursor_start, $precursor_end, $precursor_reads, $precursor_id);
    my ($precursor_of_mature, $mature_of_precursor);
    my $mature_reads;
    my $data;
    my $star = "<img src='/mirkwood/style/star.png' alt='star' style='width:15px; height:15px;' />";

    ##### Read the GFF and links each precursor with its mature

    open (my $GFF, $gff_file) or die "ERROR while opening $gff_file : $!";

    while ( <$GFF> ){
        if ( ! /^#/ ){
            chomp;

            @field = split( /\t/xms );

            if ( $field[2] eq "miRNA" and $field[8] =~ /ID=([^;]+).*Derives_from=([^;]+)/ ){
                $precursor_of_mature->{ $1 } = $2;
            }
        }
    }

    close $GFF;


    ##### Read the BED file
    open (my $BED, '<', $bed_file) or die "ERROR while opening $bed_file : $!";

    while ( <$BED> ){

        chomp;

        @field = split( /\t/xms );

        if ($field[14] =~ /ID=([^;]+).*Name=([^;]+)/ ){
            $id = $1;
            $name = $2;
        }

        if ( $field[8] eq "miRNA_primary_transcript" ){
            $precursor_id = $id;
            $data->{$precursor_id}{'identifier'}      = $id;
            $data->{$precursor_id}{'precursor_name'}  = $name;
            $data->{$precursor_id}{'name'}  = "$field[0]__$field[9]-$field[10]";
            $data->{$precursor_id}{'length'} = $field[10] - $field[9] + 1;
            $data->{$precursor_id}{'position'} = "$field[9]-$field[10]";
            $data->{$precursor_id}{'start_position'} = 1;
            $data->{$precursor_id}{'end_position'}   = $field[10] - $field[9];            
            $data->{$precursor_id}{'precursor_reads'}{"$field[1]-$field[2]"} = $field[4];
        }
        elsif ( $field[8] eq "miRNA" ){
            $precursor_id = $precursor_of_mature->{$id};
            $data->{$precursor_id}{'matures'}{$id}{'mature_name'}  = $name;
            $data->{$precursor_id}{'matures'}{$id}{'mature_start'} = $field[9];
            $data->{$precursor_id}{'matures'}{$id}{'mature_end'}   = $field[10];
            $data->{$precursor_id}{'matures'}{$id}{'mature_reads'}{"$field[1]-$field[2]"} = $field[4];
        }

        $data->{$precursor_id}{'chromosome'} = $field[0];
        $data->{$precursor_id}{'strand'}     = $field[5];

    }

    close $BED;

    ##### Treat data by precursor
    $html = '<table><tbody id="cases"><tr>';
    $html .= '<th>Name</th>';
    $html .= '<th>Chromosome</th>';
    $html .= '<th>Strand</th>';
    $html .= '<th>Position Precursor</th>';
    $html .= '<th>Number of reads</th>';
    $html .= '<th>Score</th>';
    
    $html .= "</tr>\n";
    
    my @ids = sort { $data->{$a}{'chromosome'} <=> $data->{$b}{'chromosome'}
                        ||
                     $data->{$a}{'precursor_start'} <=> $data->{$b}{'precursor_start'}
                        ||
                     $data->{$a}{'precursor_end'} <=> $data->{$a}{'precursor_end'}   } keys%{$data};
              
    foreach $precursor_id ( @ids ){

        $precursor_name  = $data->{$precursor_id}{'precursor_name'};
        $chromosome      = $data->{$precursor_id}{'chromosome'};
        $strand          = $data->{$precursor_id}{'strand'};
        $precursor_start = $data->{$precursor_id}{'start_position'};
        $precursor_end   = $data->{$precursor_id}{'end_position'};
        $precursor_reads = 0;
        $mature_reads = 0;
        $score = '';

        ##### Count number of reads
        foreach (keys %{$data->{$precursor_id}{'precursor_reads'}}){
            $precursor_reads += $data->{$precursor_id}{'precursor_reads'}{$_};
        }
        foreach my $mature_id ( keys %{$data->{$precursor_id}{'matures'}} ){
            foreach my $read ( keys %{$data->{$precursor_id}{'matures'}{$mature_id}{'mature_reads'}} ){
                $mature_reads += $data->{$precursor_id}{'matures'}{$mature_id}{'mature_reads'}{$read};
            }
        }

        ##### Calculate score
        $data->{$precursor_id}{'score'} = 0;
        if ( $precursor_reads >= 10 ){
           $data->{$precursor_id}{'score'}++; 
        }
        if ( $mature_reads >= ( $precursor_reads / 2 ) ){
            $data->{$precursor_id}{'score'}++;
        }
        $score = '<center>';
        for (my $i = 0; $i < $data->{$precursor_id}{'score'}; $i++){
            $score .= $star;
        }
        $score .= '</center>';    

        ### Create a Candidate object
        my $candidate = miRkwood::Candidate->new( $data->{$precursor_id} );
        miRkwood::CandidateHandler->serialize_candidate_information($candidates_dir, $candidate);

        ### Create individual card with reads cloud
        miRkwood::CandidateHandler::print_reads_clouds( $data->{$precursor_id}, $genome_file, $reads_dir );

        ##### Print the HTML table
        $html .= '<tr>';
        $html .= '<td><a href=' . miRkwood::Utils::make_mirbase_link($precursor_id) . ">$precursor_name</a></td>";
        $html .= "<td>$chromosome</td>";
        $html .= "<td>$strand</td>";
        $html .= "<td>$data->{$precursor_id}{'position'}</td>";
        $html .= '<td><a href=' . "./getCandidate.pl?jobId=$jobId&id=$precursor_id&type=reads" . ">$precursor_reads</a></td>";
        $html .= "<td>$score</td>";     
        $html .= "</tr>\n";

    }  
          
    $html .= '</tbody></table>';
    
    debug( "Known miRNAS stored in YAML format", miRkwood->DEBUG() );
    
    return $html;  
    
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

Retrieve the results in the given directory

Usage:
my %results = miRkwood::Results->deserialize_results( $candidates_dir );

=cut

sub deserialize_results {
	my ( $self, @args ) = @_;
	my $candidates_dir = shift @args;
	my %myResults      = ();
	opendir DIR, $candidates_dir;    #ouverture répertoire job
	my @files;
	@files = readdir DIR;
	closedir DIR;
	foreach my $file (@files)        # parcours du contenu
	{
		my $full_file = File::Spec->catfile( $candidates_dir, $file );
		if (   $file ne "."
			&& $file ne ".." )
		{
			my $candidate;
			if (
				eval {
					$candidate =
					  miRkwood::Candidate->new_from_serialized($full_file);
				}
			  )
			{
				my $identifier = $candidate->get_identifier();
				$myResults{$identifier} = $candidate;
			}
		}
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
	my $results = $self->get_basic_structure_for_jobID($jobId);
	my $size    = scalar @{$results};
	return $size;
}

1;
