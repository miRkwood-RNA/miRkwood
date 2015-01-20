package miRkwood::CandidateHandler;

# ABSTRACT: Code to manipulate around Candidate objects

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use miRkwood::Candidate;
use miRkwood::Utils;
use YAML::XS;

=method retrieve_candidate_information

Check correctness and get the result for a given candidate

Arguments:
- $job - the job directory
- $id - the candidate identifier

Returns:
 A miRkwood::Candidate instance
=cut

sub retrieve_candidate_information {
    my ( $self, @args ) = @_;
    my $job = shift @args;
    my $id = shift @args;
    my $candidate_filepath = $self->get_candidate_filepath($job, $id);
    if ( !-e $candidate_filepath ){
        die("Unvalid candidate information in $candidate_filepath ");
    }else{
        return miRkwood::Candidate->new_from_serialized($candidate_filepath);
    }
}

=method get_candidate_filepath

Get the candidate filepath given its identifier and the job directory

=cut

sub get_candidate_filepath {
    my ( $self, @args ) = @_;
    my $job = shift @args;
    my $id = shift @args; 
    if ( -e File::Spec->catfile($job, 'candidates/known', $self->make_candidate_filename($id)) ) {
        return File::Spec->catfile($job, 'candidates/known', $self->make_candidate_filename($id));
    }
    elsif ( -e File::Spec->catfile($job, 'candidates/new', $self->make_candidate_filename($id)) ) {
        return File::Spec->catfile($job, 'candidates/new', $self->make_candidate_filename($id));
    } 
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

Serialize the given candidate on disk

Arguments:
- $serialization_path - the filepath to serialize to
- %candidate - the candidate

=cut

sub serialize_candidate_information {
    my ( $self, @args ) = @_;
    my $serialization_path = shift @args;
    my $candidate_object = shift @args;
    my $candidate_base_filename = $self->make_candidate_filename($candidate_object->{'identifier'});
    my $serialization_file = File::Spec->catfile($serialization_path, $candidate_base_filename);
    my %converted_hash = %{$candidate_object};
    return YAML::XS::DumpFile($serialization_file, %converted_hash);
}

=method print_reads_clouds

Create a txt file with reads cloud
Display miRBase tags :
- for new candidates with alignments to miRBase
- for known candidates (corresponding to miRBase miRNAs)

=cut

sub print_reads_clouds {

    my @args = @_;
    my $mirna = shift @args;
    my $genome_db = shift @args;
    my $output_dir = shift @args;

    my $i;
    my $mature_id;
    my @positions_tags = ();

    my $cfg = miRkwood->CONFIG();
    my $genome_name = $cfg->param('job.plant');

    my $reads = {};

    my $precursor_id    = $mirna->{'identifier'};
    my $name            = ( $mirna->{'precursor_name'} or $mirna->{'name'} );
    my $strand          = $mirna->{'strand'};
    my $precursor_start = $mirna->{'start_position'};
    my $precursor_end   = $mirna->{'end_position'};
    my $chromosome      = $mirna->{'name'};

    $reads              = ( $mirna->{'precursor_reads'} or $mirna->{'reads'} );

    if ( !defined($reads) or $reads eq {} ){
        debug( "Cannot print the reads cloud for candidate $mirna->{'identifier'}", miRkwood->DEBUG() );  
        return;
    }

    if ( $mirna->{'name'} =~ /([^_]+)__(\d+)-(\d+)/ ){
        $chromosome = $1;
    }

    $reads = miRkwood::Utils::truncate_reads_out_of_candidate( $reads, $precursor_start, $precursor_end );

    my $reference = substr( $genome_db->{$chromosome}, $precursor_start -1, $precursor_end - $precursor_start);

    my $cloud_file = "$output_dir/$precursor_id.txt";
    open (my $OUT, '>', $cloud_file) or die "ERROR while creating $cloud_file : $!";

    ### Print the header
    print $OUT "$name\n";
    print $OUT "Genome : $genome_name\n";
    print $OUT "Locus  : $chromosome:$precursor_start-$precursor_end\n";
    print $OUT "Strand : $strand\n";
    print $OUT "\n$reference\n";

    ### Print miRBase tags

    if ( defined($mirna->{'alignment'}) and $mirna->{'alignment'} > 0  ){   # case cli bam pipeline with alignment to mirbase

        @positions_tags = sort {
            ( miRkwood::Utils::get_element_of_split( $a, '-', 0 )
                  <=> miRkwood::Utils::get_element_of_split( $b, '-', 0 ) )
              || ( miRkwood::Utils::get_element_of_split( $a, '-', 1 )
                <=> miRkwood::Utils::get_element_of_split( $b, '-', 1 ) )
        } keys %{$mirna->{'alignments'}};

    }

    elsif ( exists( $mirna->{'matures'} ) ){    # case bedpipeline for known mirnas

        foreach $mature_id (keys %{$mirna->{'matures'}}){
            my $mature_start = $mirna->{'matures'}{$mature_id}{'mature_start'};
            my $mature_end = $mirna->{'matures'}{$mature_id}{'mature_end'};
            my $relative_tag_start = $mature_start - $precursor_start +1;
            my $relative_tag_end = $relative_tag_start + $mature_end - $mature_start;
            push @positions_tags, "$relative_tag_start-$relative_tag_end";
        }

        @positions_tags = sort { miRkwood::Utils::get_element_of_split($a, '-', 0) <=> miRkwood::Utils::get_element_of_split($b, '-', 0)
                                ||
                                 miRkwood::Utils::get_element_of_split($a, '-', 1) <=> miRkwood::Utils::get_element_of_split($b, '-', 1)
                                } @positions_tags;

    }

    if ( @positions_tags ){

        my $placed = [];
        my $is_placed;

        foreach my $position (@positions_tags){
            if ( scalar(@{$placed}) == 0 ){
                push @{$placed->[0]}, $position;
                next;
            }

            $is_placed = 0;
            $i = 0;
            while ( $is_placed == 0 && $i < scalar(@$placed) ){
                if ( ! miRkwood::Utils::mirbase_tags_overlapping( $placed->[$i][-1], $position ) ){
                    push @{$placed->[$i]}, $position;
                    $is_placed = 1;
                }
                $i++;
            }
            if ( $i == scalar(@$placed) && $is_placed == 0){
                push @{$placed->[$i]}, $position;
            }
        }

        my $tag = "";

        foreach my $line ( @$placed ){

            $tag = "";

            for ($i = 0; $i < miRkwood::Utils::get_element_of_split( $line->[0], '-', 0) -1; $i++){
                $tag .= " ";
            }
            $tag .= miRkwood::Utils::create_mirbase_tag( 
                    miRkwood::Utils::get_element_of_split( $line->[0], '-', 0), 
                    miRkwood::Utils::get_element_of_split( $line->[0], '-', 1) );

            for ($i = 1; $i < scalar(@$line); $i++){
                for (my $j = 0; $j < miRkwood::Utils::get_element_of_split( $line->[$i], '-', 0) - miRkwood::Utils::get_element_of_split( $line->[$i-1], '-', 1) -1; $j++){
                    $tag .= " ";
                }
                $tag .= miRkwood::Utils::create_mirbase_tag(
                        miRkwood::Utils::get_element_of_split( $line->[$i], '-', 0),
                        miRkwood::Utils::get_element_of_split( $line->[$i], '-', 1) );    
            }
            print $OUT "$tag\n";

        }

    }

    ### Print the reads
    my @sorted_positions = sort { miRkwood::Utils::get_element_of_split( $a, "-", 0 ) <=> miRkwood::Utils::get_element_of_split( $b, "-", 0 )
                                    ||
                                  miRkwood::Utils::get_element_of_split( $a, "-", 1 ) <=> miRkwood::Utils::get_element_of_split( $b, "-", 1 )  
                              } (keys %$reads);

    foreach my $position (@sorted_positions){

        my $read_start = miRkwood::Utils::get_element_of_split( $position, "-", 0 );
        my $read_end = miRkwood::Utils::get_element_of_split( $position, "-", 1 );
        my $read_length = $read_end - $read_start;

        for ($i = 0; $i < $read_start - $precursor_start; $i++){
            print $OUT ".";
        }
        for ($i = 0; $i < $read_length; $i++){
            print $OUT "*";
        }
        for ($i = 0; $i < $precursor_end - $read_end; $i++){
            print $OUT ".";
        }
        print $OUT " length=$read_length depth=$reads->{$position}\n";

    }

    close $OUT;
    return;

}

1;
