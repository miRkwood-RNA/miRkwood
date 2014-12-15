package miRkwood::CandidateHandler;

# ABSTRACT: Code to manipulate around Candidate objects

use strict;
use warnings;

use miRkwood::Candidate;
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

sub print_reads_cloud {
    my ( $self, @args ) = @_;
    my $path = shift @args;
    my $genome_file = shift @args;
    my $candidate = shift @args;
    
    my $cloud_file = File::Spec->catfile($path, $candidate->{'identifier'} . ".txt");
    open(OUT, ">$cloud_file") or die "ERROR while creating file $cloud_file.\n";
    
    my $chromosome = "";
    my $start_cluster = -1;
    my $end_cluster = -1;
    my $locus = $candidate->{'name'};
    
    if ( $locus =~ /([^:]+)__(\d+)-(\d+)/ ){
        $chromosome = $1;
        $start_cluster = $2;
        $end_cluster = $3;
    }

    my $reference = miRkwood::Utils::get_sequence_from_positions( $genome_file, 
        $chromosome, $start_cluster, $end_cluster );

    ### Print the header
    
    my $i;
    my $len_reference = length($reference);
    my $start_structure = -1;
    if ( $candidate->{'position'} =~ /(\d+)-\d+/ ){
        $start_structure = $1;
    }
    
    print OUT "Genome:$genome_file\n";
    print OUT "Locus: $chromosome:$start_cluster-$end_cluster\n\n";
    print OUT "$reference\n";
    
    for ($i = 0; $i < $start_structure -1; $i++){
        print OUT " ";
    }
    print OUT "$candidate->{'sequence'}";
    for ($i = 0; $i < $len_reference - length($candidate->{'sequence'}) - $start_structure +1; $i++ ){
        print OUT " ";
    }
    print OUT " optimal stemloop predicted by RNAstemloop\n";    
    
    for ($i = 0; $i < $start_structure -1; $i++){
        print OUT " ";
    }
    print OUT "$candidate->{'structure_stemloop'}\n";
    
    # If there are alignments with miRBase, print mirbase tags
    if ( $candidate->{'alignment'} > 0 ){
    
        my @positions = sort {
            ( miRkwood::Utils::get_element_of_split( $a, '-', 0 )
                  <=> miRkwood::Utils::get_element_of_split( $b, '-', 0 ) )
              || ( miRkwood::Utils::get_element_of_split( $a, '-', 1 )
                <=> miRkwood::Utils::get_element_of_split( $b, '-', 1 ) )
        } keys %{$candidate->{'alignments'}};
        
        
        my $placed = [];
        my $is_placed;
        
        foreach my $position ( @positions ){
            
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
            for ($i = 0; $i < miRkwood::Utils::get_element_of_split( $line->[0], '-', 0) + $candidate->{'start_position'} -2; $i++){
                $tag .= " ";
            }
            
            $tag .= miRkwood::Utils::create_mirbase_tag( miRkwood::Utils::get_element_of_split( $line->[0], '-', 0 ), 
                                                         miRkwood::Utils::get_element_of_split( $line->[0], '-', 1 ) );            
            
            for ($i = 1; $i < scalar(@$line); $i++){
                for (my $j = 0; $j < miRkwood::Utils::get_element_of_split( $line->[$i], '-', 0) - miRkwood::Utils::get_element_of_split( $line->[$i-1], '-', 1) -1; $j++){
                    $tag .= " ";
                }
                
                $tag .= miRkwood::Utils::create_mirbase_tag( miRkwood::Utils::get_element_of_split( $line->[$i], '-', 0 ), 
                                                             miRkwood::Utils::get_element_of_split( $line->[$i], '-', 1 ) );                           
            }
            print OUT "$tag\n";
            
        }       
        
    }
    
    
    ### Print the reads
    my $symbol = ".";
    
    my $reads = $candidate->{'reads'};
    my @sorted_positions = sort { $a <=> $b } (keys %$reads);
    
    foreach my $position (@sorted_positions){
        
        my @sorted_sequences = sort (keys %{$reads->{$position}});
        
        foreach my $sequence (@sorted_sequences) {
            
            my $relative_position = $position - $start_cluster;
            
            for ($i = 0; $i < $relative_position; $i++){
                print OUT $symbol;
            }
            print OUT $sequence;
            for ($i = 0; $i < $len_reference - length($sequence) - $relative_position; $i++){
                print OUT $symbol;
            }
            print OUT " length=" . length($sequence) . " depth=$reads->{$position}{$sequence}{'count'}"
               . " strand=(" . $reads->{$position}{$sequence}{'strand'} . ")\n";
        }
        
    } 
    
    close OUT;   
    
}



1;
