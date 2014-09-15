package miRkwood::SequenceJob;

# ABSTRACT: Job processing a given sequence

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use miRkwood::CandidateJob;
use miRkwood::SequenceSubJob;
use miRkwood::Utils;

=method new

Constructor

my $sequence_job = miRkwood::SequenceJob->new($sequence_dir, $name, $sequence);

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($directory, $identifier, $name, $sequence) = @args;
    my $self = {
        directory => $directory,
        identifier => $identifier,
        name => $name,
        sequence => $sequence,
    };
    bless $self, $class;
    return $self;
}

sub get_directory {
    my ($self, @args)  = @_;
    return $self->{'directory'};
}

sub run {
    my ($self, @args) = @_;
    my $candidates_array = $self->get_raw_candidates();
    my %candidates_hash = $self->process_raw_candidates($candidates_array);
    return $self->process_candidates( \%candidates_hash );
}

=method get_raw_candidates

Get the candidates for the sequence

 Usage : $self->get_raw_candidates();

=cut

sub get_raw_candidates{
    my ($self, @args) = @_;
    my $candidates1 = $self->get_forward_strand_candidates();
    my $candidates2 = $self->get_other_strand_candidates();
    my $candidates_array = $self->merge_and_sort_candidates_array($candidates1, $candidates2);
    return $candidates_array;
}

sub get_forward_strand_candidates {
    my ($self, @args) = @_;
    my $sequence_subjob = miRkwood::SequenceSubJob->new($self->get_directory(), $self->{'name'}, $self->{'sequence'}, '+');
    my $candidates = $sequence_subjob->get_raw_candidates_for_sequence();
}

=method get_other_strand_candidates

Get the candidates for the sequence
Return an array of candidates.

 Usage : my $candidates = $self->get_other_strand_candidates();

=cut

sub get_other_strand_candidates {
    my ($self, @args) = @_;
    my $candidates = [];
    my $cfg = miRkwood->CONFIG();
    if ( $cfg->param('options.strands') ) {
        debug( 'Processing the other strand', miRkwood->DEBUG() );
        my $reversed_sequence = miRkwood::Utils::reverse_complement($self->{'sequence'});
        my $sequence_subjob2 = miRkwood::SequenceSubJob->new($self->get_directory(), $self->{'name'}, $reversed_sequence, '-');
        $candidates = $sequence_subjob2->get_raw_candidates_for_sequence();
    }
    return $candidates;
}

=method merge_and_sort_candidates_array

Merge both given candidates list
Works also if one of the list is empty

=cut

sub merge_and_sort_candidates_array {
    my ($self, @args) = @_;
    my ($candidates1, $candidates2) = @args;
    my @candidates_array = sort { $a->{start_position} <=> $b->{start_position} } ( @{$candidates1}, @{$candidates2} );
    return \@candidates_array;
}

sub process_raw_candidates{
    my ($self, @args) = @_;
    my $candidates_array = shift @args;

    my @candidates_array = @{$candidates_array};
    my $cfg = miRkwood->CONFIG();
    if ( $cfg->param('options.mfe') ) {
        debug('Select only sequences with MFEI < -0.6', miRkwood->DEBUG() );
        @candidates_array = grep { mfei_below_threshold($_, -0.6) } @candidates_array;
    }

    my %candidates_hash;
    if (@candidates_array) {
        debug('Merging candidates', miRkwood->DEBUG() );
        %candidates_hash = $self->merge_candidates( \@candidates_array );
    }
    else {
        %candidates_hash = ();
    }
    return %candidates_hash;
}

=method mfei_below_threshold

Return whether a given candidate has its mfei above a given threshold

=cut

sub mfei_below_threshold {
    my @args      = @_;
    my $candidate = shift @args;
    my $threshold = shift @args;
    my %current_candidate = %{ $candidate };
    my $mfei = $current_candidate{'mfei'};
    return $mfei < $threshold;
}

=method merge_candidates

Process the candidates and try merging them.
We assume the candidates array is already sorted by growing position

=cut

sub merge_candidates {
    my ($self, @args) = @_;
    my $candidates_array = shift @args;
    my @candidates_array = @{$candidates_array};
    my $nb_candidates = scalar @candidates_array;

    my @merged_candidates   = ();
    my %final_hash          = ();
    my %reference_candidate = %{ $candidates_array[0] };
    my %best_candidate      = %reference_candidate;
    my %current_candidate;
    for my $candidate_index ( 1 .. $#candidates_array ) {
        %current_candidate = %{ $candidates_array[$candidate_index] };
        my $start = $current_candidate{'start_position'};
        my $end   = $current_candidate{'end_position'};
        my ( $ref_start, $ref_end ) =
          ( $reference_candidate{'start_position'}, $reference_candidate{'end_position'} );
        if ( is_included( $start, $end, $ref_start, $ref_end ) ) {
            if ( $best_candidate{'mfei'} <= -0.8 ) {
                push @merged_candidates, {%current_candidate};
            }
            else {
                if ( $current_candidate{'mfei'} < $best_candidate{'mfei'} ) {
                    push @merged_candidates, {%best_candidate};
                    %best_candidate = %current_candidate;
                }
                else {
                    push @merged_candidates, {%current_candidate};
                }
            }

        }
        elsif ( is_overlapping( $start, $end, $ref_start, $ref_end ) ) {
            if ( $current_candidate{'mfei'} < $best_candidate{'mfei'} ) {
                push @merged_candidates, {%best_candidate};
                %best_candidate = %current_candidate;
            }
            else {
                push @merged_candidates, {%current_candidate};
            }
        }
        else {
            my $final_name = $best_candidate{'name'};
            $final_hash{$final_name}                 = {};
            $final_hash{$final_name}{'max'}          = {%best_candidate};
            $final_hash{$final_name}{'alternatives'} = [@merged_candidates];
            @merged_candidates                       = undef;
            @merged_candidates                       = ();
            %reference_candidate                     = %current_candidate;
            %best_candidate                          = %reference_candidate;
        }

    }    #foreach
    my $final_name = $best_candidate{'name'};
    $final_hash{$final_name}                 = {};
    $final_hash{$final_name}{'max'}          = {%best_candidate};
    $final_hash{$final_name}{'alternatives'} = [@merged_candidates];
    return %final_hash;
}


=method is_overlapping

Test whether one candidate is overlapping with the other
(based on their positions).

=cut

sub is_overlapping {
    my @args      = @_;
    my $start     = shift @args;
    my $end       = shift @args;
    my $ref_start = shift @args;
    my $ref_end   = shift @args;
    ( $start ne '' && $end ne '' && $ref_start ne '' && $ref_end ne '' ) or die('Not enough values provided');
    $ref_start <= $start or die("Positions should be ordered : $ref_start <= $start");
    return ( $start < ( $ref_start + $ref_end ) / 2 );
}

=method is_included

Test whether one candidate is included into the other
(based on their positions).

=cut

sub is_included {
    my @args      = @_;
    my $start     = shift @args;
    my $end       = shift @args;
    my $ref_start = shift @args;
    my $ref_end   = shift @args;
    ( $start ne '' && $end ne '' && $ref_start ne '' && $ref_end ne '' ) or die('Not enough values provided');
    $ref_start <= $start or die("Positions should be ordered : $ref_start <= $start");
    return ( $end <= $ref_end );
}

=method process_candidates

Run the rest of the pipeline on the candidates

=cut

sub process_candidates {
    my ($self, @args) = @_;
    my (%candidates_hash) = %{ shift @args };
    my @candidates_result;
    my $candidate_identifier = 0;
    foreach my $key ( sort keys %candidates_hash ) {
        $candidate_identifier++;
        my $candidate_dir = $self->create_candidate_directory($candidate_identifier);
        my $candidate = $candidates_hash{$key};
        push @candidates_result, $self->run_pipeline_on_candidate($candidate_identifier, $candidate);
    }
    return \@candidates_result;
}

sub run_pipeline_on_candidate {
    my ($self, @args) = @_;
    my $candidate_identifier = shift @args;
    my $candidate_structure = shift @args;

    my $candidate_dir = $self->create_candidate_directory($candidate_identifier);
    my $sequence_identifier = $self->{'identifier'};
    my $candidate_full_identifier = "$sequence_identifier-$candidate_identifier";
    my $candidate_ref = $candidate_structure->{'max'};
    my $alternatives = $candidate_structure->{'alternatives'};
    my $candidatejob = miRkwood::CandidateJob->new($candidate_dir, $self->{'name'},
                                                   $candidate_full_identifier,
                                                   $candidate_ref, $alternatives);
    return $candidatejob->run();
}

sub create_candidate_directory {
    my ($self, @args) = @_;
    my $candidate_identifier = shift @args;
    my $candidate_dir =
          File::Spec->catdir( $self->get_directory(), $candidate_identifier );
    mkdir $candidate_dir;
    return $candidate_dir;
}

1;