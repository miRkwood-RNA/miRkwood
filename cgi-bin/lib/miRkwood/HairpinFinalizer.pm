package miRkwood::HairpinFinalizer;

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use miRkwood;

use constant MFEI_THRESHOLD => -0.6;


=method process_hairpin_candidates

  Process 'raw' hairpin to 
  - filter on MFEI
  - merge overlapping candidates
  
  Usage : miRkwood::HairpinFinalizer::process_hairpin_candidates( @candidates_array );

=cut
sub process_hairpin_candidates{
    my (@args) = @_;
    my $candidates_array = shift @args;
    my @candidates_array = @{$candidates_array};
    my $cfg = miRkwood->CONFIG();
    debug('  - Process hairpin candidates', miRkwood->DEBUG() );
    if ( $cfg->param('options.mfei') ) {
        debug('     Select only sequences with MFEI < ' . MFEI_THRESHOLD, miRkwood->DEBUG() );
        @candidates_array = grep { mfei_below_threshold($_, MFEI_THRESHOLD) } @candidates_array;
    }

    my %candidates_hash;
    if (@candidates_array) {
        debug('     Merging candidates', miRkwood->DEBUG() );
        %candidates_hash = merge_candidates( \@candidates_array );
    }
    else {
        %candidates_hash = ();
    }
    return \%candidates_hash;
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
    my (@args) = @_;
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

=method process_mirna_candidates

  Run the rest of the pipeline on the candidates

=cut

sub process_mirna_candidates {
    my (@args) = @_;
    my $workspace_dir = shift @args;
    my (%candidates_hash) = %{ shift @args };
    my $sequence_identifier = shift @args;
    my $sequence_name = shift @args;    
    my @candidates_result;
    my $candidate_identifier = 0;
    debug('  - Process miRNA candidates', miRkwood->DEBUG() );
    foreach my $key ( sort keys %candidates_hash ) {
        $candidate_identifier++;
        debug( "     - Process candidate $candidate_identifier", 1);
        my $candidate = $candidates_hash{$key};
        push @candidates_result, run_pipeline_on_candidate($workspace_dir,
                                                           $candidate_identifier,
                                                           $candidate,
                                                           $sequence_identifier,
                                                           $sequence_name);
    }
    return \@candidates_result;
}


sub run_pipeline_on_candidate {
    my (@args) = @_;
    my $workspace_dir = shift @args;
    my $candidate_identifier = shift @args;
    my $candidate = shift @args;
    my $sequence_identifier = shift @args;
    my $sequence_name = shift @args;

    # Create candidate directory
    my $candidate_dir = create_candidate_directory($workspace_dir, $candidate, $candidate_identifier);

    # Run CandidateJob on the Candidate
    my $candidate_full_identifier = "$sequence_identifier-$candidate_identifier";
    my $candidate_ref = $candidate->{'max'};
    my $alternatives = $candidate->{'alternatives'};
    my $candidatejob = miRkwood::CandidateJob->new($candidate_dir,
                                                   $sequence_name,
                                                   $candidate_full_identifier,
                                                   $candidate_ref,
                                                   $alternatives);
    return $candidatejob->run();
}


sub create_candidate_directory {
    my (@args) = @_;
    my $workspace_dir = shift @args;
    my $candidate = shift @args;
    my $candidate_identifier = shift @args;
    my $chromosome = '';
    if ( $candidate->{'max'}{'name'} =~ /([^_]+)__/ ){
        $chromosome = $1;
    }
    my $cluster_directory = File::Spec->catdir( $workspace_dir, $chromosome, $candidate->{'max'}{'cluster'} . $candidate->{'max'}{'strand'} ); # TO CHECK : ok for WebBAM but for genomic... ?
    my $candidate_dir = File::Spec->catdir( $cluster_directory, $candidate_identifier );
    mkdir $candidate_dir;
    return $candidate_dir;
}

1;

