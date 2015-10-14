package miRkwood::PrecursorBuilder;

# ABSTRACT: methods to get a miRNA precursor Candidate from an Hairpin Candidate

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use miRkwood;
use miRkwood::Utils;
use miRkwood::Paths;


=method new
  
  Constructor
  my $precursorBuilder = miRkwood::PrecursorBuilder->new( $workspace_dir, $chromosome_id, $chromosome_name );
  
=cut
sub new {
    my ( $class, @args ) = @_;
    my ( $workspace_dir, $genome_db, $chromosome_id, $chromosome_name ) = @args;
    my $self = {
        'workspace_dir'   => $workspace_dir,
        'genome_db'       => $genome_db,
        'chromosome_id'   => $chromosome_id,
        'chromosome_name' => $chromosome_name,
    };
    bless $self, $class;
    return $self;
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
    return \%final_hash;
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
    my ( $self, @args ) = @_;
    my (%candidates_hash) = %{ shift @args };
    my @candidates_result;
    my $candidate_identifier = 0;
    my $cfg = miRkwood->CONFIG();
    #~ debug('     - Process miRNA candidates' . ' [' . localtime() . ']', miRkwood->DEBUG() );

    foreach my $key ( sort keys %candidates_hash ) {
        $candidate_identifier++;
        debug( "         Start run_pipeline_on_candidate on $candidate_identifier" . ' [' . localtime() . ']', miRkwood->DEBUG());
        my $candidate = $candidates_hash{$key};
        my $final_candidate = $self->run_pipeline_on_candidate( $candidate_identifier, $candidate );
        if ( $cfg->param('options.filter_bad_hairpins') && $final_candidate->{'quality'} eq '0' && $final_candidate->{'alignment'} eq '0' ){
            print STDERR "Candidate $final_candidate->{'identifier'} is an orphan hairpin.\n";
        }
        else {
            push @candidates_result, $final_candidate;
        }
        debug( "         End of run_pipeline_on_candidate on $candidate_identifier" . ' [' . localtime() . ']', miRkwood->DEBUG());
    }

    miRkwood::Utils::display_var_sizes_in_log_file( '..... PrecursorBuilder : process_mirna_candidates' );

    return \@candidates_result;
}


sub run_pipeline_on_candidate {
    my ( $self, @args ) = @_;
    my $candidate_identifier = shift @args;
    my $candidate = shift @args;

    # Create candidate directory
    my $candidate_dir = $self->create_candidate_directory( $candidate, $candidate_identifier );

    # Run CandidateJob on the Candidate
    my $candidate_full_identifier = "$self->{'chromosome_id'}__$candidate->{'max'}{'start_position'}-$candidate->{'max'}{'end_position'}-$candidate_identifier";
    my $candidate_ref = $candidate->{'max'};
    my $alternatives = $candidate->{'alternatives'};
    my $candidatejob = miRkwood::CandidateJob->new($candidate_dir,
                                                   $self->{'chromosome_name'},
                                                   $candidate_full_identifier,
                                                   $candidate_ref,
                                                   $alternatives,
                                                   $self->{'genome_db'});

    miRkwood::Utils::display_var_sizes_in_log_file( '..... PrecursorBuilder : run_pipeline_on_candidate' );

    return $candidatejob->run();
}


sub create_candidate_directory {
    my ( $self, @args ) = @_;
    my $candidate = shift @args;
    my $candidate_identifier = shift @args;

    my $cluster_directory = miRkwood::Paths::get_workspace_candidate_dir( $self->{'workspace_dir'},
                                                                          $self->{'chromosome_id'},
                                                                          $candidate->{'max'}{'cluster'},
                                                                          $candidate->{'max'}{'strand'} ); # TO CHECK : ok for WebBAM but for genomic... ?

    my $candidate_dir = miRkwood::Paths::create_folder(File::Spec->catdir( $cluster_directory, $candidate_identifier ) );
    return $candidate_dir;
}

1;

