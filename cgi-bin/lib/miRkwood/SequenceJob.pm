package miRkwood::SequenceJob;

# ABSTRACT: Job processing a given sequence

use strict;
use warnings;

use Log::Message::Simple qw[msg error debug];

use miRkwood::Programs;
use miRkwood::CandidateJob;
use miRkwood::SequenceSubJob;
use miRkwood::Utils;

=method new

Constructor

my $sequence_job = miRkwood::SequenceSubJob->new($sequence_dir, $name, $sequence, '+');

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($directory, $name, $sequence) = @args;
    my $self = {
        directory => $directory,
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
    my $candidates_array =
      $self->get_raw_candidates();

    my %candidates_hash = $self->process_raw_candidates($candidates_array);

    $self->create_directories( \%candidates_hash );
    return;
}

=method get_raw_candidates

Get the candidates for the sequence

 Usage : $self->get_raw_candidates();

=cut

sub get_raw_candidates{
    my ($self, @args) = @_;
    my $sequence_subjob = miRkwood::SequenceSubJob->new($self->get_directory(), $self->{'name'}, $self->{'sequence'}, '+');
    my $candidates = $sequence_subjob->get_raw_candidates_for_sequence();

    my @candidates_array1 = @{$candidates};
    my @candidates_array;

    my $cfg = miRkwood->CONFIG();
    if ( $cfg->param('options.strands') ) {
        debug( 'Processing the other strand', miRkwood->DEBUG() );
        my $reversed_sequence =
          miRkwood::Utils::reverse_complement($self->{'sequence'});
        my $sequence_subjob2 = miRkwood::SequenceSubJob->new($self->get_directory(), $self->{'name'}, $reversed_sequence, '-');
        my $candidates2 = $sequence_subjob2->get_raw_candidates_for_sequence();
        my @candidates_array2 = @{$candidates2};
        @candidates_array = sort { $a->{start} <=> $b->{start} } ( @candidates_array1, @candidates_array2 );
    }
    else {
        @candidates_array = sort { $a->{start} <=> $b->{start} } ( @candidates_array1 );
    }
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
        my $start = $current_candidate{'start'};
        my $end   = $current_candidate{'end'};
        my ( $ref_start, $ref_end ) =
          ( $reference_candidate{'start'}, $reference_candidate{'end'} );
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


=method create_directories

Create the necessary directories.

=cut

sub create_directories {
    my ($self, @args) = @_;
    my (%candidates_hash)    = %{ shift @args };

    my $candidate_counter    = 0;
    foreach my $key ( sort keys %candidates_hash ) {
        $candidate_counter++;
        my $candidate_dir =
          File::Spec->catdir( $self->get_directory(), $candidate_counter );
        mkdir $candidate_dir;
        my $candidate_ref = $candidates_hash{$key}{'max'};
        populate_candidate_directory( $candidate_dir, $candidate_ref,
            $candidates_hash{$key}{'alternatives'} );
    }
    return;
}

=method populate_candidate_directory

Populate a candidate directory with the sequence, strand & so on.

=cut

sub populate_candidate_directory {
    my @args               = @_;
    my $candidate_dir      = shift @args;
    my %candidate          = %{ shift @args };
    my @alternatives_array = @{ shift @args };

    debug( "Writing candidate information in $candidate_dir", miRkwood->DEBUG() );
    #Writing seq.txt
    my $candidate_sequence = File::Spec->catfile( $candidate_dir, 'seq.txt' );
    open( my $SEQ_FH, '>', $candidate_sequence )
      or die "Error when opening $candidate_sequence: $!";
    print {$SEQ_FH} ">$candidate{'name'}\n$candidate{'dna'}\n";
    close $SEQ_FH;

    #Writing sequence information
    my $seq_info_file = File::Spec->catfile( $candidate_dir, 'sequence_information.txt' );
    open( my $SEQ_INFO_FH, '>', $seq_info_file )
      or die "Error when opening $seq_info_file: $!";
    print {$SEQ_INFO_FH} $candidate{'strand'} . "\t" .  $candidate{'start'} . "\t" . $candidate{'end'};
    close $SEQ_INFO_FH;

    process_outRNAFold( $candidate_dir, 'optimal', $candidate{'name'},
        $candidate{'dna'}, $candidate{'structure_optimal'}, $candidate{'energy_optimal'} );
    process_outRNAFold( $candidate_dir, 'stemloop', $candidate{'name'},
        $candidate{'dna'}, $candidate{'structure_stemloop'}, $candidate{'energy_stemloop'} );

    # Writing energy file
    my $energy_file = File::Spec->catfile( $candidate_dir, 'outMFEI.txt' );
    open( my $ENERGY_FH, '>', $energy_file )
        or die "Unable to open $energy_file: $!";
    my $content = $candidate{'name'} . "\t" . $candidate{'mfei'} . "\t" . $candidate{'energy_optimal'} . "\t" . $candidate{'amfe'};
    print $ENERGY_FH $content;
    close $ENERGY_FH or die "Unable to close: $!";

    #Writing alternativeCandidates.txt
    my $alternative_candidates =
      File::Spec->catfile( $candidate_dir, 'alternativeCandidates.txt' );
    if (@alternatives_array){
        open( my $OUT2, '>>', $alternative_candidates )
          or die "Error when opening $alternative_candidates: $!";
        foreach my $alternative (@alternatives_array) {
            print $OUT2
">$alternative->{'name'}\t$alternative->{'dna'}\t$alternative->{'structure_optimal'}\t$alternative->{'mfei'}\n";
        }
        close $OUT2;
    }

    # Writing VARNA image
    my $cfg = miRkwood->CONFIG();

    if ( $cfg->param('options.varna') ) {
        my $varna_image = File::Spec->catfile( $candidate_dir, 'image.png' );
        debug( "Generating image using VARNA in $varna_image", miRkwood->DEBUG() );
        miRkwood::Programs::run_varna_on_structure( $candidate{'dna'}, $candidate{'structure_stemloop'}, $varna_image )
          or carp('Problem during image generation using VARNA');
    }

    return;
}

=method process_outRNAFold

Writing (pseudo) rnafold output

=cut

sub process_outRNAFold {
    my ( $candidate_dir, $suffix, $nameSeq, $dna, $structure, $energy ) = @_;

    my $candidate_rnafold_output =
      File::Spec->catfile( $candidate_dir, "outRNAFold_$suffix.txt" );

    open( my $OUT2, '>', $candidate_rnafold_output )
      or die "Error when opening $candidate_rnafold_output: $!";
    print {$OUT2} ">$nameSeq\n$dna\n$structure ($energy)\n";
    close $OUT2;
    return;
}

1;