package miRkwood::CandidateJob;

# ABSTRACT: Job processing a candidate

use strict;
use warnings;

use miRkwood;
use miRkwood::Candidate;
use miRkwood::PosterioriTests;
use miRkwood::Programs;
use miRkwood::Utils;
use miRkwood::MultiAlignments;

use Carp;
use Log::Message::Simple qw[msg error debug];

=method new

Constructor

Usage:
  my $candidate_job = miRkwood::CandidateJob->new($candidate_dir);

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($directory, $sequence_name, $identifier, $candidate_ref, $alternatives, $genome_db) = @args;
    if ( ! defined( $genome_db ) ){
        $genome_db = '';
    }
    my $self = {
        'sequence_name' => $sequence_name,
        'directory'     => $directory,
        'identifier'    => $identifier,
        'candidate'     => $candidate_ref,
        'alternatives'  => $alternatives,
        'genome_db'     => $genome_db,
    };
    bless $self, $class;
    return $self;
}

sub get_directory {
    my ( $self, @args ) = @_;
    return $self->{'directory'};
}

sub run{
    my ( $self, @args ) = @_;
    my $cfg = miRkwood->CONFIG();

    if ( $cfg->param('job.pipeline') eq 'smallRNAseq' ) {
        $self->{'candidate'} = miRkwood::Candidate::count_total_nb_of_reads_for_candidate( $self->{'candidate'} );
        $self->{'candidate'} = miRkwood::Candidate::find_mirna( $self->{'candidate'}, $self->{'genome_db'} );
    }

    my $candidate_test_info = $self->process_tests_for_candidate();

    my $candidate_information = $self->get_candidate_information();

    my %complete_candidate = (%{$candidate_test_info}, %{$candidate_information});
    my $candidate = miRkwood::Candidate->new(\%complete_candidate);

    if ( $cfg->param('job.pipeline') eq 'smallRNAseq' ) {
        $candidate->compute_alignment_quality_for_smallRNAseq();
        $candidate->compute_quality();
        $candidate->tag_orphan_hairpins();
    }
    else{
        $candidate->compute_alignment_quality_for_abinitio();
        $candidate->compute_quality();
    }
    return $candidate;
}

=method get_candidate_information


=cut

sub get_candidate_information {
    my ( $self, @args ) = @_;
    my $candidate_structure = $self->{'candidate'};

    $candidate_structure = $self->update_candidate_information_from_self($candidate_structure);
    $candidate_structure = $self->update_candidate_information_from_run($candidate_structure);

    return $candidate_structure;
}

sub update_candidate_information_from_self {
    my ( $self, @args ) = @_;
    my $candidate_structure = shift @args;
    $candidate_structure->{'name'} = $self->{'sequence_name'};
    $candidate_structure->{'identifier'} = $self->{'identifier'};
    $candidate_structure->{'alternatives'} = $self->convert_alternative_candidates();
    return $candidate_structure;
}

sub update_candidate_information_from_run{
    my ( $self, @args ) = @_;
    my $candidate_structure = shift @args;
    $candidate_structure->{'mfe'} = delete $candidate_structure->{'energy_stemloop'};
    $candidate_structure->{'image'} = $self->write_VARNA_if_needed();

    $candidate_structure->{'position'} = $self->make_full_position($candidate_structure);
    $candidate_structure->{'length'} = $candidate_structure->{'end_position'} - $candidate_structure->{'start_position'} + 1;
    $candidate_structure->{'%GC'} = miRkwood::Utils::restrict_num_decimal_digits(
                            miRkwood::Utils::compute_gc_content($candidate_structure->{'sequence'}),
                            3);

    my $hairpin = miRkwood::Utils::make_ASCII_viz($candidate_structure->{'sequence'}, $candidate_structure->{'structure_stemloop'});
    $candidate_structure->{'hairpin'} = $hairpin;
    return $candidate_structure;
}

sub update_known_candidate_information {
    my ( $self, @args ) = @_;
    my $candidate = shift @args;
    my $genome = shift @args;

    my $start = $candidate->{'start_position'};
    my $end   = $candidate->{'end_position'};

    $candidate->count_total_nb_of_reads_for_candidate();
    $candidate->compute_quality_for_known_miRNAs();
    $candidate->find_mirna_for_known_candidate( $genome );
    $candidate->{'sequence'} = $genome->seq( $candidate->{'chromosome'}, $start => $end );
    if ( $candidate->{'strand'} eq '-' ){
        $candidate->{'sequence'} = miRkwood::Utils::reverse_complement( $candidate->{'sequence'} );
    }
    $candidate->{'sequence'} =~ s/T/U/g;
    $candidate->{'%GC'} = miRkwood::Utils::restrict_num_decimal_digits(
                             miRkwood::Utils::compute_gc_content($candidate->{'sequence'}), 3);
    ($candidate->{'structure_optimal'}, $candidate->{'mfe'}) = $self->get_stemloop_structure_for_known_candidate();
    $candidate->{'structure_stemloop'} = get_stemloop_structure_from_optimal( $candidate->{'structure_optimal'} );
    ($candidate->{'mfei'}, $candidate->{'amfe'}) = miRkwood::Utils::compute_mfei_and_amfe( $candidate->{'sequence'}, $candidate->{'mfe'} );
    $candidate->{'mfei'} = miRkwood::Utils::restrict_num_decimal_digits( $candidate->{'mfei'}, 3);
    $candidate->{'amfe'} = miRkwood::Utils::restrict_num_decimal_digits( $candidate->{'amfe'}, 3);
    $candidate->{'image'} = $self->write_VARNA_if_needed();

    return $candidate;
}

=method get_stemloop_structure_from_optimal

  From the dot-bracket optimal structure, 
  calculates the stemloop structure
  
=cut
sub get_stemloop_structure_from_optimal {
    my (@args) = @_;
    my $optimal = shift @args;
    my @optimal = split( //, $optimal );
    my $length_structure = length( $optimal );

    my $stemloop = '';
    my $nb_open_brackets = 0;
    my $nb_closed_brackets = 0;
    my $i = 0;
    my $continue = 1;
    my $found_close = 0;

    while ( $i < $length_structure and $continue ){
        if ( $optimal[$i] eq '(' ){
            if ( !$found_close ){
                $nb_open_brackets++;
            }
            else{
                $continue = 0;
            }
        }
        elsif ( $optimal[$i] eq ')' ){
            if ( !$found_close ){
                $found_close = 1;
            }
            $nb_closed_brackets++;
        }
        $i++;
    }

    if ( $nb_open_brackets == $nb_closed_brackets ){
        $stemloop = $optimal;
    }
    else{
        my $nb_brackets_in_stemloop = $nb_open_brackets - $nb_closed_brackets;
        my $nb_characters = 0;
        my $nb_brackets = 0;
        while ( $nb_brackets < $nb_brackets_in_stemloop) {
            $stemloop .= $optimal[$nb_characters];
            if ( $optimal[$nb_characters] eq '(' ){
                $nb_brackets++;
            }
            $nb_characters++;
        }

        $nb_characters = scalar(@optimal) - 1;
        $nb_brackets = 0;
        my $right_arm = '';
        while ( $nb_brackets < $nb_brackets_in_stemloop) {
            $right_arm .= $optimal[$nb_characters];
            if ( $optimal[$nb_characters] eq ')' ){
                $nb_brackets++;
            }
            $nb_characters--;
        }

        my $nb_dots_in_middle = length($optimal) - length($stemloop) - length($right_arm);
        for (my $i = 0; $i < $nb_dots_in_middle; $i++){
            $stemloop .= '.';
        }
        $stemloop .= reverse $right_arm;

    }

    return $stemloop;

}

sub make_full_position {
    my ( $self, @args ) = @_;
    my $candidate_structure = shift @args;
    return "$candidate_structure->{'start_position'}-$candidate_structure->{'end_position'}";
}

=method write_VARNA_if_needed


=cut

sub write_VARNA_if_needed {
    my ( $self, @args ) = @_;
    my $cfg = miRkwood->CONFIG();
    if ( $cfg->param('options.varna') ) {
        my $varna_image = File::Spec->catfile( miRkwood::Paths::get_dir_images_path_from_job_dir( $cfg->param('job.directory') ),
                                                $self->{'identifier'} . '.png' );
        debug( "           Running VARNA in $varna_image" . ' [' . localtime() . ']', miRkwood->DEBUG() );
        miRkwood::Programs::run_varna_on_structure( $self->{'candidate'}{'sequence'}, $self->{'candidate'}{'structure_stemloop'}, $varna_image )
          or carp('Problem during image generation using VARNA');
      return $varna_image
    }
    return '';
}

sub get_stemloop_structure_for_known_candidate {
    my ( $self, @args ) = @_;

    #~ debug( '        Running RNAfold', miRkwood->DEBUG() );
    my $sequence_name = $self->{'candidate'}{'name'};
    my $sequence = $self->{'candidate'}{'sequence'};

    my $rnafold_output = File::Spec->catfile( $self->get_directory(), 'RNAfold.out' );
    my $temp_file = File::Spec->catfile( $self->get_directory(), 'tempFile.txt' );

    miRkwood::Programs::run_rnafold( $sequence_name, $sequence, $temp_file, $rnafold_output )
      or die("Problem when running RNAfold: $!");

    my $nameSeq = '';
    my $dna = '';
    my $stemloop_structure = '';
    my $energy = 0;

    ( $nameSeq, $dna, $stemloop_structure, $energy ) = miRkwood::Parsers::parse_RNAfold_output( $rnafold_output );

    return ($stemloop_structure, $energy);
}

=method process_tests_for_candidate

Perform the a posteriori tests for a given candidate

=cut

sub process_tests_for_candidate {
    my ( $self, @args ) = @_;

    my $cfg = miRkwood->CONFIG();
    my $result = {};
    my $mirdup_results;

    my $posteriori_tests = miRkwood::PosterioriTests->new($self->get_directory());

    if ( $cfg->param('options.randfold') ) {
        my $seq_file = $self->write_sequence_file();
        debug( "              Running test_randfold on $seq_file", miRkwood->DEBUG() );
        $result->{'shuffles'} = $posteriori_tests->test_randfold( $seq_file );
    }

    $posteriori_tests->{'sequence_name'} = $self->{'sequence_name'};
    $posteriori_tests->{'candidate'} = $self->{'candidate'};
    my $candidate_rnafold_stemploop_out = $self->write_RNAFold_stemloop_output();

    my $alignments = $posteriori_tests->test_alignment( $candidate_rnafold_stemploop_out );

    miRkwood::Candidate::store_attribute_ct( $self->{'candidate'}, $self->{'directory'} );

    if ( $cfg->param('options.align') ) {
        if ( scalar(keys%{$alignments}) > 0 ){
            miRkwood::MultiAlignments::fillTabTemp2D($alignments, $self->{'identifier'});
        }

        if ( $cfg->param('job.pipeline') eq 'smallRNAseq' ){
            $result->{'criteria_mirdup'} = 0;
            if ( defined( $self->{'candidate'}{'mirna_sequence'} ) and $self->{'candidate'}{'mirna_sequence'} ne '' ){
                my %mirdup_on_mirna = $posteriori_tests->validate_mirna_with_mirdup();
                $result->{'criteria_mirdup'} = $mirdup_on_mirna{'mirna'};
            }
            foreach (keys%{$alignments}){
                $alignments->{$_} = '';
            }
        }
        else {
            $mirdup_results = $posteriori_tests->validate_alignments_with_mirdup( $alignments );
        }

        if ( scalar(keys%{$alignments}) > 0) {
            $result->{'alignment_existence'} = 1;
            $result->{'alignments'} = $alignments;
            if ( $cfg->param('job.pipeline') eq 'abinitio' ){
                $result->{'mirdup_validation'} = $mirdup_results;
            }
        }
    } else {
        $result->{'alignment_existence'} = 0;
    }

    return $result;
}

=method write_RNAFold_stemloop_output

Writing (pseudo) rnafold output

=cut

sub write_RNAFold_stemloop_output {
    my ( $self, @args ) = @_;
    my $candidate_rnafold_output =
      File::Spec->catfile( $self->get_directory(), 'outRNAFold_stemloop.txt' );
    open( my $OUT2, '>', $candidate_rnafold_output )
      or die "Error when opening $candidate_rnafold_output: $!";
    print {$OUT2} ">$self->{'candidate'}{'name'}\n$self->{'candidate'}{'sequence'}\n$self->{'candidate'}{'structure_stemloop'} ($self->{'candidate'}{'energy_stemloop'})\n";
    close $OUT2;
    return $candidate_rnafold_output;
}

=method write_sequence_file

Writing the candidate as a FASTA sequence

=cut

sub write_sequence_file {
    my ( $self, @args ) = @_;
    my $candidate_sequence = File::Spec->catfile( $self->get_directory(), 'seq.txt' );
    open( my $SEQ_FH, '>', $candidate_sequence )
        or die "Error when opening $candidate_sequence: $!";
     print {$SEQ_FH} ">$self->{'candidate'}{'name'}\n$self->{'candidate'}{'sequence'}\n";
    close $SEQ_FH;
    return $candidate_sequence;
}

=method convert_alternative_candidates

=cut

sub convert_alternative_candidates {
    my ( $self, @args ) = @_;
    my @alternatives = @{$self->{'alternatives'}};
    my %results;
    foreach my $alternative (@alternatives) {
        my $name = delete $alternative->{'name'};
        $name .= ' (MFEI: ' . miRkwood::Utils::restrict_num_decimal_digits( $alternative->{'mfei'}, 3) . ')';
        $results{$name}{'sequence'} = $alternative->{'sequence'};
        $results{$name}{'structure_stemloop'} = $alternative->{'structure_stemloop'};
    }
    return \%results;
}

1;
