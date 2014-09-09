package miRkwood::CandidateJob;

# ABSTRACT: Job processing a candidate

use strict;
use warnings;

use miRkwood;
use miRkwood::CandidateHandler;
use miRkwood::Components;
use miRkwood::MiRdup;
use miRkwood::PosterioriTests;

use YAML::XS;
use Log::Message::Simple qw[msg error debug];

=method new

Constructor

Usage:
  my $candidate_job = miRkwood::CandidateJob->new($candidate_dir);

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($directory, $identifier, $candidate_ref, $alternatives) = @args;
    my $self = {
        directory => $directory,
        identifier => $identifier,
        candidate => $candidate_ref,
        alternatives => $alternatives,
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
    $self->populate_candidate_directory();
    $self->process_tests_for_candidate();
    return $self->make_candidate();
}

=method make_candidate


=cut

sub make_candidate {
    my ( $self, @args ) = @_;
    my $candidate_object = $self->get_candidate_information_from_run();
    return $self->update_candidate_information_from_self($candidate_object);
}

sub update_candidate_information_from_self {
    my ( $self, @args ) = @_;
    my $candidate_object = shift @args;
    $candidate_object->{'identifier'} = $self->{'identifier'};
    return $candidate_object;
}

=method populate_candidate_directory

Populate a candidate directory with the sequence, strand & so on.

=cut

sub populate_candidate_directory {
    my ( $self, @args ) = @_;
    my %candidate          = %{ $self->{'candidate'} };
    my @alternatives_array = @{ $self->{'alternatives'} };

    #Writing seq.txt
    my $candidate_sequence = File::Spec->catfile( $self->get_directory(), 'seq.txt' );
    open( my $SEQ_FH, '>', $candidate_sequence )
      or die "Error when opening $candidate_sequence: $!";
    print {$SEQ_FH} ">$candidate{'name'}\n$candidate{'dna'}\n";
    close $SEQ_FH;

    #Writing sequence information
    my $seq_info_file = File::Spec->catfile( $self->get_directory(), 'sequence_information.txt' );
    open( my $SEQ_INFO_FH, '>', $seq_info_file )
      or die "Error when opening $seq_info_file: $!";
    print {$SEQ_INFO_FH} $candidate{'strand'} . "\t" .  $candidate{'start'} . "\t" . $candidate{'end'};
    close $SEQ_INFO_FH;

    $self->process_outRNAFold( 'optimal', $candidate{'name'},
        $candidate{'dna'}, $candidate{'structure_optimal'}, $candidate{'energy_optimal'} );
    $self->process_outRNAFold( 'stemloop', $candidate{'name'},
        $candidate{'dna'}, $candidate{'structure_stemloop'}, $candidate{'energy_stemloop'} );

    # Writing energy file
    my $energy_file = File::Spec->catfile( $self->get_directory(), 'outMFEI.txt' );
    open( my $ENERGY_FH, '>', $energy_file )
        or die "Unable to open $energy_file: $!";
    my $content = $candidate{'name'} . "\t" . $candidate{'mfei'} . "\t" . $candidate{'energy_optimal'} . "\t" . $candidate{'amfe'};
    print $ENERGY_FH $content;
    close $ENERGY_FH or die "Unable to close: $!";

    #Writing alternativeCandidates.txt
    my $alternative_candidates =
      File::Spec->catfile( $self->get_directory(), 'alternativeCandidates.txt' );
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
        my $varna_image = File::Spec->catfile( $self->get_directory(), 'image.png' );
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
    my ( $self, @args ) = @_;
    my ( $suffix, $nameSeq, $dna, $structure, $energy ) = @args;

    my $candidate_rnafold_output =
      File::Spec->catfile( $self->get_directory(), "outRNAFold_$suffix.txt" );

    open( my $OUT2, '>', $candidate_rnafold_output )
      or die "Error when opening $candidate_rnafold_output: $!";
    print {$OUT2} ">$nameSeq\n$dna\n$structure ($energy)\n";
    close $OUT2;
    return;
}


=method get_candidate_information_from_run

Arguments:
- $job_dir - the job directory
- $seq_dir - the sequence directory
- $can_dir - the candidate directory

=cut

sub get_candidate_information_from_run {
    my ( $self, @args ) = @_;

    my $cfg    = miRkwood->CONFIG();
    my $candidate = $self->make_candidate_from_directory();

    if ( $cfg->param('options.varna') ) {
        $candidate->{'image'} = File::Spec->catfile($self->get_directory(), 'image.png');
    } else {
        $candidate->{'image'} = '';
    }

    $candidate->{'position'} = "$candidate->{'start_position'}-$candidate->{'end_position'}";
    $candidate->{'length'} = $candidate->{'end_position'} - $candidate->{'start_position'} +1;
    $candidate->{'%GC'} = miRkwood::Utils::restrict_num_decimal_digits(
                            miRkwood::Utils::compute_gc_content($candidate->{'DNASequence'}),
                            3);

    my $alternative_candidates_file = File::Spec->catfile($self->get_directory(), 'alternativeCandidates.txt');
    if (-e $alternative_candidates_file){
        my %alternatives = miRkwood::Parsers::parse_alternative_candidates_file($alternative_candidates_file);
        $candidate->{'alternatives'} = \%alternatives;
    }

    my $hairpin = miRkwood::Utils::make_ASCII_viz($candidate->{'DNASequence'}, $candidate->{'Vienna'});
    $candidate->{'hairpin'} = $hairpin;
#    my %sequence;
#    $sequence{$candidate{'name'}} = $candidate{'DNASequence'};
#    my $tmp_file = File::Spec->catfile($full_candidate_dir, "mirdup_prediction.txt");
#    $candidate{'mirdup_prediction'} = \%{miRkwood::MiRdup->predict_with_mirdup($tmp_file, \%sequence)};

    return $candidate;
}

=method parse_candidate_information

Get the results for a given candidate

Arguments:
- $full_candidate_dir - the prefixed path to the candidate results

=cut

sub parse_candidate_information {
    my ( $self, @args ) = @_;

    my %result = ();

    my $seq_info_file =
      File::Spec->catfile( $self->get_directory(), 'sequence_information.txt' );
    if ( -e $seq_info_file )    # si fichier existe
    {
        my @res = miRkwood::Components::get_sequence_information($seq_info_file);
        ($result{'strand'}, $result{'start_position'}, $result{'end_position'}) = @res;
    }

    my $randfold_output =
      File::Spec->catfile( $self->get_directory(), 'randfold.out' );
    if ( -e $randfold_output )    # si fichier existe
    {
        $result{'shuffles'} = miRkwood::Parsers::parse_pvalue($randfold_output);
    }

    #Récupération valeur MFEI
    my $mfei_out =
      File::Spec->catfile( $self->get_directory(), 'outMFEI.txt' );
    if ( -e $mfei_out )                 # si fichier existe
    {
        my @mfeis = miRkwood::Parsers::parse_mfei($mfei_out);
        $result{'mfei'} = $mfeis[0];
        $result{'mfe'} = $mfeis[1];
        $result{'amfe'} = $mfeis[2];
    }

    #Récupération séquence et format Vienna
    my $rnafold_stemloop_out = File::Spec->catfile( $self->get_directory(),
                                       'outRNAFold_stemloop.txt' );
    if ( -e $rnafold_stemloop_out )                  # si fichier existe
    {
        my @res = miRkwood::Components::get_data_from_rnafold_out($rnafold_stemloop_out);
        my $devnull;
        ($result{'name'}, $devnull, $result{'DNASequence'}, $result{'Vienna'}) = @res;
    }

    #Récupération séquence et format Vienna
    my $rnafold_optimal_out = File::Spec->catfile( $self->get_directory(),
                                                   'outRNAFold_optimal.txt' );
    if ( -e $rnafold_optimal_out )                  # si fichier existe
    {
        my @vienna_res = miRkwood::Parsers::parse_RNAfold_output($rnafold_optimal_out);

        $result{'Vienna_optimal'} = $vienna_res[2];
    }

    #Récupération alignement avec mirBase
    my $alignments_results_file = File::Spec->catfile($self->get_directory(), 'merged_alignments.yml');
    my $mirdup_results_file = File::Spec->catfile($self->get_directory(), 'mirdup_results.yml');
    $result{'alignment_existence'} = ( -e $alignments_results_file && ! -z $alignments_results_file );
    if ($result{'alignment_existence'}){
        my %mirdup_results = YAML::XS::LoadFile($mirdup_results_file) or die("Error when parsing YAML file $mirdup_results_file");
        my %alignments = YAML::XS::LoadFile($alignments_results_file);
        $result{'alignments'} = \%alignments;
        $result{'mirdup_validation'} = \%mirdup_results;
    }

    return \%result;
}

=method make_candidate_from_directory


Arguments:
- $full_candidate_dir - the prefixed path to the candidate results

=cut

sub make_candidate_from_directory {
    my ( $self, @args ) = @_;
    my $candidate_information = $self->parse_candidate_information($self->get_directory());
    my $candidate = miRkwood::Candidate->new($candidate_information);
    $candidate->compute_alignment_quality();
    $candidate->compute_quality();
    return $candidate;
}



=method process_tests_for_candidate

Perform the a posteriori tests for a given candidate

=cut

sub process_tests_for_candidate {
    my ( $self, @args ) = @_;
    my $dir = $self->get_directory();

    my $cfg = miRkwood->CONFIG();

    ####calcul p-value randfold
    if ( $cfg->param('options.randfold') ) {
        my $seq_file = File::Spec->catfile( $self->get_directory(), 'seq.txt' );
        debug( "Running test_randfold on $seq_file", miRkwood->DEBUG() );
        miRkwood::PosterioriTests::test_randfold( $self->get_directory(),
            $seq_file );
    }

    if ( $cfg->param('options.align') ) {
        my $candidate_rnafold_stemploop_out =
            File::Spec->catfile( $self->get_directory(), 'outRNAFold_stemloop.txt' );
        debug( "Running test_alignment on $candidate_rnafold_stemploop_out", miRkwood->DEBUG() );
        my $file_alignement =
          miRkwood::PosterioriTests::test_alignment( $self->get_directory(),
            $candidate_rnafold_stemploop_out );
        $self->post_process_alignments($candidate_rnafold_stemploop_out, $file_alignement );
    }

    return;
}

=method post_process_alignments


=cut

sub post_process_alignments {
    my ( $self, @args ) = @_;
    my $candidate_rnafold_stemploop_out = shift @args;
    my $file_alignement                 = shift @args;

    my @res =
      miRkwood::Components::get_data_from_rnafold_out(
        $candidate_rnafold_stemploop_out);
    my ( $name, $position, $DNASequence, $Vienna ) = @res;
    my %alignments;

    if (-z $file_alignement){
        return;
    }
    if (
        !eval {
            %alignments =
              miRkwood::Components::parse_custom_exonerate_output(
                $file_alignement);
        }
      )
    {
        # Catching exception
        carp("Exception when parsing exonerate output $file_alignement");
        return;
    }
    else {
        %alignments =
          miRkwood::Components::merge_alignments( \%alignments );
        my $tmp_file =
          File::Spec->catfile( $self->get_directory(), "mirdup_validation.txt" );
        my %mirdup_results =
          miRkwood::MiRdup->validate_with_mirdup( $tmp_file, $name,
            $DNASequence, $Vienna, keys %alignments );
        my $mirdup_results_file =
          File::Spec->catfile( $self->get_directory(), 'mirdup_results.yml' );
        YAML::XS::DumpFile( $mirdup_results_file, %mirdup_results );

        my $alignments_results_file =
          File::Spec->catfile( $self->get_directory(), 'merged_alignments.yml' );
        YAML::XS::DumpFile( $alignments_results_file, %alignments );
    }
}

1;
