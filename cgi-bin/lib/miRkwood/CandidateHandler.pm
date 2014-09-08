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

=method serialize_candidate

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


=method serialize_candidate_from_run

Retrieve the information from a run and serialize a candidate on disk

Arguments:
- $job_dir - the job directory
- $seq_dir - the sequence directory
- $can_dir - the candidate directory
- $serialization_dir - The directory in which to serialize

=cut

sub serialize_candidate_from_run{
    my ( $self, @args ) = @_;
    my $job_dir = shift @args;
    my $seq_dir = shift @args;
    my $can_dir = shift @args;
    my $serialization_dir = shift @args;
    my $candidate_object = $self->get_candidate_information_from_run( $job_dir, $seq_dir, $can_dir );
    return $self->serialize_candidate_information($serialization_dir, $candidate_object);
}

=method get_candidate_information_from_run

Arguments:
- $job_dir - the job directory
- $seq_dir - the sequence directory
- $can_dir - the candidate directory

=cut

sub get_candidate_information_from_run {
    my ( $self, @args ) = @_;
    my $job_dir = shift @args;
    my $seq_dir = shift @args;
    my $can_dir = shift @args;

    my $cfg    = miRkwood->CONFIG();
    my $full_candidate_dir = miRkwood::Paths->get_candidate_paths($job_dir,  $seq_dir, $can_dir);

    my $candidate = $self->make_candidate_from_directory($full_candidate_dir);
    $candidate->{'identifier'} = "$seq_dir-$can_dir";
#    $candidate{'name'} = $seq_dir;    #récupération nom séquence

    if ( $cfg->param('options.varna') ) {
        $candidate->{'image'} = File::Spec->catfile($full_candidate_dir, 'image.png');
    } else {
        $candidate->{'image'} = '';
    }

    $candidate->{'position'} = "$candidate->{'start_position'}-$candidate->{'end_position'}";
    $candidate->{'length'} = $candidate->{'end_position'} - $candidate->{'start_position'} +1;
    $candidate->{'%GC'} = miRkwood::Utils::restrict_num_decimal_digits(
                            miRkwood::Utils::compute_gc_content($candidate->{'DNASequence'}),
                            3);

    my $alternative_candidates_file = File::Spec->catfile($full_candidate_dir, 'alternativeCandidates.txt');
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
    my $full_candidate_dir = shift @args;

    my %result = ();

    my $seq_info_file =
      File::Spec->catfile( $full_candidate_dir, 'sequence_information.txt' );
    if ( -e $seq_info_file )    # si fichier existe
    {
        my @res = miRkwood::Components::get_sequence_information($seq_info_file);
        ($result{'strand'}, $result{'start_position'}, $result{'end_position'}) = @res;
    }

    my $randfold_output =
      File::Spec->catfile( $full_candidate_dir, 'randfold.out' );
    if ( -e $randfold_output )    # si fichier existe
    {
        $result{'shuffles'} = miRkwood::Parsers::parse_pvalue($randfold_output);
    }

    #Récupération valeur MFEI
    my $mfei_out =
      File::Spec->catfile( $full_candidate_dir, 'outMFEI.txt' );
    if ( -e $mfei_out )                 # si fichier existe
    {
        my @mfeis = miRkwood::Parsers::parse_mfei($mfei_out);
        $result{'mfei'} = $mfeis[0];
        $result{'mfe'} = $mfeis[1];
        $result{'amfe'} = $mfeis[2];
    }

    #Récupération séquence et format Vienna
    my $rnafold_stemloop_out = File::Spec->catfile( $full_candidate_dir,
                                       'outRNAFold_stemloop.txt' );
    if ( -e $rnafold_stemloop_out )                  # si fichier existe
    {
        my @res = miRkwood::Components::get_data_from_rnafold_out($rnafold_stemloop_out);
        my $devnull;
        ($result{'name'}, $devnull, $result{'DNASequence'}, $result{'Vienna'}) = @res;
    }

    #Récupération séquence et format Vienna
    my $rnafold_optimal_out = File::Spec->catfile( $full_candidate_dir,
                                                   'outRNAFold_optimal.txt' );
    if ( -e $rnafold_optimal_out )                  # si fichier existe
    {
        my @vienna_res = miRkwood::Parsers::parse_RNAfold_output($rnafold_optimal_out);

        $result{'Vienna_optimal'} = $vienna_res[2];
    }

    #Récupération alignement avec mirBase
    my $alignments_results_file = File::Spec->catfile($full_candidate_dir, 'merged_alignments.yml');
    my $mirdup_results_file = File::Spec->catfile($full_candidate_dir, 'mirdup_results.yml');
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
    my $full_candidate_dir = shift @args;
    my $candidate_information = $self->parse_candidate_information($full_candidate_dir);
    my $candidate = miRkwood::Candidate->new($candidate_information);
    $candidate->compute_alignment_quality();
    $candidate->compute_quality();
    return $candidate;
}

1;