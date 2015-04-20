package miRkwood::BEDPipeline;

# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use parent 'miRkwood::Pipeline';

use Log::Message::Simple qw[msg error debug];

use miRkwood::Utils;
use miRkwood::ClusterBuilder;
use miRkwood::HairpinBuilder;
use miRkwood::PrecursorBuilder;

use Bio::DB::Fasta;

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($job_dir, $bed_file, $genome_file) = @args;
    #~ my %genome_db = miRkwood::Utils::multifasta_to_hash( $genome_file );
    my $self = {
        job_dir     => $job_dir,
        initial_bed => $bed_file,       # this is for the non filtered BED provided by the user
        bed_file    => '',              # this is for the final BED after all filtrations have been done
        mirna_bed   => '',              # this is for known miRNAs (ie miRNAs from miRBase). created during the filtration step. Maybe change that name ?
        genome_file => $genome_file,
        genome_db   => Bio::DB::Fasta->new($genome_file),
        sequences   => undef
    };
    bless $self, $class;
    return $self;
}


sub run_pipeline {
    my ($self, @args) = @_;

    $self->init_pipeline();

    debug( 'miRkwood : start processing.' . ' [' . gmtime() . ']', miRkwood->DEBUG() );

    $self->calculate_reads_coverage();

    $self->filter_BED();

    # Look for known miRNAs
    if ( $self->{'mirna_bed'} ne '' ){
        debug( 'Treat known miRNAs.' . ' [' . gmtime() . ']', miRkwood->DEBUG() );
        $self->treat_known_mirnas();
    }
    else{
        debug( 'No BED for known miRNAs.', miRkwood->DEBUG() );
    }

    # Look for new miRNAs
    debug( 'Treat new miRNAs.' . ' [' . gmtime() . ']', miRkwood->DEBUG() );
    
# Start new way of doing (chromosome by chromosome)
# Comment out the following lines to go back to the previous version    
    $self->list_chrom_in_bed();
    debug( scalar(@{$self->{'chromosomes_in_bed'}}) . ' chromosome(s) to consider', miRkwood->DEBUG() );
    foreach my $chromosome ( @{$self->{'chromosomes_in_bed'}} ){
        debug( "- Considering chromosome $chromosome" . ' [' . gmtime() . ']', miRkwood->DEBUG() );
        $self->init_sequences_per_chr( $chromosome );
        $self->run_pipeline_on_sequences_per_chr( $chromosome );
    }
# End new way of doing

# Un-comment out the two following lines to go back to previous version
    #~ $self->init_sequences();
    #~ $self->run_pipeline_on_sequences();

    $self->mark_job_as_finished();

    debug( 'Writing finish file' . ' [' . gmtime() . ']', miRkwood->DEBUG() );

    return;
}


=method create_additional_directories

=cut
sub create_additional_directories {
    my ($self, @args) = @_;
    mkdir $self->get_new_candidates_dir();
    mkdir $self->get_known_candidates_dir();
    mkdir $self->get_reads_dir();
    mkdir $self->get_known_reads_dir();
    mkdir $self->get_new_reads_dir();
    return;
}


=method calculate_reads_coverage

  On the BED given by the user, calculate the read coverage, ie the 
  number X such as on average one read starts every X nt.

=cut
sub calculate_reads_coverage {
    my ($self, @args) = @_;
    my $line;
    my $size_genome = 0;
    my $nb_tot_reads = 0;
    my $chromosome = '';
    my $end_position = 0;
    open (my $BED, $self->{'initial_bed'}) or die "ERROR while opening $self->{'initial_bed'} : $!";
    while( <$BED> ){
        chomp;
        if ( $_ ne '' ){
            my @fields = split( /\t/);
            $nb_tot_reads += $fields[4];
            if ( $chromosome ne '' and $chromosome ne $fields[0] ){
                $size_genome += $end_position;
            }
            $chromosome = $fields[0];
            $end_position = $fields[2];
            $line = $_;
        }
    }
    close $BED;
    chomp $line;
    my @fields = split( /\t/, $line);
    $size_genome += $fields[2];
    $self->{'average_coverage'} = int( $size_genome / $nb_tot_reads );
    debug( "Average reads coverage for this BED file : 1 read every $self->{'average_coverage'} nt.", miRkwood->DEBUG() );
    return;
}


=method filter_BED

Method to call module BEDHandler.pm and filter the given BED file
of CDS, other RNA, multimapped reads, and known miRNAs.
Initialize attribute 'bed_file' with the filtered BED, or with the
initial BED file if no filter had been done.

=cut

sub filter_BED {
    my ($self, @args) = @_;
    my $cfg                = miRkwood->CONFIG();
    my $species            = $cfg->param('job.plant');
    my $filter_CDS         = $cfg->param('options.filter_CDS');
    my $filter_tRNA_rRNA   = $cfg->param('options.filter_tRNA_rRNA');
    my $filter_multimapped = $cfg->param('options.filter_multimapped');
    my $localBED           = $self->{'initial_bed'};
    my $filteredBED        = '';
    my $mirnaBED           = '';

    if ( $species ne '' ){
        ($filteredBED, $mirnaBED) = miRkwood::BEDHandler->filterBEDfile_for_model_organism( $self->get_job_dir(),
                                                                                            $localBED,
                                                                                            $species,
                                                                                            $filter_CDS,
                                                                                            $filter_tRNA_rRNA,
                                                                                            $filter_multimapped );
    }
    else{
        ($filteredBED, $mirnaBED) = miRkwood::BEDHandler->filterBEDfile_for_user_sequence( $self->get_job_dir(),
                                                                                           $localBED,
                                                                                           $filter_CDS,
                                                                                           $filter_tRNA_rRNA,
                                                                                           $filter_multimapped );
    }

    if ( ! defined($filteredBED) or $filteredBED eq '' ){
        $self->{'bed_file'} = $self->{'initial_bed'};
        $self->{'mirna_bed'} = '';
    }
    else{
        $self->{'bed_file'} = $filteredBED;
        $self->{'mirna_bed'} = $mirnaBED;
    }

    return;
}


=method list_chrom_in_bed

  Read the BED file and store the name of all chromosomes
  
=cut
sub list_chrom_in_bed {
    my ($self) = @_;

    my %list_chromosomes;
    open ( my $BED, '<', $self->{'bed_file'} ) or die "ERROR while opening $self->{'bed_file'} : $!";
    while ( <$BED> ){
        my @fields = split( /\t/ );
        $list_chromosomes{ $fields[0] } = 1;
    }
    close $BED;
    my @sorted_list = sort( keys( %list_chromosomes ) );
    $self->{'chromosomes_in_bed'} = \@sorted_list;

}


=method init_sequences

=cut
# SEB BEGIN
sub init_sequences {
    my ($self, @args) = @_;
    debug( 'Extracting sequences from genome using BED clusters', miRkwood->DEBUG() );
    my $clustering = miRkwood::ClusterBuilder->new($self->{'genome_db'}, $self->{'bed_file'});   
    $self->{'sequences'} = $clustering->build_loci( $self->{'average_coverage'} );
    $self->{'parsed_reads'} = $clustering->get_parsed_bed();
    miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : init_sequences()');
    return;
}
# SEB END


=method init_sequences_per_chr

=cut
# SEB BEGIN
sub init_sequences_per_chr {
    my ($self, @args) = @_;
    my $chromosome = shift @args;
    debug( '   Constructing loci' . ' [' . gmtime() . ']', miRkwood->DEBUG() );
    my $clustering = miRkwood::ClusterBuilder->new($self->{'genome_db'}, $self->{'bed_file'});
    $self->{'sequences'} = $clustering->build_loci_per_chr( $chromosome, $self->{'average_coverage'} );
    $self->{'parsed_reads'} = $clustering->get_parsed_bed();
    miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : init_sequences_per_chr()');
    return;
}
# SEB END


sub run_pipeline_on_sequences_per_chr {
    my ($self, @args) = @_;
    my $chromosome = shift @args;

    $self->{'basic_candidates'} = [];

    my $sequences_count = scalar( @{$self->{'sequences'}} );

    debug( "   $sequences_count sequences to process", miRkwood->DEBUG() );

    $self->compute_candidates_per_chr( $chromosome );

    $self->serialize_basic_candidates( 'basic_candidates' );

    return;
}


=method compute_candidates

=cut

sub compute_candidates {
    my ($self, @args) = @_;

    my $sequence_identifier = 0;

    # Look for new miRNAs
    foreach my $chr (keys %{ $self->{'sequences'} }) {
        debug( "- Considering chromosome $chr", miRkwood->DEBUG() );
        my $hairpinBuilder = miRkwood::HairpinBuilder->new($self->{'genome_db'}, $self->get_workspace_path(), $self->{'parsed_reads'});
        my $loci_for_chr = $self->{'sequences'}->{$chr};
        my @hairpin_candidates_for_chr = ();
        foreach my $locus (@{$loci_for_chr}) {
            debug( "  - Considering sequence $sequence_identifier", miRkwood->DEBUG() );
            $sequence_identifier++;
            push @hairpin_candidates_for_chr, @{ $hairpinBuilder->build_hairpins($locus) };

            #~ miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : compute_candidates() (boucle sur chrom and sequences)' );  # /!\ WARNING : comment this line for web pipeline
        }
        undef $hairpinBuilder;

        my @sorted_hairpin_candidates_for_chr = sort { $a->{'start_position'} <=> $b->{'start_position'} } @hairpin_candidates_for_chr;

        undef @hairpin_candidates_for_chr;

        if ( scalar(@sorted_hairpin_candidates_for_chr) ){
            my $precursorBuilderJob = miRkwood::PrecursorBuilder->new( $self->get_workspace_path(), $chr, $chr );

            # Merge candidates
            my $candidates_hash = miRkwood::PrecursorBuilder::merge_candidates( \@sorted_hairpin_candidates_for_chr );

            # Posteriori tests and update candidate information
            my $final_candidates_hash = $precursorBuilderJob->process_mirna_candidates( $candidates_hash );

            undef $candidates_hash;

            $self->serialize_candidates($final_candidates_hash);

            miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : compute_candidates() (sort candidates)' );

            undef $final_candidates_hash;
        }
    }

    return;
}


sub compute_candidates_per_chr {
    my ($self, @args) = @_;
    my $chromosome = shift @args;

    my $sequence_identifier = 0;

    # Look for new miRNAs
    my $hairpinBuilder = miRkwood::HairpinBuilder->new($self->{'genome_db'}, $self->get_workspace_path(), $self->{'parsed_reads'});
    my @hairpin_candidates_for_chr = ();
    foreach my $locus ( @{$self->{'sequences'}} ) {
        debug( "  - Considering sequence $sequence_identifier" . ' [' . gmtime() . ']', miRkwood->DEBUG() );
        $sequence_identifier++;
        push @hairpin_candidates_for_chr, @{ $hairpinBuilder->build_hairpins($locus) };

        #~ miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : compute_candidates() (boucle sur chrom and sequences)' );  # /!\ WARNING : comment this line for web pipeline
    }
    undef $hairpinBuilder;

    my @sorted_hairpin_candidates_for_chr = sort { $a->{'start_position'} <=> $b->{'start_position'} } @hairpin_candidates_for_chr;

    undef @hairpin_candidates_for_chr;

    if ( scalar(@sorted_hairpin_candidates_for_chr) ){
        my $precursorBuilderJob = miRkwood::PrecursorBuilder->new( $self->get_workspace_path(), $chromosome, $chromosome );

        # Merge candidates
        my $candidates_hash = miRkwood::PrecursorBuilder::merge_candidates( \@sorted_hairpin_candidates_for_chr );

        # Posteriori tests and update candidate information
        my $final_candidates_hash = $precursorBuilderJob->process_mirna_candidates( $candidates_hash );

        undef $candidates_hash;

        $self->serialize_candidates($final_candidates_hash);

        miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : compute_candidates() (sort candidates)' );

        undef $final_candidates_hash;
    }

    return;
}


sub serialize_candidates {
    my ($self, @args) = @_;
    my $candidates = shift @args;
    my @candidates_array = @{$candidates};

    foreach my $candidate (@candidates_array ) {

        #~ $candidate = $candidate->get_reads_from_bed_file($self->{'bed_file'});
        miRkwood::CandidateHandler::print_reads_clouds( $candidate, $self->get_new_reads_dir() );
        miRkwood::CandidateHandler->serialize_candidate_information( $self->get_new_candidates_dir(), $candidate );

        push $self->{'basic_candidates'}, $candidate->get_basic_informations();
    }
    return;
}

=method treat_known_mirnas

 Usage : $self->treat_known_mirnas();
 
=cut

sub treat_known_mirnas {
    my ($self, @args) = @_;

    $self->{'basic_known_candidates'} = [];

    $self->store_known_mirnas_as_candidate_objects();
    $self->serialize_basic_candidates( 'basic_known_candidates' );

    undef $self->{'basic_known_candidates'};

    return;

}

sub store_known_mirnas_as_candidate_objects {
    my ($self, @args) = @_;
    my $job_dir        = $self->{'job_dir'};
    my $bed_file       = $self->{'mirna_bed'};
    my $reads_dir      = miRkwood::Paths::get_known_reads_dir_from_job_dir( $job_dir );
    my $candidates_dir = miRkwood::Paths::get_known_candidates_dir_from_job_dir( $job_dir );
    my $cfg            = miRkwood->CONFIG();
    my $species        = $cfg->param('job.plant');
    my $gff_file       = File::Spec->catfile( miRkwood::Paths->get_data_path(), "miRBase/${species}_miRBase.gff3");

    my @field;
    my ($id, $name);
    my $precursor_reads;
    my $precursor_of_mature;
    my $mature_reads;
    my $data;

    if ( ! -r $gff_file ){
        debug("No miRBase file for $species.", miRkwood->DEBUG() );
        return;
    }

    ##### Read the GFF and links each precursor with its mature

    open (my $GFF, '<', $gff_file) or die "ERROR while opening $gff_file : $!";

    while ( <$GFF> ){
        if ( ! /^#/ ){
            chomp;

            @field = split( /\t/xms );

            if ( $field[2] eq 'miRNA' and $field[8] =~ /ID=([^;]+).*Derives_from=([^;]+)/ ){
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
        my $precursor_id = '';
        my $read_start = $field[1] + 1; # 1-based
        my $read_end   = $field[2];     # 1-based       

        if ($field[14] =~ /ID=([^;]+).*Name=([^;]+)/ ){
            $id = $1;
            $name = $2;
        }

        if ( $field[8] eq 'miRNA_primary_transcript' ){
            $precursor_id = $id;

            $data->{$precursor_id}{'identifier'}      = $id;
            $data->{$precursor_id}{'precursor_name'}  = $name;
            $data->{$precursor_id}{'name'}  = $field[0];
            $data->{$precursor_id}{'length'} = $field[10] - $field[9] + 1;
            $data->{$precursor_id}{'start_position'} = $field[9];
            $data->{$precursor_id}{'end_position'}   = $field[10];
            $data->{$precursor_id}{'position'} = $data->{$precursor_id}{'start_position'} . '-' . $data->{$precursor_id}{'end_position'};
            $data->{$precursor_id}{'reads'}{"$read_start-$read_end"} = $field[4];
        }
        elsif ( $field[8] eq 'miRNA' ){
            $precursor_id = $precursor_of_mature->{$id};
            $data->{$precursor_id}{'matures'}{$id}{'mature_name'}  = $name;
            $data->{$precursor_id}{'matures'}{$id}{'mature_start'} = $field[9];
            $data->{$precursor_id}{'matures'}{$id}{'mature_end'}   = $field[10];
            $data->{$precursor_id}{'matures'}{$id}{'mature_reads'}{"$read_start-$read_end"} = $field[4];
        }

        $data->{$precursor_id}{'chromosome'} = $field[0];
        $data->{$precursor_id}{'strand'}     = $field[5];

    }

    close $BED;

    ##### Treat data by precursor
    foreach my $precursor_id ( keys%{$data} ){

        $precursor_reads = 0;
        $mature_reads = 0;

        ##### Count number of reads
        foreach (keys %{$data->{$precursor_id}{'reads'}}){
            $precursor_reads += $data->{$precursor_id}{'reads'}{$_};
        }
        foreach my $mature_id ( keys %{$data->{$precursor_id}{'matures'}} ){
            foreach my $read ( keys %{$data->{$precursor_id}{'matures'}{$mature_id}{'mature_reads'}} ){
                $mature_reads += $data->{$precursor_id}{'matures'}{$mature_id}{'mature_reads'}{$read};
            }
        }

        ##### Calculate score
        $data->{$precursor_id}{'quality'} = 0;
        if ( $precursor_reads >= 10 ){
           $data->{$precursor_id}{'quality'}++;
        }
        if ( $mature_reads >= ( $precursor_reads / 2 ) ){
            $data->{$precursor_id}{'quality'}++;
        }

        ### Create a Candidate object
        my $candidate = miRkwood::Candidate->new( $data->{$precursor_id} );
        my $candidate_dir = File::Spec->catdir( $self->get_workspace_path, $precursor_id );
        mkdir $candidate_dir;

        my $candidatejob = miRkwood::CandidateJob->new( $candidate_dir,
                                                        $candidate->{'name'},
                                                        $precursor_id,
                                                        $candidate,
                                                        [] );

        $candidate = $candidatejob->update_known_candidate_information( $candidate, $self->{'genome_db'} );

        miRkwood::CandidateHandler->serialize_candidate_information($candidates_dir, $candidate);

        ### Store basic information (used for the HTML table) for this Candidate
        push $self->{'basic_known_candidates'}, $candidate->get_basic_informations();

        ### Create individual card with reads cloud
        miRkwood::CandidateHandler::print_reads_clouds( $candidate, $reads_dir );

    }

    miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : store_known_mirnas_as_candidate_objects' );

    return;

}


1;
