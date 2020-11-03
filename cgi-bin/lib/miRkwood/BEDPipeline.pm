package miRkwood::BEDPipeline;

# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use parent 'miRkwood::Pipeline';

use File::Which;
use Log::Message::Simple qw[msg error debug];
use miRkwood::Paths;
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

    my $cfg = miRkwood->CONFIG();

    debug( 'miRkwood : start processing.' . ' [' . localtime() . ']', miRkwood->DEBUG() );

    my $start_time = time();

    $self->calculate_reads_coverage();

    $self->filter_BED();

    # Look for known miRNAs
    debug( 'Treat known miRNAs.' . ' [' . localtime() . ']', miRkwood->DEBUG() );
    $self->treat_known_mirnas();

    # Compress BED files
    my @list_of_BED_files = qw{_miRNAs};
    my $annotation_gff    = $cfg->param( 'options.annotation_gff' );
    my @annotation_gff    = split( /\&/, $annotation_gff );
    foreach my $gff ( @annotation_gff ){
        if ( $gff =~ /_([^_]+)[.](gtf|gff3?|dat)/ ){
            push @list_of_BED_files, '_'.$1;
        }
    }
    push @list_of_BED_files, '_multimapped';

    my $bed_sizes_file = File::Spec->catfile( $self->get_job_dir(), miRkwood::Paths::get_bed_size_file_name() );
    open (my $FH, '>', $bed_sizes_file) or die "ERROR while creating $bed_sizes_file : $!";
    print $FH "#file\treads\tunique reads\n";
    close $FH;

    miRkwood::BEDHandler::store_reads_nb_in_BED_file( $self->{'initial_bed'}, $bed_sizes_file );
    miRkwood::BEDHandler::store_reads_nb_in_BED_file( $self->{'bed_file'}, $bed_sizes_file );

    foreach my $BED_type ( @list_of_BED_files ){
        my $BED_file = miRkwood::Paths::get_bed_file ( $self->get_job_dir(), $BED_type, 'bed' );
        miRkwood::BEDHandler::store_reads_nb_in_BED_file( $BED_file, $bed_sizes_file );
        miRkwood::BEDHandler::zipBEDfile( $BED_file );
    }

    # Look for new miRNAs
    debug( 'Treat new miRNAs.' . ' [' . localtime() . ']', miRkwood->DEBUG() );

    $self->list_chrom_in_bed();

    debug( scalar(@{$self->{'chromosomes_in_bed'}}) . ' chromosome(s) to consider', miRkwood->DEBUG() );

    $self->{'basic_candidates'} = [];

    my $bed_name = $cfg->param('job.bed');
    $self->{'orphan_clusters'} = File::Spec->catfile( $self->{'job_dir'}, miRkwood::Paths::get_orphan_clusters_file_name( $bed_name ) );
    $self->{'orphan_hairpins'} = File::Spec->catfile( $self->{'job_dir'}, miRkwood::Paths::get_orphan_hairpins_file_name( $bed_name ) );

    foreach my $chromosome ( @{$self->{'chromosomes_in_bed'}} ){
        debug( "- Start chromosome $chromosome" . ' [' . localtime() . ']', miRkwood->DEBUG() );

        debug( "  Start init_sequences_per_chr for $chromosome" . ' [' . localtime() . ']', miRkwood->DEBUG() );
        $self->init_sequences_per_chr( $chromosome );

        debug( "  Start run_pipeline_on_sequences_per_chr for $chromosome" . ' [' . localtime() . ']', miRkwood->DEBUG() );
        $self->run_pipeline_on_sequences_per_chr( $chromosome );

        $self->clean_workspace_per_chr( $chromosome );
    }

    $self->count_alignments_per_miRNA();
    $self->count_precursors_per_miRNA();

	# Count the number of orphan clusters and orphan hairpins and write it in a file
    my $orphan_regions_file = File::Spec->catfile( $self->get_job_dir(), miRkwood::Paths::get_orphan_regions_file_name() );
    open(my $OR, '>', $orphan_regions_file) or die "ERROR while opening $orphan_regions_file: $!";
    my $nb_orphan_clusters = miRkwood::BEDHandler::count_alignments_nb_in_BED_file( $self->{'orphan_clusters'} );
    my $nb_orphan_hairpins = miRkwood::BEDHandler::count_alignments_nb_in_BED_file( $self->{'orphan_hairpins'} );
    print $OR "Orphan clusters:$nb_orphan_clusters\n";
    print $OR "Orphan hairpins:$nb_orphan_hairpins\n";
    close $OR;

    miRkwood::BEDHandler::zipBEDfile( $self->{'orphan_clusters'} );
    miRkwood::BEDHandler::zipBEDfile( $self->{'orphan_hairpins'} );

    $self->serialize_basic_candidates( 'basic_candidates' );

    # Delete workspace
    my $workspace_path = miRkwood::Paths->get_workspace_path( $cfg->param('job.directory') );
    system("rm -Rf $workspace_path");

    $self->mark_job_as_finished();

    debug( 'Writing finish file' . ' [' . localtime() . ']', miRkwood->DEBUG() );

    my $total_time = time() - $start_time;
    my $day  = int( $total_time / 86_400 );
    my $hour = int( ($total_time % 86_400 ) / 3_600 );
    my $min  = int( ( ($total_time % 86_400 ) % 3_600 ) / 60 );
    my $sec  = int( ( ($total_time % 86_400 ) % 3_600 ) % 60 );

    debug( "Done in $day day $hour h $min min $sec sec.", miRkwood->DEBUG() );

    return;
}

=method count_alignments_per_miRNA

  Read the non filtered BED file and the basic_candidates hash
  file to count how many times each read corresponding to a miRNA
  is aligned.
  Modify the basic_candidates hash to add a feature 'nb_alignments_for_miRNA'
  for each candidate.

=cut
sub count_alignments_per_miRNA {
    my ($self, @args) = @_;
    my $nb_alignments_per_read;

    open( my $BED, '<', $self->{'bed_file'}) or die "ERROR while opening $self->{'bed_file'}: $!";
    while ( <$BED> ){
        chomp;
        my @fields = split( /\t/ );
        $fields[3] =~ s/T/U/g;
        if ( ! defined( $nb_alignments_per_read->{ $fields[3] } ) ){
            $nb_alignments_per_read->{ $fields[3] } = 0;
        }
        $nb_alignments_per_read->{ $fields[3] }++;
    }
    close $BED;
    foreach my $candidate ( @{$self->{'basic_candidates'}} ){
        if ( $candidate->{'mirna_sequence'} ne '' ){
            if ( defined( $nb_alignments_per_read->{ $candidate->{'mirna_sequence'} } ) && ($nb_alignments_per_read->{ $candidate->{'mirna_sequence'} } != 0) ){
                $candidate->{'nb_alignments_for_miRNA'} = $nb_alignments_per_read->{ $candidate->{'mirna_sequence'} };
            }
            else{
                $candidate->{'nb_alignments_for_miRNA'} = 1;
            }
            $candidate->{'weight'} = $candidate->{'mirna_depth'} / $candidate->{'nb_alignments_for_miRNA'};
            $candidate->{'weight'} = miRkwood::Utils::restrict_num_decimal_digits($candidate->{'weight'}, 3);
        }
        else {
            $candidate->{'nb_alignments_for_miRNA'} = '';
            $candidate->{'weight'} = '';
        }
    }
    return $self;
}

=method count_precursors_per_miRNA

  

=cut
sub count_precursors_per_miRNA {
    my ($self, @args) = @_;
    my $list_ID_foreach_miRNA;

    foreach my $candidate ( @{$self->{'basic_candidates'}} ){
        if ( $candidate->{'mirna_sequence'} ne '' ){
            push @{$list_ID_foreach_miRNA->{ $candidate->{'mirna_sequence'} } }, $candidate->{'identifier'};
        }
    }
    foreach my $candidate ( @{$self->{'basic_candidates'}} ){
        @{ $candidate->{'list_id_with_same_mirna'} } = ();
        if ( $candidate->{'mirna_sequence'} ne '' ){
            my @list_id = @{$list_ID_foreach_miRNA->{ $candidate->{'mirna_sequence'} } };
            @list_id = miRkwood::Utils::delete_element_in_array( $candidate->{'identifier'}, \@list_id );
            @{ $candidate->{'list_id_with_same_mirna'} } = @list_id;
        }
    }
    return $self;
}


=method create_additional_directories

=cut
sub create_additional_directories {
    my ($self, @args) = @_;
    mkdir miRkwood::Paths::get_new_candidates_dir_from_job_dir( $self->{'job_dir'} );
    mkdir miRkwood::Paths::get_known_candidates_dir_from_job_dir( $self->{'job_dir'} );
    mkdir miRkwood::Paths::get_dir_reads_path_from_job_dir( $self->{'job_dir'} );
    mkdir miRkwood::Paths::get_new_reads_dir_from_job_dir( $self->{'job_dir'} );
    mkdir miRkwood::Paths::get_known_reads_dir_from_job_dir( $self->{'job_dir'} );
    return;
}


=method calculate_reads_coverage

  On the BED given by the user, calculate the read coverage, ie the 
  number X such as on average one read starts every X nt.

=cut
sub calculate_reads_coverage {
    my ($self, @args) = @_;
    my $line;
    my $data = {};
    open (my $BED, $self->{'initial_bed'}) or die "ERROR while opening $self->{'initial_bed'} : $!";
    while( <$BED> ){
        chomp;
        if ( $_ ne '' ){
            my @fields = split( /\t/);
            if ( ! defined($data->{ $fields[0] }{'start_position'}) || $data->{ $fields[0] }{'start_position'} eq '' ){
                $data->{ $fields[0] }{'start_position'} = $fields[1];
            }
            $data->{ $fields[0] }{'end_position'} = $fields[2];
            $data->{ $fields[0] }{'nb_tot_reads'} += $fields[4];
            $line = $_;
        }
    }
    close $BED;
    chomp $line;
    my @fields = split( /\t/, $line);
    $data->{ $fields[0] }{'end_position'} = $fields[2];

    foreach my $chromosome ( sort (keys%{ $data })){
        $data->{$chromosome}{'size_genome'} = $data->{$chromosome}{'end_position'} - $data->{$chromosome}{'start_position'};
        $self->{'average_coverage'}{$chromosome} = int( $data->{$chromosome}{'size_genome'} / $data->{$chromosome}{'nb_tot_reads'} );
        if ( $self->{'average_coverage'}{$chromosome} == 0 ){
            $self->{'average_coverage'}{$chromosome} = 1;
        }
        debug( "Average reads coverage for chromosome $chromosome : 1 read every $self->{'average_coverage'}{$chromosome} nt.", miRkwood->DEBUG() );
    }

    return $self;
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
    my $localBED           = $self->{'initial_bed'};
    my $filteredBED        = '';
    my $mirnaBED           = '';

    if ( ! which( 'intersectBed' ) ){
        debug('[WARNING] bedtools are not installed. Cannot filter the BED file.', miRkwood->DEBUG());
    }
    else{
        ($filteredBED, $mirnaBED) = miRkwood::BEDHandler->filterBEDfile( $localBED );
    }

    if ( ! defined($filteredBED) || $filteredBED eq '' ){
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

    my @list_chrom_in_genome = $self->{'genome_db'}->get_all_primary_ids;
    my %hash_chrom_in_genome;
    @hash_chrom_in_genome{@list_chrom_in_genome} = 0..$#list_chrom_in_genome;

    my @final_list;
    foreach my $key ( @sorted_list ){
        if ( exists $hash_chrom_in_genome{$key} ){
            push @final_list, $key;
        }
    }

    $self->{'chromosomes_in_bed'} = \@final_list;

    return;
}


=method init_sequences_per_chr

=cut
sub init_sequences_per_chr {
    my ($self, @args) = @_;
    my $chromosome = shift @args;
    my $clustering = miRkwood::ClusterBuilder->new($self->{'genome_db'}, $self->{'bed_file'});
    $self->{'sequences'} = $clustering->build_loci_per_chr( $chromosome, $self->{'average_coverage'} );
    $self->{'parsed_reads'} = $clustering->get_parsed_bed();
    miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : init_sequences_per_chr()');
    return;
}


sub run_pipeline_on_sequences_per_chr {
    my ($self, @args) = @_;
    my $chromosome = shift @args;

    my $sequences_count = scalar( @{$self->{'sequences'}} );

    debug( "   $sequences_count sequences to process", miRkwood->DEBUG() );

    $self->compute_candidates_per_chr( $chromosome );

    return;
}

sub compute_candidates_per_chr {
    my ($self, @args) = @_;
    my $chromosome = shift @args;
    my $distance_min = 2_000;

    my $sequence_identifier = 0;

    my $hairpinBuilder = miRkwood::HairpinBuilder->new($self->{'genome_db'}, $self->get_workspace_path(), $self->{'parsed_reads'});
    my @hairpin_candidates_for_chr = ();
    my $previous_end = 0;

    foreach my $locus ( @{$self->{'sequences'}} ) {

        debug( "    Start sequence $sequence_identifier" . ' [' . localtime() . ']', miRkwood->DEBUG() );
        $sequence_identifier++;

        if ( ( $locus->{'begin'} - $previous_end ) < $distance_min ){
            push @hairpin_candidates_for_chr, @{ $hairpinBuilder->build_hairpins($locus) };
        }
        else {
            # Treat previous bunch
            debug( '    * Treat previous bunch', miRkwood->DEBUG());
            undef $hairpinBuilder;
            my @sorted_hairpin_candidates_for_chr = sort { $a->{'start_position'} <=> $b->{'start_position'} } @hairpin_candidates_for_chr;
            undef @hairpin_candidates_for_chr;

            if ( scalar(@sorted_hairpin_candidates_for_chr) ){
                my $precursorBuilderJob = miRkwood::PrecursorBuilder->new( $self->get_workspace_path(), $self->{'genome_db'}, $chromosome, $chromosome );

                # Merge candidates
                debug( "     - Merge candidates (current locus : $locus->{'begin'}-$locus->{'end'})" . ' [' . localtime() . ']', miRkwood->DEBUG() );
                my $candidates_hash = miRkwood::PrecursorBuilder::merge_candidates( \@sorted_hairpin_candidates_for_chr );
                undef @sorted_hairpin_candidates_for_chr;

                # Posteriori tests and update candidate information
                debug("     - Process miRNA candidates (current locus : $locus->{'begin'}-$locus->{'end'})" . ' [' . localtime() . ']', miRkwood->DEBUG() );
                my $final_candidates_hash = $precursorBuilderJob->process_mirna_candidates( $candidates_hash );
                undef $candidates_hash;

                debug("     - serialize_candidates (current locus : $locus->{'begin'}-$locus->{'end'})" . ' [' . localtime() . ']', miRkwood->DEBUG() );
                $self->serialize_candidates($final_candidates_hash);

            }

            # Re-initialize some variables
            @hairpin_candidates_for_chr = ();
            $hairpinBuilder = miRkwood::HairpinBuilder->new($self->{'genome_db'}, $self->get_workspace_path(), $self->{'parsed_reads'});

            # Treat current locus
            debug( '    * Treat current locus', miRkwood->DEBUG());
            push @hairpin_candidates_for_chr, @{ $hairpinBuilder->build_hairpins($locus) };
        }

        $previous_end = $locus->{'end'};
        miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : compute_candidates_per_chr_new() (boucle sur loci)' );

    }   # end foreach locus

    # Treat last bunch
    my @sorted_hairpin_candidates_for_chr = sort { $a->{'start_position'} <=> $b->{'start_position'} } @hairpin_candidates_for_chr;
    undef @hairpin_candidates_for_chr;

    if ( scalar(@sorted_hairpin_candidates_for_chr) ){
        my $precursorBuilderJob = miRkwood::PrecursorBuilder->new( $self->get_workspace_path(), $self->{'genome_db'}, $chromosome, $chromosome );

        # Merge candidates
        debug( '     - Merge candidates (current locus : last one)' . ' [' . localtime() . ']', miRkwood->DEBUG() );
        my $candidates_hash = miRkwood::PrecursorBuilder::merge_candidates( \@sorted_hairpin_candidates_for_chr );
        undef @sorted_hairpin_candidates_for_chr;

        # Posteriori tests and update candidate information
        debug('     - Process miRNA candidates (current locus : last one)' . ' [' . localtime() . ']', miRkwood->DEBUG() );
        my $final_candidates_hash = $precursorBuilderJob->process_mirna_candidates( $candidates_hash );
        undef $candidates_hash;

        $self->serialize_candidates($final_candidates_hash);
        miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : compute_candidates_per_chr_new() (boucle sur loci)' );

    }

    return;

}


sub serialize_candidates {
    my ($self, @args) = @_;
    my $candidates = shift @args;
    my @candidates_array = @{$candidates};

    foreach my $candidate (@candidates_array ) {
        miRkwood::CandidateHandler::print_reads_clouds( $candidate, miRkwood::Paths::get_new_reads_dir_from_job_dir( $self->{'job_dir'} ) );
        $candidate->{'mfei'} = miRkwood::Utils::restrict_num_decimal_digits( $candidate->{'mfei'}, 3);
        $candidate->{'amfe'} = miRkwood::Utils::restrict_num_decimal_digits( $candidate->{'amfe'}, 3);
        miRkwood::CandidateHandler->serialize_candidate_information( miRkwood::Paths::get_new_candidates_dir_from_job_dir( $self->{'job_dir'} ), $candidate );
        push @{ $self->{'basic_candidates'} }, $candidate->get_basic_informations();
    }
    return;
}

=method treat_known_mirnas

 Usage : $self->treat_known_mirnas();
 
=cut

sub treat_known_mirnas {
    my ($self, @args) = @_;

    $self->{'basic_known_candidates'} = [];

    if ( $self->{'mirna_bed'} ne '' ){
        $self->store_known_mirnas_as_candidate_objects();
    }
    else{
        debug( '   No BED for known miRNAs.', miRkwood->DEBUG() );
    }
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
    my $gff_file       = $cfg->param('options.mirbase_gff');

    my @field;
    my ($id, $name);
    my $mature_informations_for_precursor;
    my $data;

    if ( ! -r $self->{'mirna_bed'} || ( ! -s $self->{'mirna_bed'} ) ){
        debug( "No reads corresponding to known miRNAs for $species", miRkwood->DEBUG() );
        return;
    }
    if ( ! -r $gff_file ){
        debug("No miRBase file for $species.", miRkwood->DEBUG() );
        return;
    }

    ##### Read the GFF and links each precursor with its mature
    # /!\ If a precursor does not have any mature, it will not appear in this hash
    open (my $GFF, '<', $gff_file) or die "ERROR while opening $gff_file : $!";

    while ( <$GFF> ){
        if ( ! /^#/ ){
            chomp;

            @field = split( /\t/xms );

            if ( $field[2] eq 'miRNA' and $field[8] =~ /ID=([^;]+).*Name=([^;]+).*Derives_from=([^;]+)/ ){
                $mature_informations_for_precursor->{$3}{'mature_id'} = $1;
                $mature_informations_for_precursor->{$3}{'mature_name'} = $2;
                $mature_informations_for_precursor->{$3}{'mature_start'} = $field[3];
                $mature_informations_for_precursor->{$3}{'mature_end'} = $field[4];
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

            $data->{$precursor_id}{'identifier'} = $id;
            $data->{$precursor_id}{'precursor_name'} = $name;
            $data->{$precursor_id}{'mirbase_id'} = $id;
            $data->{$precursor_id}{'name'}  = $field[0];
            $data->{$precursor_id}{'length'} = $field[10] - $field[9] + 1;
            $data->{$precursor_id}{'start_position'} = $field[9];
            $data->{$precursor_id}{'end_position'}   = $field[10];
            $data->{$precursor_id}{'position'} = $data->{$precursor_id}{'start_position'} . '-' . $data->{$precursor_id}{'end_position'};
            $data->{$precursor_id}{'reads'}{"$read_start-$read_end"} = $field[4];

            my $mature_id = $mature_informations_for_precursor->{$precursor_id}{'mature_id'};
            if ( defined($mature_id) and $mature_id ne '' ){
                $data->{$precursor_id}{'matures'}{$mature_id}{'mature_name'}  = $mature_informations_for_precursor->{$precursor_id}{'mature_name'};
                $data->{$precursor_id}{'matures'}{$mature_id}{'mature_start'} = $mature_informations_for_precursor->{$precursor_id}{'mature_start'};
                $data->{$precursor_id}{'matures'}{$mature_id}{'mature_end'}   = $mature_informations_for_precursor->{$precursor_id}{'mature_end'};
            }
        }
        elsif ( $field[8] eq 'miRNA' ){
            if ( $field[14] =~ /Derives_from=([^;]+)/ ){
                $precursor_id = $1;
                $data->{$precursor_id}{'matures'}{$id}{'mature_name'}  = $name;
                $data->{$precursor_id}{'matures'}{$id}{'mature_start'} = $field[9];
                $data->{$precursor_id}{'matures'}{$id}{'mature_end'}   = $field[10];
                $data->{$precursor_id}{'matures'}{$id}{'mature_reads'}{"$read_start-$read_end"} = $field[4];
            }

        }

        $data->{$precursor_id}{'chromosome'} = $field[0];
        $data->{$precursor_id}{'strand'}     = $field[5];

    }

    close $BED;

    ##### Treat data by precursor
    foreach my $precursor_id ( keys%{$data} ){

        ### Create a Candidate object
        my $candidate = miRkwood::Candidate->new( $data->{$precursor_id} );
        my $candidate_dir = miRkwood::Paths::create_folder( File::Spec->catdir( $self->get_workspace_path, $precursor_id ) );

        my $candidatejob = miRkwood::CandidateJob->new( $candidate_dir,
                                                        $candidate->{'name'},
                                                        $precursor_id,
                                                        $candidate,
                                                        [] );

        $candidate = $candidatejob->update_known_candidate_information( $candidate, $self->{'genome_db'} );

        miRkwood::CandidateHandler->serialize_candidate_information($candidates_dir, $candidate);

        ### Store basic information (used for the HTML table) for this Candidate
        push @{ $self->{'basic_known_candidates'} }, $candidate->get_basic_informations();

        ### Create individual card with reads cloud
        miRkwood::CandidateHandler::print_reads_clouds( $candidate, $reads_dir );

    }

    miRkwood::Utils::display_var_sizes_in_log_file( '..... BEDPipeline : store_known_mirnas_as_candidate_objects' );

    return;

}

sub clean_workspace_per_chr {
    my ($self, @args) = @_;
    my $chromosome = shift @args;
    my %orphan_clusters;

    my $workspace_chr_path = miRkwood::Paths::get_workspace_chromosome_dir( $self->get_workspace_path(), $chromosome );
    my @list_clusters = glob "$workspace_chr_path/*";

    foreach my $cluster ( @list_clusters ){
        my $candidate_directory_found = 0;
        my @content = glob "$cluster/*";
        foreach my $file ( @content ){
            if ( -d $file ){
                $candidate_directory_found = 1;
            }
        }
        if ( ! $candidate_directory_found ){
            my $cluster_start = 0;
            my $cluster_end = 0;
            my $strand = '';
            if ( $cluster =~ /.*\/(\d+)-(\d+)([+-])/ ){
                $cluster_start = $1;
                $cluster_end = $2;
                $strand = $3;
            }

            my $reads = miRkwood::HairpinBuilder::get_contained_reads( $self->{'parsed_reads'}, $chromosome, $cluster_start, $cluster_end, $strand);
            my $nb_reads = 0;
            foreach my $key ( keys%{$reads} ){
                $nb_reads += $reads->{$key};
            }
            $orphan_clusters{ $cluster_start } = $cluster_end."\t".'orphan_cluster'."\t".$nb_reads."\t".$strand;

            system( "rm -Rf $cluster" );
        }
    }

    open ( my $FILE, '>>', $self->{'orphan_clusters'} ) or die "ERROR while opening $self->{'orphan_clusters'} : $!";
    foreach ( sort {$a <=> $b} keys%orphan_clusters ){
        print $FILE "$chromosome\t$_\t$orphan_clusters{$_}\n";
    }
    close $FILE;

    return;
}


1;
