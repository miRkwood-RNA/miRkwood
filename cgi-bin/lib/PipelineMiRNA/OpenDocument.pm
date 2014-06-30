package PipelineMiRNA::OpenDocument;

# ABSTRACT: Exporting pipeline results as OpenDocument

use strict;
use warnings;

use File::Spec;
use File::Copy;
use PipelineMiRNA::Utils;
use PipelineMiRNA::WebPaths;
use PipelineMiRNA::Results;
use PipelineMiRNA::Candidate;
use ODF::lpOD;

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my $jobId = shift @args;
    my $sequences_to_export = shift @args;
    my $self = bless {
        jobId => $jobId,
        sequences => $sequences_to_export
    }, $class;
    $self->{doc} = odf_document->create('text')
      or die 'Error when initialising ODF document';
    $self->prepare_document();
    return $self;
}

=method prepare_document

Return a prepared ODF document object
with styles already set.

=cut

sub prepare_document {
    my ( $self, @args ) = @_;
    my $elt;

    # Default paragraph style creation
    $elt = $self->{doc}->insert_style(
        odf_create_style(
            'paragraph',
            margin_top    => '2mm',
            margin_bottom => '2mm',
            orphans       => 2,
            widows        => 2
        ),
        default => TRUE
    );
    $elt->set_properties(
        area     => 'text',
        language => 'none',
        country  => 'none'
    );

    # Basic paragraph style creation
    odf_style->create(
        'paragraph',
        name        => 'Basic',
        margin_top  => '0mm',
        margin_left => '0mm',
        area   => 'text',
        size   => '10pt',
    )->register($self->{doc});

    # None paragraph style creation
    odf_style->create(
        'paragraph',
        name   => 'None',
        parent => 'Basic',
        area   => 'text',
        size   => '0pt',
    )->register($self->{doc});

    # Other MiRBase sequences
    my $s0 = odf_style->create(
        'paragraph',
        name => 'MirbaseSequences',
        parent => 'Basic');
    $s0->set_properties(
        area=>'paragraph',
        margin_left=>'4mm',
        margin_top=>'2mm',
        margin_bottom=>'4mm',
    );
    $s0->set_properties(
        area   => 'text',
        size   => '9pt',
    );
    $self->{doc}->register_style($s0);

    # Monospace
    $self->{doc}->set_font_declaration('Monospace');
    my $s = odf_style->create('paragraph', name => 'Monospace');
    $s->set_properties(
        area=>'paragraph',
        margin_left=>'4mm',
        margin_right=>'4mm',
        margin_top=>'4mm',
    );
    $s->set_properties(
        area   => 'text',
        size   => '10pt',
        color  => 'black',
        font   => 'Monospace',
        );
    $self->{doc}->register_style($s);

    # Vienna (Monospace style)
    odf_style->create(
        'paragraph',
        name  => 'Vienna',
        parent => 'Monospace',
        area   => 'text',
        size   => '8pt',
    )->register($self->{doc});

    # Hairpin (Monospace style)
    odf_style->create(
        'paragraph',
        name  => 'Hairpin',
        parent => 'Monospace',
        area   => 'text',
        size   => '8pt',
    )->register($self->{doc});

    # Alignment (Monospace style)
    odf_style->create(
        'paragraph',
        name  => 'Alignment',
        parent => 'Monospace',
        area   => 'text',
        size   => '9pt',
    )->register($self->{doc});

    # StandardBold (Monospace style)
    odf_style->create(
        'paragraph',
        name  => 'StandardBold',
        parent => 'Standard',
        area   => 'text',
        weight => 'bold',
    )->register($self->{doc});

   # Level 1 Heading style creation
    odf_style->create(
        'paragraph',
        name           => 'Level 1 Heading',
        keep_with_next => 'always',
        margin_top     => '1cm',
        margin_bottom  => '4mm',
        area   => 'text',
        size   => '14pt',
        weight => 'bold',
    )->register($self->{doc});

    # Level 2 Heading style creation
    odf_style->create(
        'paragraph',
        name   => 'Level 2 Heading',
        parent => 'Level 1 Heading',
        area   => 'text',
        size   => '12pt',
        weight => 'bold',
        style  => 'italic',
    )->register($self->{doc});

    # Level 3 Heading style creation
    odf_style->create(
        'paragraph',
        name           => 'Level 3 Heading',
        keep_with_next => 'always',
        margin_top     => '2mm',
        margin_bottom  => '1mm',
        area   => 'text',
        size   => '11pt',
        weight => 'bold',
    )->register($self->{doc});

    # Level 4 Heading style creation
    odf_style->create(
        'paragraph',
        name  => 'Level 4 Heading',
        parent => 'Level 3 Heading',
        area   => 'text',
        size   => '10pt',
    )->register($self->{doc});

    # Level 5 Heading style creation
    odf_style->create(
        'paragraph',
        name  => 'Level 5 Heading',
        parent => 'Level 4 Heading',
        area   => 'text',
    )->register($self->{doc});

    # top title style
    odf_style->create(
        'paragraph',
        name          => 'Top title',
        align         => 'center',
        margin_top    => '0cm',
        margin_bottom => '1cm',
        area   => 'text',
        size   => '20pt',
        weight => 'bold',
        color  => 'navy blue'
      )->register($self->{doc});

    # Graphic style
    odf_style->create(
            'graphic',
            name        => 'Classic',
            align       => 'center',
            margin_top  => '5mm'
    )->register($self->{doc});

    odf_style->create(
        'paragraph',
        name        => 'Centré',
        align       => 'center',
        margin_top  => '5mm'
    )->register($self->{doc});
    return;
}

=method generate_report

Generate an ODF document for the given jobID

Write it on disk and return the server path to it.

=cut

sub generate_report {
    my ( $self, @args ) = @_;
    my $jobId = $self->{jobId};
    my @sequences_to_export = @{$self->{sequences}};

    my $jobPath = PipelineMiRNA::Results->jobId_to_jobPath($jobId);

    my $images_dir = $self->get_images_dir();
    mkdir $images_dir;

    my %results = PipelineMiRNA::Results->get_structure_for_jobID($jobId);

    my $doc = $self->{doc};

    # Main context access
    my $context = $doc->body;

    # Metadata access
    my $meta = $doc->meta;

    $meta->set_generator('PipelineMiRNA');
    $meta->set_title('PipelineMiRNA results');

    # make sure that the document body is empty
    $context->clear;

    # put the main title
    $context->append_element(
        odf_create_heading(
            level => 1,
            text  => 'miRkwood results',
            style => 'Top title'
        )
    );

    my @keys = sort { ( $results{$a}{'name'} cmp
                        $results{$b}{'name'} )
                      ||
                      ( $results{$a}{'start_position'} <=>
                        $results{$b}{'start_position'} )
                 } keys %results;

    foreach my $key (@keys) {
        if ( $key ~~ \@sequences_to_export )
        {
            my $candidate = $results{$key};
            $self->add_candidate($context, $candidate, $key);
        }   # if key in tab
    }    #  while each %results

    # save the generated document and quit
    $doc->save( target => $self->get_ODF_absolute_path(), pretty => TRUE );
    return;
}

=method add_candidate

=cut

sub add_candidate {
    my ( $self, @args ) = @_;
    my $context   = shift @args;
    my $candidate = shift @args;
    my $key = shift @args;
    my $images_dir = $self->get_images_dir();
    my ( $start, $end ) = (${$candidate}{'start_position'}, ${$candidate}{'end_position'} );
    $context->append_element(
        odf_create_heading(
            level => 1,
            style => 'Level 1 Heading',
            text  => "Sequence: ${$candidate}{'name'}, ${$candidate}{'position'}",
        )
    );

    my $para =
      $context->append_element( odf_create_paragraph() );

    my $list = $para->insert_element(odf_create_list, position => NEXT_SIBLING);

    my $size = length ${$candidate}{'DNASequence'};
    my $vienna_seq =
      PipelineMiRNA::Candidate->make_Vienna_viz( ${$candidate}{'Vienna'},
        ${$candidate}{'DNASequence'} );

    $list->add_item(text => "Name: ${$candidate}{'name'}", style => 'Basic');
    $list->add_item(text => "Position: ${$candidate}{'position'} ($size nt)", style => 'Basic');
    $list->add_item(text => "Strand: ${$candidate}{'strand'}", style => 'Basic');
    $list->add_item(text => "G+C content: ${$candidate}{'%GC'}%", style => 'Basic');

    my $subtext = qw{};
    if(${$candidate}{'Vienna'} ne ${$candidate}{'Vienna_optimal'}){
        $subtext .= qw{}
    } else {
        $subtext.= 'This stem-loop structure is the MFE structure'
    }
    my $para1 = $context->append_element(
    odf_create_paragraph(
        text    => $vienna_seq,
        style   =>'Vienna'
        )
    );
    my $para2 = $context->append_element(
    odf_create_paragraph(
        text    => $subtext,
        style   =>'Basic'
        )
    );
    $para2->set_span(filter  => 'structure', style   => 'StandardBold');


    # Copying the image
    my $img_path      = ${$candidate}{'image'};
    my $img_full_path = $img_path;
    my $new_img_path = File::Spec->catfile($images_dir, "$key.png");
    copy($img_full_path, $new_img_path)
        or die "Copy of $img_full_path to $images_dir failed: $!";

    my ( $lien_image, $taille_image ) =
      $self->{doc}->add_image_file($new_img_path);

    my $factor = 0.5;
    my @width = split( 'pt', shift $taille_image);
    my $width = ($width[0] * $factor) . 'pt';
    my @height = split( 'pt', shift $taille_image);
    my $height = ($height[0] * $factor) . 'pt';
    my $new_size = [$width, $height];

    $para->append_element(
        odf_frame->create(
            image => $lien_image,
            name  => "Structure_${$candidate}{'name'}_${$candidate}{'position'}",
            title => 'Structure',
            description => 'Structure',
            size => $new_size,
        )
    );

   # Section Thermodynamics structure
    $self->add_thermodynamics_section($context, $candidate);

    # Section Mirbase alignments
    $self->add_ODF_alignments($context, $candidate);
    return;
}


=method add_thermodynamics_section

=cut

sub add_thermodynamics_section {
    my ( $self, @args ) = @_;
    my $context   = shift @args;
    my $candidate = shift @args;
    $context->append_element(
        odf_create_heading(
            level => 3,
            style => 'Level 3 Heading',
            text  => 'Thermodynamics stability',
        )
    );

    my $para = $context->append_element( odf_create_paragraph(style =>'None') );
    my $list = $para->insert_element(odf_create_list, position => NEXT_SIBLING);

    # TODO: maybe we do not have those ; infer that from  run_options config file
    $list->add_item(text => "MFE: ${$candidate}{'mfe'} kcal/mol", style => 'Basic');
    $list->add_item(text => "AMFE: ${$candidate}{'amfe'}", style => 'Basic');
    $list->add_item(text => "MFEI: ${$candidate}{'mfei'}", style => 'Basic');
    return;
}

=method add_ODF_alignments

=cut

sub add_ODF_alignments {
    my ( $self, @args ) = @_;
    my $context   = shift @args;
    my %candidate = %{ shift @args };
    $context->append_element(
        odf_create_heading(
            level => 3,
            style => 'Level 3 Heading',
            text  => 'Conserved mature miRNA',
        )
    );

    if ( !$candidate{'alignment_existence'} ) {
        $context->append_element(
            odf_create_paragraph(
                text  => 'No alignment has been found.',
                style => 'Basic'
            )
        );
        return 0;
    }
    my %alignments     = %{ $candidate{'alignments'} };
    my %mirdup_results = %{ $candidate{'mirdup_validation'} };

    my $predictionCounter = 0;

    # Sorting by position
    my @keys = sort { ( PipelineMiRNA::Utils::get_element_of_split($a, '-', 0)  <=>
                        PipelineMiRNA::Utils::get_element_of_split($b, '-', 0)
                      ) ||
                      ( PipelineMiRNA::Utils::get_element_of_split($a, '-', 1)  <=>
                        PipelineMiRNA::Utils::get_element_of_split($b, '-', 1))
                    } keys %alignments;
    foreach my $position (@keys) {
        my ( $left, $right ) = split( /-/, $position );

        # MiRdup
        my $mirdup_key = $candidate{'name'} . '__' . $position;
        my $mirdup_prediction;
        if ( $mirdup_results{$mirdup_key} ) {
            $mirdup_prediction = 'This prediction is validated by miRdup.';
        }
        else {
            $mirdup_prediction = 'This prediction is not validated by miRdup.';
        }

        # Hairpin
        my $hairpin_with_mature =
            PipelineMiRNA::Utils::make_hairpin_with_mature($candidate{'hairpin'},
                                                           $left, $right,
                                                           length $candidate{'DNASequence'},
                                                           'ascii');

        $predictionCounter += 1;

        # Sorting the hit list by descending value of the 'score' element
        my @hits =
          reverse sort { $a->{'score'} <=> $b->{'score'} } @{ $alignments{$position} };
        my $title = "Prediction $predictionCounter: $position";

        $context->append_element(
            odf_create_heading(
                level => 4,
                style => 'Level 4 Heading',
                text  => $title,
            )
        );
        $context->append_element(
            odf_create_paragraph(
                text  => $hairpin_with_mature,
                style => 'Hairpin'
            )
        );

        $context->append_element(
            odf_create_paragraph(
                text  => $mirdup_prediction,
                style => 'Basic'
            )
        );

        $context->append_element(
            odf_create_heading(
                level => 5,
                style => 'Level 5 Heading',
                text  => 'Alignments',
            )
        );

        foreach my $hit (@hits){
            my $alignment = $hit->{'alignment'};
            my $names = $hit->{'name'} . q{ } . $hit->{'def_query'};

            my $name;
            my @splitted = split('\|', $names);

            my $spacing = 15;
            my ($top, $middle, $bottom) = split(/\n/, $alignment);
            $top    = sprintf "%-${spacing}s %3s %s %s", 'query', $hit->{'begin_target'}, $top,   $hit->{'end_target'};
            $middle = sprintf "%-${spacing}s %3s %s %s", qw{},    qw{},                   $middle, qw{};

            my $mirbase_title = qw{};
            if( (scalar @splitted) > 1 ) {
                $mirbase_title = 'miRBase sequences: ';
            }else{
                $mirbase_title = 'miRBase sequence: ';
            }

            $name = 'miRBase';

            my @mirbase_links;
            my @mirbase_ids;
            foreach my $seq (@splitted){
                $seq =~ s/^\s+//;
                $seq =~ s/\s+$//;
                if ($seq =~ 'revcomp'){
                } else {
                    my @splitted_one = split(/ /, $seq);
                    my $local_name = $splitted_one[0];
                    my $mirbase_id = $splitted_one[1];
                    my $mirbase_link = PipelineMiRNA::Utils::make_mirbase_link($mirbase_id);
                    push @mirbase_links, $mirbase_link;
                    push @mirbase_ids, $local_name;
                }
            }
            my $para_seqs = odf_create_paragraph(
                text    => $mirbase_title . join(', ', @mirbase_ids),
                style   => 'MirbaseSequences'
            );
            foreach my $i (0..$#mirbase_ids){
                my $mirbase_id = $mirbase_ids[$i];
                my $mirbase_link = $mirbase_links[$i];
                $para_seqs->set_hyperlink(
                    filter  => $mirbase_id,
                    url     => $mirbase_link,
                    name    => "MiRBase entry for $mirbase_id"
                );
            }
            $bottom = sprintf "%-${spacing}s %3s %s %s", $name,   $hit->{'begin_query'},  $bottom, $hit->{'end_query'};

            $context->append_element(
            odf_create_paragraph(
                text    => "$top\n$middle\n$bottom",
                style   =>'Alignment'
                )
            );
            $context->append_element($para_seqs);
        } # foreach @hits
    } # foreach @keys
    return;
}


=method get_ODF_server_path

Return the server path to the ODF document

=cut

sub get_ODF_server_path{
    my ( $self, @args ) = @_;
    my $ODT_abspath = $self->get_ODF_absolute_path();
    return PipelineMiRNA::WebPaths->filesystem_to_relative_path($ODT_abspath);
}

=method get_ODF_absolute_path

Return the absolute path to the ODF document

=cut

sub get_ODF_absolute_path {
    my ( $self, @args ) = @_;
    my $jobPath = PipelineMiRNA::Results->jobId_to_jobPath($self->{jobId});
    my $ODT_filename = $self->get_ODF_filename();
    return File::Spec->catfile( $jobPath, $ODT_filename );
}

=method get_ODF_filename

Return the document filename

=cut

sub get_ODF_filename {
    my ( $self, @args ) = @_;
    return "Prediction_report_$self->{jobId}.odt";
}

sub get_images_dir {
    my ( $self, @args ) = @_;
    my $jobPath = PipelineMiRNA::Results->jobId_to_jobPath($self->{jobId});
    return File::Spec->catdir($jobPath, 'images');
}

=method get_report

Gererates the report
and return the server path to it.

=cut

sub get_report {
    my ( $self, @args ) = @_;
    $self->generate_report();
    return $self->get_ODF_server_path();;
}

1;
