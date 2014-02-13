package PipelineMiRNA::OpenDocument;

# ABSTRACT: Exporting pipeline results as OpenDocument

use strict;
use warnings;

use File::Spec;
use File::Copy;
use PipelineMiRNA::Utils;
use PipelineMiRNA::Paths;
use PipelineMiRNA::Results;
use PipelineMiRNA::Candidate;
use ODF::lpOD;

=method prepare_document

Return a prepared ODF document object
with styles already set.

=cut

sub prepare_document {
    my ( $self, @args ) = @_;
    my $doc = odf_document->create('text')
      or die 'Error when initialising ODF document';

    my $elt;

    # Default paragraph style creation
    $elt = $doc->insert_style(
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
    $elt = $doc->insert_style(
        odf_create_style(
            'paragraph',
            name        => 'Basic',
            margin_top  => '0mm',
            margin_left => '0mm'
        )
    );

    # Other MiRBase sequences
    my $s0 = odf_style->create(
        'paragraph',
        name => "MirbaseSequences",
        parent => "Basic");
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
    $doc->register_style($s0);

    # Monospace
    $doc->set_font_declaration("Monospace");
    my $s = odf_style->create('paragraph', name => "Monospace");
    $s->set_properties(
        area=>'paragraph',
        margin_left=>'4mm',
        margin_right=>'4mm',
        margin_top=>'4mm',
    );
    $s->set_properties(
            area   => 'text',
            size   => '11pt',
            color  => 'black',
            font   => 'Monospace',
        );
    $doc->register_style($s);

    # Hairpin (Monospace style)
    odf_style->create(
        'paragraph',
        name  => "Hairpin",
        parent => "Monospace",
        area   => 'text',
        size   => '8pt',
    )->register($doc);

    # Alignment (Monospace style)
    odf_style->create(
        'paragraph',
        name  => "Alignment",
        parent => "Monospace",
        area   => 'text',
        size   => '9pt',
    )->register($doc);

    # StandardBold (Monospace style)
    odf_style->create(
        'paragraph',
        name  => "StandardBold",
        parent => "Standard",
        area   => 'text',
        weight => 'bold',
    )->register($doc);

    # Level 2 Heading style creation
    $doc->insert_style(
        odf_create_style(
            'paragraph',
            name           => 'Level 2 Heading',
            keep_with_next => 'always',
            margin_top     => '1cm',
            margin_bottom  => '4mm'
        )
    )->set_properties(
        area   => 'text',
        size   => '16pt',
        weight => 'bold',
        style  => 'italic',
        color  => 'navy blue'
    );

    # Level 3 Heading style creation
    $doc->insert_style(
        odf_create_style(
            'paragraph',
            name           => 'Level 3 Heading',
            keep_with_next => 'always',
            margin_top     => '1cm',
            margin_bottom  => '4mm'
        )
    )->set_properties(
        area   => 'text',
        size   => '14pt',
        weight => 'bold',,
        color  => 'navy blue'
    );

    # top title style
    $doc->insert_style(
        odf_create_style(
            'paragraph',
            name          => 'Top title',
            align         => 'center',
            margin_top    => '0cm',
            margin_bottom => '1cm'
        )
      )->set_properties(
        area   => 'text',
        size   => '200%',
        weight => 'bold',
        color  => 'navy blue'
      );

    # Graphic style
    $doc->insert_style(
        odf_create_style(
                'graphic',
                name        => "Classic",
                align       => 'center',
                margin_top  => '5mm'
                ),
        automatic       => TRUE
        );

    $doc->insert_style(
        odf_create_style(
            'paragraph',
            name        => "Centré",
            align       => 'center',
            margin_top  => '5mm'
        )
    );

    return $doc;
}

=method generate_report

Generate an ODF document for the given jobID

Write it on disk and return the server path to it.

=cut

sub generate_report {
    my ( $self, @args ) = @_;
    my $jobId = shift @args;
    my @sequences_to_export = shift @args;

    my $jobPath = PipelineMiRNA::Results->jobId_to_jobPath($jobId);

    my $images_dir = File::Spec->catdir($jobPath, 'images');
    mkdir $images_dir;

    my ($ODT_abspath, $ODT_serverpath) = $self->get_ODF_path($jobId);

    my %results = PipelineMiRNA::Results->get_structure_for_jobID($jobId);
    use Data::Dumper;

    my $doc = $self->prepare_document();

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
            text  => 'PipelineMiRNA results',
            style => 'Top title'
        )
    );

    my @keys    = sort keys %results;

    foreach my $key (@keys) {
        if ( $key ~~ \@sequences_to_export )
        {
            my $candidate = $results{$key};
            my ( $start, $end ) = split( m/[-]/xms, ${$candidate}{'position'} );
            my $candidate_name = PipelineMiRNA::Candidate->get_name($candidate);
            $context->append_element(
                odf_create_heading(
                    level => 1,
                    text  => $candidate_name,
                )
            );


            my $para =
              $context->append_element( odf_create_paragraph() );

            my $list = $para->insert_element(
                    odf_create_list, position => NEXT_SIBLING
                    );

            my $size = length ${$candidate}{'DNASequence'};
            my $vienna_seq =
              PipelineMiRNA::Candidate->make_Vienna_viz( ${$candidate}{'Vienna'},
                ${$candidate}{'DNASequence'} );

            $list->add_item(text => "Name: ${$candidate}{'name'}", style => "Standard");
            $list->add_item(text => "Position: ${$candidate}{'position'} ($size nt)", style => "Standard");
            $list->add_item(text => "Strand:", style => "Standard");

            my $subtext = "";
            if(${$candidate}{'Vienna'} ne ${$candidate}{'Vienna_optimal'}){
                $subtext .= ""
            } else {
                $subtext.= "(This stem-loop structure is the MFE structure)"
            }
            my $item = $list->add_item(text => "Sequence and stem-loop structure: ", style => "Standard");
#            $item->add(text => "Sequence and stem-loop structure:\n$vienna_seq \n$subtext", style => "Monospace");
            $item->append_element(
            odf_create_paragraph(
                text    => $vienna_seq,
                style   =>'Monospace'
                )
            );
            my $para1 = $item->append_element(
            odf_create_paragraph(
                text    => $subtext,
                style   =>'Standard'
                )
            );
            $para1->set_span(filter  => "structure", style   => 'StandardBold');

            # Section secondary structure

            $context->append_element(
                odf_create_heading(
                    level => 3,
                    text  => "Secondary structure",
                )
            );

            # Copying the image
            my $img_path      = ${$candidate}{'image'};
            my $img_full_path = $img_path;
            my $new_img_path = File::Spec->catfile($images_dir, "$key.png");
            copy($img_full_path, $new_img_path)
                or die "Copy of $img_full_path to $images_dir failed: $!";

            my ( $lien_image, $taille_image ) =
              $doc->add_image_file($new_img_path);

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

           # Section secondary structure
            $context->append_element(
                odf_create_heading(
                    level => 3,
                    text  => "Thermodynamics stability",
                )
            );

            $para =
              $context->append_element( odf_create_paragraph() );

            $list = $para->insert_element(
                    odf_create_list, position => NEXT_SIBLING
                    );

            # TODO: maybe we do not have those ; infer that from  run_options config file
            $list->add_item(text => "MFE: ${$candidate}{'mfe'} kcal/mol", style => "Standard");
            $list->add_item(text => "AMFE: ${$candidate}{'amfe'}", style => "Standard");
            $list->add_item(text => "MFEI: ${$candidate}{'mfei'}", style => "Standard");


            # Section Mirbase alignments
            $self->add_ODF_alignments($context, $candidate);
        }   # if key in tab
    }    #  while each %results

    # save the generated document and quit
    $doc->save( target => $ODT_abspath, pretty => TRUE );
    return $ODT_abspath;
}

=method add_ODF_alignments

=cut

sub add_ODF_alignments{
    my ( $self, @args ) = @_;
    my $context = shift @args;
    my %candidate = %{shift @args};
    my %alignments = %{$candidate{'alignments'}};
    my %mirdup_results = %{$candidate{'mirdup_validation'}};

    $context->append_element(
        odf_create_heading(
            level => 3,
            text  => "miRBase alignments",
        )
    );

    my @TOC;
    my $predictionCounter = 0;

    # Sorting by position
    my @keys = sort { ( PipelineMiRNA::Utils::get_element_of_split($a, '-', 0)  <=>
                        PipelineMiRNA::Utils::get_element_of_split($b, '-', 0)
                      ) ||
                      ( PipelineMiRNA::Utils::get_element_of_split($a, '-', 1)  <=>
                        PipelineMiRNA::Utils::get_element_of_split($b, '-', 1))
                    } keys %alignments;
    foreach my $position (@keys) {
        my ($left, $right) = split(/-/, $position);

        # MiRdup
#        my $mirdup_key = $dir . '__' . $position;
#        my $mirdup_prediction;
#        if ( $mirdup_results{$mirdup_key} ){
#            $mirdup_prediction = 'This prediction is validated by MiRdup';
#        } else {
#            $mirdup_prediction = 'This prediction is not validated by MiRdup';
#        }

        $predictionCounter += 1;

        # Sorting the hit list by descending value of the 'score' element
        my @hits = sort { $b->{'score'} <=> $a->{'score'} } @{$alignments{$position}};
        my $title = "Prediction $predictionCounter: $position";

        $context->append_element(
            odf_create_heading(
                level => 4,
                text  => $title,
            )
        );
        $context->append_element(
            odf_create_paragraph(
                text    => $candidate{'hairpin'},
                style   =>'Hairpin'
            )
        );

        $context->append_element(
            odf_create_heading(
                level => 5,
                text  => 'Alignments',
            )
        );

        foreach my $hit (@hits){
            my $alignment = $hit->{'alignment'};
            my $names = $hit->{'name'} . q{ } . $hit->{'def_query'};

            my $name;
            my $html_name;
            my @splitted = split('\|', $names);

            my $spacing = 15;
            my ($top, $middle, $bottom) = split(/\n/, $alignment);
            $top    = sprintf "%-${spacing}s %3s %s %s", 'query', $hit->{'begin_target'}, $top,   $hit->{'end_target'};
            $middle = sprintf "%-${spacing}s %3s %s %s", '',      '',                     $middle, '';

            my $title = '';
            if( (scalar @splitted) > 1 ) {
                $title = 'miRBase sequences: ';
            }else{
                $title = 'miRBase sequence: ';
            }

            $name = "miRBase";
            my @sequences;

            my @mirbase_links;
            my @mirbase_ids;
            foreach my $seq (@splitted){
                $seq =~ s/^\s+//;
                $seq =~ s/\s+$//;
                if ($seq =~ 'revcomp'){
                } else {
                    my @splitted_one = split(/ /, $seq);
                    my $name = $splitted_one[0];
                    my $mirbase_id = $splitted_one[1];
                    my $mirbase_link = PipelineMiRNA::Utils::make_mirbase_link($mirbase_id);
                    push @mirbase_links, $mirbase_link;
                    push @mirbase_ids, $mirbase_id;
                }
            }
            my $para_seqs = odf_create_paragraph(
                text    => $title . join(', ', @mirbase_ids),
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


=method get_ODF_path

Return the paths to the ODT document

=cut

sub get_ODF_path{
    my ( $self, @args ) = @_;
    my $jobId = shift @args;
    my $ODT_filename = "Prediction_report_$jobId.odt";
    my $jobPath = PipelineMiRNA::Results->jobId_to_jobPath($jobId);
    my $ODT_abspath = File::Spec->catfile( $jobPath, $ODT_filename );
    my $ODT_serverpath = PipelineMiRNA::Paths->filesystem_to_relative_path($ODT_abspath);
    return ($ODT_abspath, $ODT_serverpath);
}

=method get_report

Gererates the report if it does not exist already,
and return the server path to it.

=cut

sub get_report {
    my ( $self, @args ) = @_;
    my $jobId = shift @args;
    my @sequences_to_export = shift @args;
    my ($ODT_abspath, $ODT_serverpath) = $self->get_ODF_path($jobId);
#    if (! -e $ODT_abspath){
#    }
    my $path = $self->generate_report($jobId, \@sequences_to_export);
    return $ODT_serverpath;
}

1;
