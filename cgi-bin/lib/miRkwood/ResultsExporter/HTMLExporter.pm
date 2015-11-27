package miRkwood::ResultsExporter::HTMLExporter;

# ABSTRACT: Class for exporting results in HTML table

use strict;
use warnings;

use parent 'miRkwood::ResultsExporter::ResultsExporter';

use miRkwood::Candidate;
use miRkwood::Utils;


sub get_header {
    my ( $self, @args ) = @_;
    my $output = '<tr>';
    for my $header ( ('name'), $self->get_headers() ) {
        $output .= "<th>$header</th>\n";
    }
    $output .= "</tr>\n";
    return $output;
}

sub export_candidate {
    my ( $self, @args ) = @_;
    my $candidate = shift @args;
    my $output = "<tr>\n";
    my $non_clickable_cell = "class='non_clickable_cell'";
    my $clickable_cell = "class='clickable_cell'";
    my $read_cell = "class='read_cell'";
    my $star_cell = "class='star_cell'";
    my $alignment_cell = "class='alignment_cell'";
    my $anchor   = "${$candidate}{'name'}-${$candidate}{'position'}";
    my $contents = "<a href='#$anchor' class='nodecoration'>${$candidate}{'name'}</a>";

    $output .= "<td $clickable_cell id='table_$anchor'>$contents</td>\n";

    my $cfg    = miRkwood->CONFIG();

    my $mirna_type = 'novel_miRNA';
    if ( defined(${$candidate}{'mirbase_id'}) ){
        $mirna_type = 'known_miRNA';
    }
    my $reads_path = File::Spec->catdir( File::Spec->updir(), File::Spec->updir(), miRkwood::Paths::get_reads_dir_name(), $mirna_type);
    my $reads_file = File::Spec->catfile( $reads_path, ${$candidate}{'identifier'} . '.txt' );

    my @headers = $self->get_headers();
    for my $header ( @headers ) {
        $contents = ${$candidate}{$header};
        if ($header eq 'quality'){
            $contents = "<td $non_clickable_cell><font color='#FFBE00'>";
            for (my $i = 0; $i < ${$candidate}{'quality'}; $i++){
                $contents .= '&#x2605;';
            }
            $contents .= '</font></td>';
        }
        elsif ($header eq 'reads_distribution'){
            if ( ${$candidate}{'reads_distribution'} == 0){
                $contents = "<td $non_clickable_cell> </td>";
            }
            else{
                $contents = "<td $star_cell><a href='$reads_file' class='nodecoration'><font color='#FFBE00'>";
                for (my $i = 0; $i < ${$candidate}{'reads_distribution'}; $i++){
                    $contents .= '&#x2605;';
                }
                $contents .= '</font></a></td>';
            }
        }
        elsif ($header eq 'alignment'){
            $contents = "<td ";
            my $alignment_file = '';
            if ( $cfg->param('job.pipeline') eq 'abinitio' ){
                $alignment_file = File::Spec->catfile( File::Spec->updir(),
                                                         miRkwood::Paths::get_alignments_dir_name(),
                                                         ${$candidate}{'identifier'}. '_aln.txt' );
            }
            elsif ( $cfg->param('job.pipeline') eq 'smallRNAseq' ){
                $alignment_file = File::Spec->catfile( File::Spec->updir(),
                                                         File::Spec->updir(),
                                                         miRkwood::Paths::get_alignments_dir_name(),
                                                         ${$candidate}{'identifier'}. '_aln.txt' );
            }
            if ( ${$candidate}{'alignment'} > 0){
                $contents .= "$alignment_cell><a href='$alignment_file' class='nodecoration'><font color='#41BE47'>";
                for (my $i = 0; $i < ${$candidate}{'alignment'}; $i++){
                    $contents .= '&#x2713;';
                }
                $contents .= '</font></a></td>';
            }
            else {
                $contents .= "$non_clickable_cell></td>";
            }
        }
        elsif ($header eq 'mfei'){
            $contents = miRkwood::Utils::restrict_num_decimal_digits($contents, 3);
            if ( $contents < -0.8 ){
                $contents = '<font color="#E13EA5">' . $contents . '</font>';
            }
            $contents = "<td $non_clickable_cell>$contents</td>";
        }
        elsif ($header eq 'mfe' or $header eq 'amfe'){
            $contents = "<td $non_clickable_cell>" . miRkwood::Utils::restrict_num_decimal_digits($contents, 3) . '</td>';
        }
        elsif ( $header eq 'position'){
            my ($start, $end) = split( /-/, $contents);
            $start = miRkwood::Utils::make_numbers_more_readable( $start );
            $end = miRkwood::Utils::make_numbers_more_readable( $end );
            $contents = "<td $clickable_cell><a href='#$anchor' class='nodecoration'>$start-$end</a></td>";
        }
        elsif ( $header eq 'nb_reads' ){
            if ( ${$candidate}{'criteria_nb_reads'} ){
                $contents = "<font color='#12D0E5'>${$candidate}{$header}</font>";
            }
            $contents = "<td $read_cell><a href='$reads_file' class='nodecoration'>$contents</a></td>";
        }
        else {
            $contents = "<td $non_clickable_cell>$contents</td>";
        }
        if ( !defined $contents ) {
            $contents = "<td $non_clickable_cell></td>";
        }
        $output .= "$contents\n";
    }
    $output .= "</tr>\n\n";
    return $output;
}

sub perform_export{
    my ( $self, @args ) = @_;

    my %results = %{$self->{'results'}};
    my @keys = $self->get_sorted_keys();
    my $nb_results = scalar( @keys );

    my $output = '';

    $output .= "<table>\n<tbody>";

    $output .= $self->get_header();

    foreach my $key (@keys) {
        if ( $self->is_sequence_to_export($key)){
            my $candidate = $results{$key};
            $output .= $self->export_candidate($candidate);
        }
    }
    $output .= "</tbody>\n</table>\n";
    return $output;
}

1;
