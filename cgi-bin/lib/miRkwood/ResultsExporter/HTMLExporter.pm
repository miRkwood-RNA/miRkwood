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
    my $output   = '<tr>';
    my $anchor   = "${$candidate}{'name'}-${$candidate}{'position'}";
    my $contents = "<a href='#$anchor'>${$candidate}{'name'}</a>";
    $output .= "<td id='table_$anchor'>$contents</td>\n";

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
            $contents = '<td><center><font color="#FF8000">';
            for (my $i = 0; $i < ${$candidate}{'quality'}; $i++){
                $contents .= '&#x2605;';
            }
            $contents .= '</font></center></td>';
        }
        elsif ($header eq 'reads_distribution'){
            $contents = "<td onclick=\"location.href='$reads_file'\" style='cursor:pointer;'><center><font color='#FF8000'>";
            for (my $i = 0; $i < ${$candidate}{'reads_distribution'}; $i++){
                $contents .= '&#x2605;';
            }
            $contents .= '</font></center></td>';
        }
        elsif ($header eq 'alignment'){
            $contents = '<td><center><font color="#008000">';
            for (my $i = 0; $i < ${$candidate}{'alignment'}; $i++){
                $contents .= '&#x2713;';
            }
            $contents .= '</font></center></td>';
        }
        elsif ($header eq 'mfei'){
            $contents = miRkwood::Utils::restrict_num_decimal_digits($contents, 3);
            if ( $contents < -0.8 ){
                $contents = '<font color="#FF00FF">' . $contents . '</font>';
            }
            $contents = "<td>$contents</td>";
        }
        elsif ($header eq 'mfe' or $header eq 'amfe'){
            $contents = '<td>' . miRkwood::Utils::restrict_num_decimal_digits($contents, 3) . '</td>';
        }
        elsif ( $header eq 'position'){
            $contents = "<td><a href='#$anchor'>${$candidate}{$header}</a></td>\n";
        }
        elsif ( $header eq 'nb_reads' ){
            if ( ${$candidate}{'criteria_nb_reads'} ){
                $contents = "<font color='#FF00FF'>${$candidate}{$header}</font>";
            }
            $contents = "<td onclick=\"location.href='$reads_file'\" style='cursor:pointer;'>$contents</td>";
        }
        else {
            $contents = "<td>$contents</td>";
        }
        if ( !defined $contents ) {
            $contents = '<td></td>';
        }
        $output .= "$contents\n";
    }
    $output .= "\n</tr>\n";
    return $output;
}

sub perform_export{
    my ( $self, @args ) = @_;

    my %results = %{$self->{'results'}};
    my @keys = $self->get_sorted_keys();
    my $nb_results = scalar( @keys );

    my $output = '';
    $output .= "<h3>$nb_results candidates found.</h3>\n";

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
