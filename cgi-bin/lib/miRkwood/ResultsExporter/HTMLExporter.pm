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
    my @headers = $self->get_headers();
    for my $header ( @headers ) {
        my $td_content = '';
        $contents = ${$candidate}{$header};
        if ($header eq 'quality'){
            $contents = '<center><font color="#FF8000">';
            for (my $i = 0; $i < ${$candidate}{'quality'}; $i++){
                $contents .= '&#x2605;';
            }
            $contents .= '</font></center>';
        }
        elsif ($header eq 'reads_distribution'){
            $contents = '<center><font color="#FF8000">';
            for (my $i = 0; $i < ${$candidate}{'reads_distribution'}; $i++){
                $contents .= '&#x2605;';
            }
            $contents .= '</font></center>';
        }
        elsif ($header eq 'alignment'){
            $contents = '<center><font color="#008000">';
            for (my $i = 0; $i < ${$candidate}{'alignment'}; $i++){
                $contents .= '&#x2713;';
            }
            $contents .= '</font></center>';
        }
        elsif ($header eq 'mfei'){
            $contents = miRkwood::Utils::restrict_num_decimal_digits($contents, 3);
            if ( $contents < -0.8 ){
                $contents = '<font color="#FF00FF">' . $contents . '</font>';
            }
        }
        elsif ($header eq 'mfe' or $header eq 'amfe'){
            $contents = miRkwood::Utils::restrict_num_decimal_digits($contents, 3);
        }
        elsif ( $header eq 'position'){
            $contents = "<a href='#$anchor'>${$candidate}{$header}</a>\n";
        }
        elsif ( $header eq 'nb_reads' ){
            my $mirna_type = 'novel_miRNA';
            if ( defined(${$candidate}{'mirbase_id'}) ){
                $mirna_type = 'known_miRNA';
            }
            my $reads_path = File::Spec->catdir( File::Spec->updir(), File::Spec->updir(), miRkwood::Paths::get_reads_dir_name(), $mirna_type);
            my $reads_file = File::Spec->catfile( $reads_path, ${$candidate}{'identifier'} . '.txt' );
            $contents = "<a href='$reads_file' class='nodecoration'>${$candidate}{$header}</a>";
            if ( ${$candidate}{'criteria_nb_reads'} ){
                $contents = "<a href='$reads_file' class='nodecoration'><font color='#FF00FF'>${$candidate}{$header}</font></a>";
            }
        }
        if ( !defined $contents ) {
            $contents = q{};
        }
        $output .= "<td>$contents</td>\n";
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
