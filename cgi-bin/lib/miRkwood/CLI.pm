package miRkwood::CLI;

# ABSTRACT: Code for the command line interface

use strict;
use warnings;
use File::Spec;

use miRkwood;
use miRkwood::Paths;
use miRkwood::Results;
use miRkwood::ResultsExporterMaker;

=method process_results_dir_for_offline

Process the results in the given directory for offline use.

Usage:
  miRkwood::CLI::process_results_dir_for_offline($folder);

=cut

sub process_results_dir_for_offline {
    my @args              = @_;
    my $abs_output_folder = shift @args;
    my $pipeline_type     = shift @args;
    my $mirna_type        = shift @args;

# In debug mode (without executing the pipeline), we need to set the config file
#miRkwood->CONFIG_FILE(miRkwood::Paths->get_job_config_path( $output_folder )); -> not sure this is still working...

    my $candidates_dir = '';

    my $final_results_folder = miRkwood::Paths::get_results_folder_for_CLI_from_job_dir( $abs_output_folder, $pipeline_type, $mirna_type );

    if ( $pipeline_type eq 'smallRNAseq' ){
        if ( $mirna_type eq 'known_miRNA' ){
            $candidates_dir = miRkwood::Paths::get_known_candidates_dir_from_job_dir( $abs_output_folder );
        }
        else{
            $candidates_dir = miRkwood::Paths::get_new_candidates_dir_from_job_dir( $abs_output_folder );
        }
    }
    else{
        $candidates_dir = miRkwood::Paths::get_dir_candidates_path_from_job_dir( $abs_output_folder );
    }

    my $tmp_pieces_folder = miRkwood::Paths::create_folder(
            File::Spec->catdir( $final_results_folder, miRkwood::Paths::get_pieces_folder_basename_for_CLI() )
        );

    my %results = miRkwood::Results->deserialize_results($candidates_dir);

    my $html = make_html_from_results( \%results, $abs_output_folder, $pipeline_type, $mirna_type );

    my $html_output = 'results.html';
    if ( $pipeline_type ne 'abinitio' ){
        $html_output = "results_$mirna_type.html";
    }

    my $html_page = File::Spec->catfile( $final_results_folder, $html_output );
    open( my $HTML, '>', $html_page )
      or die("Cannot open $html_page: $!");
    print {$HTML} $html
      or die("Cannot write in $html_page: $!");
    close($HTML)
      or die("Cannot close $html_page: $!");
    return;
}

=method make_html_from_results

Given a reference to a results hash, makes the HTML
Usage:
  my $html = make_html_from_results( \%results, $output_folder );

=cut

sub make_html_from_results {
    my @args    = @_;
    my $results           = shift @args;
    my $abs_output_folder = shift @args;
    my $pipeline_type     = shift @args;
    my $mirna_type        = shift @args;
    my $pieces_folder = miRkwood::Paths::get_pieces_folder_basename_for_CLI();

    my %results = %{$results};
    my ($css) = get_page_css();
    my $page = '<h2>Overview of results</h2>';
    $page .= make_all_exports( \%results, $abs_output_folder, $pieces_folder, $pipeline_type, $mirna_type );

    my $exporter = miRkwood::ResultsExporterMaker->make_html_results_exporter( $pipeline_type, $mirna_type );
    $exporter->initialize('', \%results);
    $page .= $exporter->perform_export();

    my @keys = sort {
        ( $results{$a}->{'name'} cmp $results{$b}->{'name'} )
          || (
            $results{$a}->{'start_position'} <=> $results{$b}->{'start_position'} )
    } keys %results;

    foreach my $key ( @keys ){
        $page .= make_candidate_page( $results{ $key }, $pieces_folder, $abs_output_folder, $pipeline_type, $mirna_type );
    }
    my $html = get_simple_results_page( $page, $css );

    return $html;
}

=method make_all_exports

Given a reference to a results hash, generates the various
exports in the given output directory.

Usage:
  my $html = make_all_exports(\%results, $output_folder);

=cut

sub make_all_exports {
    my (@args)        = @_;
    my $results_ref       = shift @args;
    my $abs_output_folder = shift @args;
    my $pieces_folder     = shift @args;
    my $pipeline_type     = shift @args;
    my $mirna_type        = shift @args;
    my $id_job = '';
    my $exporter;
    my $html = '<h3>Get results as</h3> <ul>';

    my $final_results_folder = miRkwood::Paths::get_results_folder_for_CLI_from_job_dir( $abs_output_folder, $pipeline_type, $mirna_type );

    $exporter = miRkwood::ResultsExporterMaker->make_csv_results_exporter( $pipeline_type, $mirna_type );
    $exporter->initialize($id_job, $results_ref);
    $exporter->export_on_disk( $final_results_folder );
    my $csv_file = File::Spec->catfile($final_results_folder, $exporter->get_filename());
    $html .= '<li><a href="' . $exporter->get_filename() . '">tab-delimited format (csv)</a></li>';

    $exporter = miRkwood::ResultsExporterMaker->make_fasta_results_exporter( $mirna_type );
    $exporter->initialize($id_job, $results_ref);
    $exporter->export_on_disk( $final_results_folder );
    my $fasta_file = File::Spec->catfile($final_results_folder, $exporter->get_filename());
    $html .= '<li><a href="' . $exporter->get_filename() . '">Fasta</a></li>';

    $exporter = miRkwood::ResultsExporterMaker->make_dotbracket_results_exporter( $mirna_type );
    $exporter->initialize($id_job, $results_ref);
    $exporter->export_on_disk( $final_results_folder );
    my $dotbracket_file = File::Spec->catfile($final_results_folder, $exporter->get_filename());
    $html .= '<li><a href="' . $exporter->get_filename() . '">dot-bracket format (plain sequence + secondary structure)</a></li>';

    $exporter = miRkwood::ResultsExporterMaker->make_gff_results_exporter( $mirna_type );
    $exporter->initialize($id_job, $results_ref);
    $exporter->export_on_disk( $final_results_folder );
    my $gff_file = File::Spec->catfile($final_results_folder, $exporter->get_filename());
    $html .= '<li><a href="' . $exporter->get_filename() . '">gff format</a></li>';

    $html .= '</ul>';
    return $html;
}

=method make_candidate_page

Given a candidate hash, make the HTML page

Usage:
  my $html = make_candidate_page( \$candidate, $pieces_folder, $output_folder );

=cut

sub make_candidate_page {
    my (@args)        = @_;
    my $candidate     = shift @args;
    my $pieces_folder = shift @args;
    my $abs_output_folder = shift @args;
    my $pipeline_type = shift @args;
    my $mirna_type    = shift @args;

    my $output_folder = miRkwood::Paths::get_results_folder_for_CLI_from_job_dir( $abs_output_folder, $pipeline_type, $mirna_type );

    my $size = length $candidate->{'sequence'};
    my $cfg    = miRkwood->CONFIG();
    my $candidate_name = $candidate->get_shortened_name();

    my $star = '<font color=\'#FF8000\'>&#x2605;</font>';
    my $coche = '<font color=\'#008000\'>&#x2713;</font>';
    my $arrow = '<font color=\'#008000\'>&#x25B2;</font>';
    #~ my $arrow = '<font color=\'#008000\'>&uarr;</font>';

    ### make files in pieces folder
    my $candidate_fasta_file =
      File::Spec->catfile( $pieces_folder, "$candidate_name.fa" );
    open( my $FASTA_FILE,
        '>', File::Spec->catfile( $output_folder, $candidate_fasta_file ) )
      or die("Cannot open file $candidate_fasta_file: $!");
    print {$FASTA_FILE}
      $candidate->candidateAsFasta()
      or die("Cannot write in file $candidate_fasta_file: $!");
    close($FASTA_FILE)
      or die("Cannot close file $candidate_fasta_file: $!");

    my $vienna_file =
      File::Spec->catfile( $pieces_folder, "$candidate_name.txt" );
    open( my $VIENNA_FILE,
        '>', File::Spec->catfile( $output_folder, $vienna_file ) )
      or die("Cannot open $vienna_file: $!");
    print {$VIENNA_FILE}
      $candidate->candidateAsVienna(0)
      or die("Cannot write in file $vienna_file: $!");
    close($VIENNA_FILE)
      or die("Cannot close file $vienna_file: $!");

    my $vienna_file_optimal =
      File::Spec->catfile( $pieces_folder, $candidate_name . '_optimal.txt' );
    open( my $VIENNA_FILE_OPT,
        '>', File::Spec->catfile( $output_folder, $vienna_file_optimal ) )
      or die("Cannot open $vienna_file_optimal: $!");
    print {$VIENNA_FILE_OPT}
      $candidate->candidateAsVienna(1)
      or die("Cannot write in file $vienna_file_optimal: $!");
    close($VIENNA_FILE_OPT)
      or die("Cannot close file $vienna_file_optimal: $!");

    my $alternatives_file = File::Spec->catfile( $pieces_folder,
        $candidate_name . '_alternatives.txt' );
    open( my $ALT_FILE, '>',
        File::Spec->catfile( $output_folder, $alternatives_file ) )
      or die("Cannot open $alternatives_file: $!");
    print {$ALT_FILE}
      $candidate->alternativeCandidatesAsVienna()
      or die("Cannot write in file $alternatives_file: $!");
    close($ALT_FILE)
      or die("Cannot close file $alternatives_file: $!");


    ### links
    my $linkFasta         = "$candidate_fasta_file";
    my $linkVienna        = "$vienna_file";
    my $linkAlternatives  = "$alternatives_file";
    my $linkViennaOptimal = "$vienna_file_optimal";


    ### HTML sections
    my $mfei = $candidate->{'mfei'};
    if ( $candidate->{'mfei'} < -0.8 ){
        $mfei = "<font color='#FF00FF'>$candidate->{'mfei'}</font>";
    }

    # Reads
    my $mirna_html = '';
    my $reads_html = '';
    if ( $pipeline_type eq 'smallRNAseq' ){
        # miRNA sequence
        my $mirna_sequence = '';
        if ( $candidate->{'mirna_length'} eq '' ||  $candidate->{'mirna_length'} eq '0' ){
            $mirna_sequence = 'None';
        }
        else{
            $mirna_sequence = "$candidate->{'mirna_sequence'} ($candidate->{'mirna_length'} nt)";
        }
        $mirna_html = "<li><b>miRNA sequence:</b> $mirna_sequence</li>";

        # Reads
        # this is not very robust. be cautious if you change the tree view
        my $reads_path = File::Spec->catdir( File::Spec->updir(), File::Spec->updir(), miRkwood::Paths::get_reads_dir_name(), $mirna_type);
        my $reads_file = File::Spec->catfile( $reads_path, $candidate->{'identifier'} . '.txt' );
        my $reads_score = '';
        my $absolute_read_cloud_path = File::Spec->catfile(
                miRkwood::Paths::get_dir_reads_path_from_job_dir( $cfg->param('job.directory') ),
                $mirna_type,
                $candidate->{'identifier'}.'.txt');
        my $read_cloud  = include_read_cloud_in_html( $absolute_read_cloud_path );
        my $read_duplex = '';

        my $nb_reads = $candidate->{'nb_reads'};
        if ( $candidate->{'criteria_nb_reads'} == 1 ){
            $nb_reads = "<font color='#FF00FF'>$candidate->{'nb_reads'}</font>";
        }

        if ( ! defined( $candidate->{'mirbase_id'} ) ){
            # reads distribution
            if ( $candidate->{'criteria_reads_mirna'} == 1 && $candidate->{'criteria_star'} == 1 ){
                $reads_score = "<li><b>Distribution of reads:</b> two islands $star$star</li>";
            }
            if ( ($candidate->{'criteria_reads_mirna'} + $candidate->{'criteria_star'}) == 1 ){
                $reads_score = "<li><b>Distribution of reads:</b> one island $star</li>";
            }
            if ( $candidate->{'criteria_reads_mirna'} == 0 && $candidate->{'criteria_star'} == 0 ){
                $reads_score = '<li><b>Distribution of reads:</b> random</li>';
            }

            # stability of duplex
            $read_duplex = '<li><b>Stability of the miRNA duplex (mirdup):</b> ';
            if ( $candidate->{'criteria_mirdup'} == 1 ){
                $read_duplex .= "yes $star</li>";
            }
            else {
                $read_duplex .= 'no </li>';
            }
        }

        $reads_html = <<"END_TXT";
<li>
    <b>Number of reads:</b> $nb_reads (<a href='$reads_file'>download<a/>)
</li>
$reads_score
$read_cloud
$read_duplex
END_TXT
    }

    # miRBase alignments
    my $alignments_html = '';
    if ( ! defined( $candidate->{'mirbase_id'} ) ){
        $alignments_html .= '<li><b>miRBase alignment:</b> ';
        if ( $candidate->{'alignment'} == 0 ){
             $alignments_html .= 'none </li>';
        }
        else {
            if ( $candidate->{'alignment'} == 2 ){
                $alignments_html .= "$coche$coche presence of alignments that cover the miRNA locus (see reads cloud above)";
            }
            elsif ( $candidate->{'alignment'} == 1 ){
                $alignments_html .= "$coche presence of alignments, which do not overlap the miRNA locus (see reads cloud above)";
            }
            else {
                $alignments_html .= 'none';
            }
            my $alignment_file = File::Spec->catfile( File::Spec->updir(),
                                                      File::Spec->updir(),
                                                      miRkwood::Paths::get_alignments_dir_name(),
                                                      $candidate->{'identifier'}. '_aln.txt' );

            $alignments_html .= $candidate->include_alignments_in_html();
        }
        $alignments_html .= '</li>';
    }

    # Alternative sequences
    my $alternatives_HTML = '';
    if( $candidate->{'alternatives'} and scalar(keys%{ $candidate->{'alternatives'} }) ){
        $alternatives_HTML = "[<a href='$linkAlternatives'>alternative sequences</a>]"
    }

    ### make page
    my $html = <<"END_TXT";
<h2 id='$candidate->{'name'}-$candidate->{'position'}'><a href='#table_$candidate->{'name'}-$candidate->{'position'}' class='nodecoration'>Results for $candidate->{'name'}, $candidate->{'position'} ($candidate->{'strand'}) $arrow</a></h2>

    <ul>
        <li>
            <b>Name: </b>$candidate->{'name'}
        </li>
        <li>
            <b>Position:</b> $candidate->{'position'} ($size nt)
        </li>
        <li>
            <b>Strand:</b> $candidate->{'strand'}
        </li>
        $mirna_html
        <li>
            <b>miRNA precursor:</b> [<a href='$linkFasta'>FASTA sequence</a>] 
                                    [<a href='$linkVienna'>stem-loop structure</a>] 
                                    [<a href='$linkViennaOptimal'>optimal MFE structure</a>] 
                                    $alternatives_HTML
        </li>
        <li>
            <b>Stability of the secondary structure of the precursor:</b> <i>MFE</i> $candidate->{'mfe'} kcal/mol | 
                                                                          <i>AMFE</i> $candidate->{'amfe'} | 
                                                                          <i>MFEI</i> $mfei
        </li>
        $reads_html
        $alignments_html
    </ul>
END_TXT


    return $html;
}

=method include_read_cloud_in_html

  Read contents of read cloud file and store it
  in 'pre' tags, except the header.

=cut
sub include_read_cloud_in_html {
    my @args = @_;
    my $read_cloud_file = shift @args;
    my $result = '<pre>';
    open(my $IN, '<', $read_cloud_file) or return '';
    my $line = <$IN>;
    while ( <$IN> ){
        chomp;
        if ( $_ !~ /^Chromosome/ && $_ !~ /^Position/ && $_ !~ /^Strand/ && $_ ne '' ){
            $result .= "$_\n";
        }
    }
    close $IN;
    $result .= '</pre>';
    return $result;
}


=method get_simple_results_page

Make a simple HTML page with the given body and CSS.

Usage:
  my $html = get_simple_results_page( $page, $css );

=cut

sub get_simple_results_page {
    my @args = @_;
    my $page = shift @args;
    my $css  = shift @args;
    my $HTML = <<"END_TXT";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <title>miRkwood - MicroRNA identification</title>
        <STYLE type="text/css">$css</STYLE>
    </head>
    <body>
        $page
    </body>
</html>
END_TXT
    return $HTML;
}


=method get_page_css

Returns the CSS needed for the webpage

=cut

sub get_page_css {
    my @args = @_;
    my $css = <<"END_TXT";
table{
border:1px solid black;
border-collapse:collapse;
width:80%;
}
th, td {
border:1px solid black;
}
span.mature {
    color: blue;
}
a.nodecoration
{ 
    text-decoration:none;
    color: black;
    display:block;
    width:100%; 
    height: 100%;
}
END_TXT
    return ($css);
}

1;
