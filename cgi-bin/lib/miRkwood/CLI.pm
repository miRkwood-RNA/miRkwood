package miRkwood::CLI;

# ABSTRACT: Code for the command line interface

use strict;
use warnings;
use File::Spec;

use miRkwood;
use miRkwood::Paths;
use miRkwood::Utils;
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

    my $abs_sequences_folder = miRkwood::Paths::create_folder(
            File::Spec->catdir( $abs_output_folder, miRkwood::Paths::get_sequences_folder_basename_for_CLI() )
        );

    my $basic_results_file = File::Spec->catfile($abs_output_folder, 'basic_candidates.yml');

    if ( $pipeline_type eq 'smallRNAseq' ){
        if ( $mirna_type eq 'known_miRNA' ){
            $candidates_dir = miRkwood::Paths::get_known_candidates_dir_from_job_dir( $abs_output_folder );
            miRkwood::Paths::create_folder( File::Spec->catdir( $abs_sequences_folder, miRkwood::Paths::get_basename_for_known_miRNA() ) );
            $basic_results_file = File::Spec->catfile($abs_output_folder, 'basic_known_candidates.yml');
        }
        else{
            $candidates_dir = miRkwood::Paths::get_new_candidates_dir_from_job_dir( $abs_output_folder );
            miRkwood::Paths::create_folder( File::Spec->catdir( $abs_sequences_folder, miRkwood::Paths::get_basename_for_novel_miRNA() ) );
        }
    }
    else{
        $candidates_dir = miRkwood::Paths::get_dir_candidates_path_from_job_dir( $abs_output_folder );
    }

    my %results = miRkwood::Results->deserialize_results($basic_results_file);

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
  my $html = make_html_from_results( \%results, $abs_output_folder, $pipeline_type, $mirna_type );

=cut

sub make_html_from_results {
    my @args    = @_;
    my $results           = shift @args;
    my $abs_output_folder = shift @args;
    my $pipeline_type     = shift @args;
    my $mirna_type        = shift @args;
    my $sequences_folder = miRkwood::Paths::get_sequences_folder_basename_for_CLI();

    my %results = %{$results};
    my ($css) = get_page_css();
    my $page = "<h2>Overview of results</h2>\n";

    my $nb_results = scalar( keys%results );
    $nb_results = miRkwood::Utils::make_numbers_more_readable( $nb_results );
    $page .= "<h3>$nb_results candidates found.</h3>\n";

    $page .= make_all_exports( \%results, $abs_output_folder, $pipeline_type, $mirna_type );
    $page .= '<br />';

    my $exporter = miRkwood::ResultsExporterMaker->make_html_results_exporter( $pipeline_type, $mirna_type );
    $exporter->initialize('', \%results);
    $page .= $exporter->perform_export();

    my @keys = sort {
        ( $results{$a}->{'name'} cmp $results{$b}->{'name'} )
          || (
            $results{$a}->{'start_position'} <=> $results{$b}->{'start_position'} )
    } keys %results;

    foreach my $key ( @keys ){
        $page .= make_candidate_page( $results{ $key }, $sequences_folder, $abs_output_folder, $pipeline_type, $mirna_type );
    }
    my $html = get_simple_results_page( $page, $css );

    return $html;
}

=method make_all_exports

Given a reference to a results hash, generates the various
exports in the given output directory.

Usage:
  my $html = make_all_exports(\%results, $abs_output_folder, $pipeline_type, $mirna_type);

=cut

sub make_all_exports {
    my (@args)        = @_;
    my $results_ref       = shift @args;
    my $abs_output_folder = shift @args;
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

    $exporter = miRkwood::ResultsExporterMaker->make_org_results_exporter( $mirna_type );
    $exporter->initialize($id_job, $results_ref);
    $exporter->export_on_disk( $final_results_folder );
    my $org_file = File::Spec->catfile($final_results_folder, $exporter->get_filename());
    $html .= '<li><a href="' . $exporter->get_filename() . '">ORG format</a></li>';

    $html .= '</ul>';
    return $html;
}

=method make_candidate_page

Given a candidate hash, make the HTML page

Usage:
  my $html = make_candidate_page( \$candidate, 
      $sequences_folder_folder,
      $abs_output_folder,
      $pipeline_type,
      $mirna_type );

=cut

sub make_candidate_page {
    my (@args)        = @_;
    my $candidate     = shift @args;
    my $sequences_folder = shift @args;
    my $abs_output_folder = shift @args;
    my $pipeline_type = shift @args;
    my $mirna_type    = shift @args;
    my $relative_sequences_folder = File::Spec->catdir( File::Spec->updir(), File::Spec->updir(), $sequences_folder);
    my $absolute_sequences_folder = File::Spec->catfile( $abs_output_folder, $sequences_folder );
    if ( $pipeline_type eq 'smallRNAseq' ){
        if ( $mirna_type eq 'known_miRNA' ){
            $relative_sequences_folder = File::Spec->catdir( 
                $relative_sequences_folder,
                miRkwood::Paths::get_basename_for_known_miRNA()
            );
            $absolute_sequences_folder = File::Spec->catdir( 
                $absolute_sequences_folder,
                miRkwood::Paths::get_basename_for_known_miRNA()
            );
        }
        else{
            $relative_sequences_folder = File::Spec->catdir(
                $relative_sequences_folder,
                miRkwood::Paths::get_basename_for_novel_miRNA()
            );
            $absolute_sequences_folder = File::Spec->catdir( 
                $absolute_sequences_folder,
                miRkwood::Paths::get_basename_for_novel_miRNA()
            );
        }
    }

    my $size = length $candidate->{'sequence'};
    my $cfg    = miRkwood->CONFIG();
    my $candidate_name = $candidate->get_shortened_name();

    my $star = '<font color=\'#FFBE00\'>&#x2605;</font>';
    my $coche = '<font color=\'#41BE47\'>&#x2713;</font>';
    my $arrow = '<font color=\'#BDBDBD\'>&uarr;</font>';

    ### make files in 'sequences' folder
    my $rel_candidate_fasta_file = File::Spec->catfile( $relative_sequences_folder, "$candidate_name.fa" );
    my $abs_candidate_fasta_file = File::Spec->catfile( $absolute_sequences_folder, "$candidate_name.fa" );
    open( my $FASTA_FILE, '>', $abs_candidate_fasta_file ) 
      or die("Cannot open file $abs_candidate_fasta_file: $!");
    print {$FASTA_FILE} $candidate->candidateAsFasta()
      or die("Cannot write in file $abs_candidate_fasta_file: $!");
    close($FASTA_FILE)
      or die("Cannot close file $abs_candidate_fasta_file: $!");

    my $rel_vienna_file = File::Spec->catfile( $relative_sequences_folder, "$candidate_name.txt" );
    my $abs_vienna_file = File::Spec->catfile( $absolute_sequences_folder, "$candidate_name.txt" );
    open( my $VIENNA_FILE, '>', $abs_vienna_file )
      or die("Cannot open $abs_vienna_file: $!");
    print {$VIENNA_FILE} $candidate->candidateAsVienna(0)
      or die("Cannot write in file $abs_vienna_file: $!");
    close($VIENNA_FILE)
      or die("Cannot close file $abs_vienna_file: $!");

    my $rel_vienna_file_optimal = '';
    my $abs_vienna_file_optimal = '';
    if ( $candidate->{'structure_optimal'} ne $candidate->{'structure_stemloop'} ){
        $rel_vienna_file_optimal = File::Spec->catfile( $relative_sequences_folder, $candidate_name . '_optimal.txt' );
        $abs_vienna_file_optimal = File::Spec->catfile( $absolute_sequences_folder, $candidate_name . '_optimal.txt' );
        open( my $VIENNA_FILE_OPT, '>', $abs_vienna_file_optimal )
          or die("Cannot open $abs_vienna_file_optimal: $!");
        print {$VIENNA_FILE_OPT}
          $candidate->candidateAsVienna(1)
          or die("Cannot write in file $abs_vienna_file_optimal: $!");
        close($VIENNA_FILE_OPT)
          or die("Cannot close file $abs_vienna_file_optimal: $!");
    }

    my $rel_alternatives_file = '';
    my $abs_alternatives_file = '';
    if( $candidate->{'alternatives'} and scalar(keys%{ $candidate->{'alternatives'} }) ){
        $rel_alternatives_file = File::Spec->catfile( $relative_sequences_folder, $candidate_name . '_alternatives.txt' );
        $abs_alternatives_file = File::Spec->catfile( $absolute_sequences_folder, $candidate_name . '_alternatives.txt' );
        open( my $ALT_FILE, '>', $abs_alternatives_file )
          or die("Cannot open $abs_alternatives_file: $!");
        print {$ALT_FILE}
          $candidate->alternativeCandidatesAsVienna()
          or die("Cannot write in file $abs_alternatives_file: $!");
        close($ALT_FILE)
          or die("Cannot close file $abs_alternatives_file: $!");
    }


    ### links
    my $linkFasta         = "$rel_candidate_fasta_file";
    my $linkVienna        = "$rel_vienna_file";
    my $linkAlternatives  = "$rel_alternatives_file";
    my $linkViennaOptimal = "$rel_vienna_file_optimal";


    ### HTML sections
    my $mfei = "<i>MFEI</i> $candidate->{'mfei'}";
    if ( $candidate->{'mfei'} < -0.8 ){
        $mfei = "<font color='#E13EA5'><i>MFEI</i> $candidate->{'mfei'}</font>";
    }

    # Reads
    my $mirna_html = '';
    my $reads_html = '';
    if ( $pipeline_type eq 'smallRNAseq' ){
        # miRNA sequence
        my $mirna_sequence = '';
        if ( $candidate->{'mirna_length'} eq '' ||  $candidate->{'mirna_length'} eq '0' ){
            $mirna_sequence = 'none';
        }
        else{
            $mirna_sequence = "$candidate->{'mirna_sequence'} ($candidate->{'mirna_length'} nt)";
            $mirna_html = "<li><b>miRNA sequence:</b> $mirna_sequence</li>\n";
            if ( $mirna_type eq 'novel_miRNA' ){
                $mirna_html .= "<li><b>miRNA depth:</b> $candidate->{'mirna_depth'} (weight: $candidate->{'weight'})</li>\n";
                $mirna_html .= "<li><b>Other precursors with the same miRNA:</b>";
                if ( defined( $candidate->{'list_id_with_same_mirna'} ) && scalar( @{ $candidate->{'list_id_with_same_mirna'} } ) ){
                    foreach ( @{ $candidate->{'list_id_with_same_mirna'} } ){
                        if ( /(.*)__(\d+)-(\d+)/ ){
                            $mirna_html .= " <a href='#$1-$2-$3'>$1__$2-$3</a>";
                        }
                    }
                    $mirna_html .= "</li>\n";
                }
                else{
                    $mirna_html .= ' none</li>';
                }
            }
        }

        # Reads
        # this is not very robust. be cautious if you change the tree view
        my $reads_path = File::Spec->catdir( File::Spec->updir(), File::Spec->updir(), miRkwood::Paths::get_reads_dir_name(), $mirna_type);
        my $reads_file = File::Spec->catfile( $reads_path, $candidate->{'identifier'} . '.txt' );
        my $reads_score = '';
        my $absolute_read_cloud_path = File::Spec->catfile(
                miRkwood::Paths::get_dir_reads_path_from_job_dir( $cfg->param('job.directory') ),
                $mirna_type,
                $candidate->{'identifier'}.'.txt');
        my $read_cloud  = include_read_cloud_in_html( $absolute_read_cloud_path, $candidate->{'length'}, $candidate->{'nb_reads'} );
        my $read_duplex = '';

        my $nb_reads = $candidate->{'nb_reads'};
        if ( $candidate->{'criteria_nb_reads'} == 1 ){
            $nb_reads = "<font color='#12D0E5'>$candidate->{'nb_reads'}</font>";
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
        else{
            $reads_score .= '<li><b>Quality:</b> ';
            if ( $candidate->{'quality'} == 0 ){
                $reads_score .= 'random distribution</li>';
            }
            if ( $candidate->{'criteria_nb_reads'} == 1 && $candidate->{'criteria_reads_mirna'} == 0 ){
                $reads_score .= "the locus contains more than 10 reads $star</li>";
            }
            if ( $candidate->{'criteria_nb_reads'} == 0 && $candidate->{'criteria_reads_mirna'} == 1 ){
                $reads_score .= "more than half of the reads intersect either with the miRNA or the miRNA* $star</li>";
            }
            if ( $candidate->{'quality'} == 2 ){
                $reads_score .= "the locus contains more than 10 reads, and more than half of them intersect either with the miRNA or the miRNA* $star$star</li>";
            }
        }

        $reads_html = <<"END_TXT";
$read_duplex
<li>
    <b>Total number of reads mapped to the precursor:</b> $nb_reads [<a href='$reads_file'>download<a/>]
</li>
$reads_score
<pre style='font-size:0.8em;'>
$read_cloud</pre>
END_TXT
    }

    # miRBase alignments
    my $alignments_html = '';
    if ( ! defined( $candidate->{'mirbase_id'} ) && $cfg->param('options.align') ){
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

            $alignments_html .= '<br /><pre>' . $candidate->include_alignments_in_html( 'html' ) . '</pre>';
        }
        $alignments_html .= '</li>';
    }

    # Alternative sequences
    my $alternatives_HTML = '';
    if ( -e $abs_alternatives_file ){
        $alternatives_HTML = "[<a href='$linkAlternatives'>alternative sequences</a>]"
    }

    # Optimal structure
    my $optimal_HTML = '';
    if ( -e $abs_vienna_file_optimal ){
        $optimal_HTML = "[<a href='$linkViennaOptimal'>optimal MFE structure</a>]";
    }

    # miRBase name
    my $mirbase_name = '';
    if ( defined( $candidate->{'mirbase_id'} ) ){
        my $mirbase_link = miRkwood::Utils::make_mirbase_link( $candidate->{'mirbase_id'} );
        $mirbase_name = "<li><b>miRbase name:</b> <a href='$mirbase_link'>$candidate->{'identifier'}</a></li>";
    }

    ### make page
    my ($start, $end) = split( /-/, $candidate->{'position'});
    $start = miRkwood::Utils::make_numbers_more_readable( $start );
    $end = miRkwood::Utils::make_numbers_more_readable( $end );
    my $html = <<"END_TXT";
<h3 id='$candidate->{'name'}-$candidate->{'position'}'><a href='#table_$candidate->{'name'}-$candidate->{'position'}' class='nodecoration'>Results for $candidate->{'name'}__$start-$end ($candidate->{'strand'}) $arrow</a></h3>

    <ul>
        $mirbase_name
        <li>
            <b>Chromosome: </b>$candidate->{'name'}
        </li>
        <li>
            <b>Position:</b> $start-$end ($size nt)
        </li>
        <li>
            <b>Strand:</b> $candidate->{'strand'}
        </li>
        <li>
            <b>G+C content:</b> $candidate->{'%GC'} &#37;
        </li>
        $mirna_html
        <li>
            <b>miRNA precursor:</b> [<a href='$linkFasta'>FASTA sequence</a>] 
                                    [<a href='$linkVienna'>stem-loop structure</a>] 
                                    $optimal_HTML 
                                    $alternatives_HTML
        </li>
        <li>
            <b>Stability of the secondary structure of the precursor:</b> <i>MFE</i> $candidate->{'mfe'} kcal/mol | 
                                                                          <i>AMFE</i> $candidate->{'amfe'} | 
                                                                          $mfei
        </li>
        $reads_html
        $alignments_html
    </ul>
END_TXT


    return $html;
}

=method include_read_cloud_in_html

  Read contents of read cloud file (except header)
  and return it.

=cut
sub include_read_cloud_in_html {
    my @args = @_;
    my $read_cloud_file = shift @args;
    my $locus_length    = shift @args;
    my $locus_nb_reads  = shift @args;
    my $threshold = int( $locus_nb_reads / $locus_length ) + 1;
    my $result = '';
    open(my $IN, '<', $read_cloud_file) or return '';
    my $line = <$IN>;
    while ( <$IN> ){
        chomp;
        if ( $_ !~ /^Chromosome/ && $_ !~ /^Position/ && $_ !~ /^Strand/ && $_ ne '' ){
            my $depth = -1;
            if ( /depth=([\d]+)/ ){
                $depth = $1;
            }
            if ( $depth == -1 || $depth > $threshold ){
                $result .= "$_\n";
            }
        }
    }
    close $IN;
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
body {
    font-family:Sans-serif;
}
table{
    border:1px solid black;
    border-collapse:collapse;
    width:80%;
    color:#505050;
    cellpadding: 0px;
}
th, td {
    border:1px solid #505050;
    text-align:center;
}
span.mature {
    color: #F3791B;
}
ul {
   list-style-type: circle;
}
td .nodecoration
{ 
    text-decoration:none;
    color: #505050;
    display:block;
    width:100%; 
    height: 100%;
}
h2 .nodecoration
{ 
    text-decoration:none;
    color: black;
    display:block;
    width:100%; 
    height: 100%;
}
h3 .nodecoration
{
    text-decoration:none;
    color: black;
    display:block;
    width:100%;
    height: 100%;
}
.non_clickable_cell:hover{
    background-color: white;
}
.clickable_cell:hover{
    background-color: #EEEEEE;
    cursor:pointer;
}
.read_cell:hover{
    background-color: #E2F9FB;
    cursor:pointer;
}
.star_cell:hover{
   background-color: #FCF3C3;
   cursor:pointer;
}
.alignment_cell:hover{
   background-color: #E4F7E6;
   cursor:pointer;
}
END_TXT
    return ($css);
}

1;
