package miRkwood::CLI;

# ABSTRACT: Code for the command line interface

use strict;
use warnings;

use miRkwood;
use miRkwood::Results;
use miRkwood::ResultsExporterMaker;

=method process_results_dir_for_offline

Process the results in the given directory for offline use.

Usage:
  miRkwood::CLI::process_results_dir_for_offline($folder);

=cut

sub process_results_dir_for_offline {
    my @args           = @_;
    my $candidates_dir = shift @args;
    my $output_folder  = shift @args;

# In debug mode (without executing the pipeline), we need to set the config file
#miRkwood->CONFIG_FILE(miRkwood::Paths->get_job_config_path( $output_folder ));

    my $tmp_pieces_folder = File::Spec->catdir( $output_folder, 'pieces' );
    if ( !-e $tmp_pieces_folder ) {
        mkdir $tmp_pieces_folder or die("Error when creating $tmp_pieces_folder");
    }

    my %results = miRkwood::Results->deserialize_results($candidates_dir);

    my $html = make_html_from_results( \%results, $output_folder );

    my $html_page = File::Spec->catfile( $output_folder, 'results.html' );
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
    my $results = shift @args;
    my %results = %{$results};
    my $output_folder = shift @args;
    my $pieces_folder = File::Spec->catdir('pieces');
    
    my ($css) = get_page_css();
    my $page = '<h2>Overview of results</h2>';
    my $exporter = miRkwood::ResultsExporterMaker->make_html_results_exporter();
    $exporter->initialize('', \%results);
    $page .= $exporter->perform_export();

    $page .= make_all_exports( \%results, $output_folder );
    while ( my ( $key, $value ) = each %results ) {
        my $candidate_html =
          make_candidate_page( $value, $pieces_folder, $output_folder );
        $page .= $candidate_html;
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
    my $results_ref   = shift @args;
    my $output_folder = shift @args;
    my $pieces_folder = File::Spec->catdir('pieces');
    my $id_job = '';

    my $exporter = miRkwood::ResultsExporterMaker->make_fasta_results_exporter();
    $exporter->initialize($id_job, $results_ref);
    $exporter->export_on_disk( $output_folder );
    my $fasta_file = File::Spec->catfile($output_folder, $exporter->get_filename());

    $exporter = miRkwood::ResultsExporterMaker->make_dotbracket_results_exporter();
    $exporter->initialize($id_job, $results_ref);
    $exporter->export_on_disk( $output_folder );
    my $dotbracket_file = File::Spec->catfile($output_folder, $exporter->get_filename());

    $exporter = miRkwood::ResultsExporterMaker->make_gff_results_exporter();
    $exporter->initialize($id_job, $results_ref);
    $exporter->export_on_disk( $output_folder );
    my $gff_file = File::Spec->catfile($output_folder, $exporter->get_filename());

    $exporter = miRkwood::ResultsExporterMaker->make_csv_results_exporter();
    $exporter->initialize($id_job, $results_ref);
    $exporter->export_on_disk( $output_folder );
    my $csv_file = File::Spec->catfile($output_folder, $exporter->get_filename());

    my $html = '<h3>Get results as</h3> <ul>';
    $html .= "<li><a href='$csv_file'>tab-delimited format (csv)</a></li>";
    $html .= "<li><a href='$fasta_file'>Fasta</a></li>";
    $html .=
"<li><a href='$dotbracket_file'>dot-bracket format (plain sequence + secondary structure)</a></li>";
    $html .= "<li><a href='$gff_file'>gff format</a></li>";
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
    my $output_folder = shift @args;

    my $size = length $candidate->{'sequence'};
    my $cfg    = miRkwood->CONFIG();
    my $candidate_name = $candidate->get_shortened_name();

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
        $candidate_name . '__alternatives.txt' );
    open( my $ALT_FILE, '>',
        File::Spec->catfile( $output_folder, $alternatives_file ) )
      or die("Cannot open $alternatives_file: $!");
    print {$ALT_FILE}
      $candidate->alternativeCandidatesAsVienna()
      or die("Cannot write in file $alternatives_file: $!");
    close($ALT_FILE)
      or die("Cannot close file $alternatives_file: $!");
      
    my $reads_file = File::Spec->catfile( 
        $output_folder ."/clusters/", $candidate->{'identifier'} . '.txt' );

    my $linkFasta         = "$candidate_fasta_file";
    my $linkVienna        = "$vienna_file";
    my $linkAlternatives  = "$alternatives_file";
    my $linkViennaOptimal = "$vienna_file_optimal";
    my $linkReads         = "$reads_file";  

    my $Vienna_HTML =
"<ul><li><b>Stem-loop structure (dot-bracket format):</b> <a href='$linkVienna'>download</a>";
    if ( $candidate->{'structure_stemloop'} ne $candidate->{'structure_optimal'} ) {
        $Vienna_HTML .=
"</li><li><b>Optimal MFE secondary structure (dot-bracket format):</b> <a href='$linkViennaOptimal'>download</a></li></ul>";
    }
    else {
        $Vienna_HTML .=
          '<br/>This stem-loop structure is the MFE structure.</li></ul>';
    }
    my $alternatives_HTML =
      '<b>Alternative candidates (dot-bracket format):</b> ';
    if ( $candidate->{'alternatives'} ) {
        $alternatives_HTML .= "<a href='$linkAlternatives'>download</a>";
    }
    else {
        $alternatives_HTML .= '<i>None</i>';
    }

    my $alignmentHTML = q{};
    if ( $candidate->{'alignment'} ) {
        $alignmentHTML =
          $candidate->make_alignments_HTML();
    }
    else {
        $alignmentHTML = 'No alignment has been found.';
    }

    my $html = <<"END_TXT";
<h2 id='$candidate->{'name'}-$candidate->{'position'}'>Results for $candidate->{'name'}, $candidate->{'position'}</h2>
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
    <li>
      <b>Sequence (FASTA format):</b> <a href='$linkFasta'>download</a>
    </li>
    <li>
      $alternatives_HTML
    </li>
    <li>
      <b>Reads:</b> <a href='$linkReads'>download<a/>
    </li>
    </ul>
<h3>Secondary structure</h3>
    $Vienna_HTML
<h3>Thermodynamics stability</h3>
    <ul>
    <li>
      <b>MFE:</b> $candidate->{'mfe'} kcal/mol
    </li>
    <li>
      <b>AMFE:</b> $candidate->{'amfe'}
    </li>
    <li>
      <b>MFEI:</b> $candidate->{'mfei'}
    </li>
    </ul>
END_TXT

    if ( $cfg->param('options.align') ){   
        
        $html .=  <<"END_TXT";
<h3>miRBase alignments</h3>

$alignmentHTML
END_TXT

    }

    return $html;
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
END_TXT
    return ($css);
}

1;
