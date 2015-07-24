#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood;
use miRkwood::WebPaths;
use miRkwood::Utils;
use miRkwood::Results;
use miRkwood::Candidate;
use miRkwood::WebTemplate;

my $cgi            = CGI->new();
my $jobId          = $cgi->param('jobID');
my $candidate_id   = $cgi->param('id');

my @css = (File::Spec->catfile(miRkwood::WebPaths->get_css_path(), 'results.css'));
my @js  = (miRkwood::WebTemplate->get_js_file());

my $job = miRkwood::Results->jobId_to_jobPath($jobId);

my $candidate;
my $html_contents;

my $cfg_path = miRkwood::Paths->get_job_config_path($job);
miRkwood->CONFIG_FILE($cfg_path);

my $cfg = miRkwood->CONFIG();

my $returnlink = '';
my $return_html = '';
my $help_page = '';

if (! eval {$candidate = miRkwood::CandidateHandler->retrieve_candidate_information($job, $candidate_id);}) {
    # Catching exception
    $html_contents = 'No results for the given identifiers';
}else{

    if ( $cfg->param('job.mode') eq 'WebBAM' ){
        if ( defined($candidate->{'precursor_name'}) ){
            $returnlink = miRkwood::WebTemplate::get_link_back_to_BED_known_results($jobId);
        }
        else{
            $returnlink = miRkwood::WebTemplate::get_link_back_to_BED_new_results($jobId);
        }
        $help_page = File::Spec->catfile( File::Spec->catdir( miRkwood::WebPaths->get_html_path(), 'smallRNAseq'), 'help.php');
    }
    else {
        $returnlink = miRkwood::WebTemplate::get_link_back_to_results($jobId);
        $help_page = File::Spec->catfile( File::Spec->catdir( miRkwood::WebPaths->get_html_path(), 'abinitio'), 'help.php');
    }
    $return_html = "<a class='returnlink' href='$returnlink'>Back to main results page</a>";

    my $image_url = $candidate->get_relative_image();

    my $size = length $candidate->{'sequence'};

    my $export_link = "./getCandidate.pl?jobId=$jobId&id=$candidate_id";

    my $linkFasta = "$export_link&type=fas";
    my $linkVienna = "$export_link&type=dot";
    my $linkAlternatives = "$export_link&type=alt";
    my $linkViennaOptimal = $linkVienna . '&optimal=1';
    my $linkReadsCloud = "$export_link&type=reads";
    my $htmlReadsCloud = '';
    my $reads_length_diagramm = $candidate->create_reads_length_diagramm();
    my $mirna = '';
    my $quality = '';
    if ( $cfg->param('job.mode') ne 'fasta' ){
        my $nb_reads = 0;
        foreach my $key (keys( %{$candidate->{'reads'}} )){
            $nb_reads += $candidate->{'reads'}{$key};
        }
        $htmlReadsCloud = <<"END_TXT";
        <h2>Reads</h2>
        <ul>
            <li><b>Number of reads:</b> $nb_reads</li>
            <li><b>Reads length distribution:</b> <br />
                $reads_length_diagramm
            </li>
            <li><b>Read cloud : </b><a href='$linkReadsCloud'>download</a>
        </ul>
END_TXT

        if ( $candidate->{'mirna_sequence'} ne '' ){
            $mirna = <<"END_TXT";
            <li>
                <b>miRNA:</b> $candidate->{'mirna_sequence'}
            </li>
            <li>
                <b>miRNA length:</b> $candidate->{'mirna_length'} nt
            </li>
END_TXT
        }

        $quality = <<"END_TXT";
        <h2>Quality <a href="$help_page#quality">[?]</a> </h2>
        <b>Quality:</b> $candidate->{'quality'}   <br /><ul>

END_TXT

        if ( $candidate->{'mirna_sequence'} ne '' ){
            $quality .= "<li><b>Existence of a miRNA:</b> Yes</li>";
        }
        else{
            $quality .= "<li><b>Existence of a miRNA:</b> No</li>";
        }
        if ( $candidate->{'mfei'} < -0.8 ){
            $quality .= "<li><b>MFEI < -0.8:</b> Yes</li>";
        }
        else{
            $quality .= "<li><b>MFEI < -0.8:</b> No</li>";
        }
        if ( $candidate->{'criteria_nb_reads'} ){
            $quality .= "<li><b>Criteria number of reads:</b> Yes</li>";
        }
        else{
            $quality .= "<li><b>Criteria number of reads:</b> No</li>";
        }
        if ( $candidate->{'criteria_star'} ){
            $quality .= "<li><b>Criteria presence of the miRNA:miRNA* duplex:</b> Yes</li>";
        }
        else{
            $quality .= "<li><b>Criteria presence of the miRNA:miRNA* duplex:</b> No</li>";
        }
        if ( $candidate->{'criteria_reads_mirna'} ){
            $quality .= "<li><b>Criteria precision of the precursor processing:</b> Yes</li>";
        }
        else{
            $quality .= "<li><b>Criteria precision of the precursor processing:</b> No</li>";
        }

        $quality .= "</ul>";

    }

    my $Vienna_HTML = "<li><b>Stem-loop structure (dot-bracket format):</b> <a href='$linkVienna'>download</a>";
    if( defined($candidate->{'structure_optimal'}) && ( $candidate->{'structure_stemloop'} ne $candidate->{'structure_optimal'}) ){
        $Vienna_HTML .= "</li><li><b>Optimal MFE secondary structure (dot-bracket format):</b> <a href='$linkViennaOptimal'>download</a></li>"
    } else {
        $Vienna_HTML .= '<br/>This stem-loop structure is the MFE structure.</li>'
    }

    my $alternatives_HTML = '<b>Alternative candidates (dot-bracket format):</b> ';
    if( $candidate->{'alternatives'} and scalar(keys%{ $candidate->{'alternatives'} }) ){
        $alternatives_HTML .= "<a href='$linkAlternatives'>download</a>"
    } else {
        $alternatives_HTML .= '<i>None</i>'
    }

    my $alignmentHTML;

    if ( !$cfg->param('options.align') || defined($candidate->{'precursor_name'}) ) {
        $alignmentHTML = qw{};
    }
    else {
        $alignmentHTML = '<h2>Conserved mature miRNA</h2>';

        if ( $candidate->{'alignment'} ) {
            $alignmentHTML .=
              $candidate->make_alignments_HTML();
        }
        else {
            $alignmentHTML .= 'No alignment has been found.';
        }
    }

    my $imgHTML = '';
    if ( $cfg->param('options.varna') ) {
        $imgHTML = "<img class='structure' id='structure' src='$image_url' height='300px' alt='$candidate->{'name'} secondary structure'>"
    }

    my $shufflesHTML = '';
    if ( $cfg->param('options.randfold') ) {
        $shufflesHTML = "<li>
          <b>Shuffles:</b> $candidate->{'shuffles'}
        </li>"
    }

    my $name = '';
    if ( defined($candidate->{'precursor_name'}) ){
        my $mirbase_link = miRkwood::Utils::make_mirbase_link( $candidate->{'identifier'} );
        $name = "<li><b>miRbase name :</b> <a href='$mirbase_link'>$candidate->{'precursor_name'}</a></li>";
    }

    $html_contents = <<"END_TXT";
            <div id = 'showInfo'>
        <ul>
        $name
        <li>
          <b>Chromosome: </b>$candidate->{'name'}
        </li>
        <li>
          <b>Position:</b> $candidate->{'position'} ($size nt)
        </li>
        <li>
          <b>Strand:</b> $candidate->{'strand'}
        </li>
        <li>
          <b>G+C content:</b> $candidate->{'%GC'}%
        </li>
        $mirna
        <li>
          <b>Sequence (FASTA format):</b> <a href='$linkFasta'>download</a>
        </li>
        $Vienna_HTML
        <li>
          $alternatives_HTML
        </li>
        </ul>
        $imgHTML
        $quality
        $htmlReadsCloud
        <h2>Thermodynamics stability</h2>
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
        $shufflesHTML
        </ul>
        $alignmentHTML

    </div><!-- showInfo -->
END_TXT
}

my $body  = <<"END_TXT";
    <body>
       <h1>Results for $candidate->{'name'}, $candidate->{'position'} ($candidate->{'strand'})</h1>
        $return_html
        $html_contents
        $return_html
    </body>
END_TXT

my $html = miRkwood::WebTemplate::get_HTML_page_for_body($body, \@css, \@js);

print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
###End###
