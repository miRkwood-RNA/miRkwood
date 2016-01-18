#!/usr/bin/perl
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Spec;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood;
use miRkwood::Paths;
use miRkwood::Utils;
use miRkwood::Pipeline;
use miRkwood::WebTemplate;
use miRkwood::BEDHandler;
use miRkwood::Results;
use miRkwood::ResultsExporterMaker;


##### Page settings
my $header_menu = miRkwood::WebTemplate::get_header_menu( 'smallrnaseq' );
my $footer      = miRkwood::WebTemplate::get_footer();
my $web_scripts = miRkwood::WebPaths->get_web_scripts();

my @css = (
	miRkwood::WebTemplate->get_server_css_file(),
	miRkwood::WebTemplate->get_css_file(),
    miRkwood::WebTemplate->get_mirkwood_css_file()
);
my @js = (
	File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'results.js' ),
	File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'graphics.js' ),
	File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'miARN.js' ),
	File::Spec->catfile( miRkwood::WebPaths->get_js_path(), 'jquery.min.js' ),
	miRkwood::WebTemplate->get_bioinfo_js_file()
);

my $help_page = File::Spec->catfile( File::Spec->catdir( miRkwood::WebPaths->get_html_path(), 'smallRNAseq'), 'help.php');

##### Parameters
my $cgi = CGI->new();
my $id_job = $cgi->param('run_id');    # get id job



##### Create page
my $valid = miRkwood::Results->is_valid_jobID($id_job);
my $html = '';
my $HTML_additional = '';
my $HTML_results = '';
my $page = '';


$HTML_additional .= '<p class="header-results" id="job_id"><b>Job ID:</b> ' . $id_job . '</p>';

if ( $valid ){

    my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($id_job);

    my $run_options_file = miRkwood::Paths->get_job_config_path($absolute_job_dir);
    miRkwood->CONFIG_FILE($run_options_file);
    my $cfg = miRkwood->CONFIG();

    my $basic_known_yaml = File::Spec->catfile( $absolute_job_dir, 'basic_known_candidates.yml');
    my $basic_yaml = File::Spec->catfile( $absolute_job_dir, 'basic_candidates.yml');

    my $initial_bed         = miRkwood::Paths::get_bed_file ( $absolute_job_dir, '', 'bed' );

    my $bed_sizes;
    my $nb_new_results                   = 0;
    my $nb_known_results                 = 0;
    my $nb_reads_known_miRNAs            = 0;
    my $nb_reads_new_miRNAs              = 0;
    my $nb_orphan_reads                  = 0;
    my $percentage_CDS_reads             = 0;
    my $percentage_other_reads           = 0;
    my $percentage_multi_reads           = 0;
    my $percentage_known_miRNAs_reads    = 0;
    my $percentage_new_miRNAs_reads      = 0;
    my $percentage_orphan_clusters_reads = 0;
    my $percentage_orphan_hairpins_reads = 0;

    if ( $cfg->param('job.title') ) {
        $HTML_additional .= "<p class='header-results' id='job_title'><b>Job title:</b> " . $cfg->param('job.title') . '</p>';
    }

	unless ( miRkwood::Results->is_job_finished($id_job) ) {
		$HTML_additional .= "<p class='warning'>Still processing...</p>";
	} else {
        ##### Summary of options
        my $basename_bed = '';
        if ( $initial_bed =~ /.*[\/\\]([^\/]+)[.]bed/ ){
            $basename_bed = $1;
        }
        $HTML_additional .= "<div class='results_summary'><ul>";
        $HTML_additional .= '<h2>Options summary</h2>';
        $HTML_additional .= '<br />';
        $HTML_additional .= "<li><em>BED file:</em> $basename_bed.bed</li>";

        # Reference species
        if ( $cfg->param('job.plant') ){
            $HTML_additional .= '<li><em>Reference species:</em> ' . $cfg->param('job.plant') . '</li>';
        }

        # Align
        if ( $cfg->param('options.align') ){
            $HTML_additional .= '<li><em>Flag conserved mature miRNAs:</em> Yes</li>';
        }
        else{
            $HTML_additional .= '<li><em>Flag conserved mature miRNAs:</em> No</li>';
        }

        # MFEI
        if ( $cfg->param('options.mfei') ){
            $HTML_additional .= '<li><em>Select only sequences with MFEI < -0.6:</em> Yes</li>';
        }
        else{
            $HTML_additional .= '<li><em>Select only sequences with MFEI < -0.6:</em> No</li>';
        }

        # Ranfold
        if ( $cfg->param('options.randfold') ){
            $HTML_additional .= '<li><em>Compute thermodynamic stability:</em> Yes</li>';
        }
        else{
            $HTML_additional .= '<li><em>Compute thermodynamic stability:</em> No</li>';
        }

        # CDS
        if ( $cfg->param('options.filter_CDS') ){
            $HTML_additional .= '<li><em>Filter CoDing Sequences:</em> Yes</li>';
        }
        else{
            $HTML_additional .= '<li><em>Filter CoDing Sequences:</em> No</li>';
        }

        # tRNA/rRNA/snoRNA
        if ( $cfg->param('options.filter_tRNA_rRNA') ){
            $HTML_additional .= '<li><em>Filter tRNA/rRNA/snoRNA:</em> Yes</li>';
        }
        else{
            $HTML_additional .= '<li><em>Filter tRNA/rRNA/snoRNA:</em> No</li>';
        }

        # Multimapped reads
        if ( $cfg->param('options.multimapped_interval') ne '[0;0]' ){
            $HTML_additional .= '<li><em>Filter multiply mapped reads:</em> Yes</li>';
        }
        else{
            $HTML_additional .= '<li><em>Filter multiply mapped reads:</em> No</li>';
        }

        # Orphan hairpins
        if ( $cfg->param('options.filter_bad_hairpins') ){
            $HTML_additional .= '<li><em>Filter low quality hairpins:</em> Yes</li>';
        }
        else{
            $HTML_additional .= '<li><em>Filter low quality hairpins:</em> No</li>';
        }

        $HTML_additional .= '</ul></div>';


        ##### Summary of results

        ### Count the number of results
        $nb_new_results   = miRkwood::Results->number_of_results_bis( $id_job, 'novel_miRNA' );
        $nb_known_results = miRkwood::Results->number_of_results_bis( $id_job, 'known_miRNA' );

        ### Count the number and percentage of reads in each category
        my $bed_sizes_file = File::Spec->catfile( $absolute_job_dir, miRkwood::Paths::get_bed_size_file_name() );
        my $species = $cfg->param('job.plant');
        if ( -e $bed_sizes_file ){
            open (my $FH, '<', $bed_sizes_file) or die "ERROR while opening $bed_sizes_file: $!";
            while ( <$FH> ){
                if ( $_ !~ /^#/ ){
                    chomp;
                    my @line = split(/\t/);
                    my $name = '';
                    if ( $line[0] eq "$basename_bed.bed" ){
                        $name = $basename_bed;
                    }
                    elsif ( $line[0] =~ /${species}_CDS[.]bed/ ){
                        $name = 'CDS';
                    }
                    elsif ( $line[0] =~ /${species}_tRNA_rRNA_snoRNA[.]bed/ ){
                        $name = 'tRNA_rRNA_snoRNA';
                    }
                    elsif ( $line[0] =~ /${basename_bed}_(.*)[.]bed/ ){
                        $name = $1;
                    }
                    $bed_sizes->{$name}{'reads'} = $line[1];
                    $bed_sizes->{$name}{'unique_reads'} = $line[2];
                }
            }
            close $FH;
        }

        if ( $bed_sizes->{$basename_bed}{'reads'} == 0 ){
            $bed_sizes->{$basename_bed}{'reads'} = 1;
        }

        $percentage_CDS_reads = $bed_sizes->{'CDS'}{'reads'} / $bed_sizes->{$basename_bed}{'reads'} * 100;
        $percentage_other_reads = $bed_sizes->{'tRNA_rRNA_snoRNA'}{'reads'} / $bed_sizes->{$basename_bed}{'reads'} * 100;
        $percentage_multi_reads = $bed_sizes->{'multimapped'}{'reads'} / $bed_sizes->{$basename_bed}{'reads'} * 100;
        $percentage_orphan_clusters_reads = $bed_sizes->{'orphan_clusters'}{'reads'} / $bed_sizes->{$basename_bed}{'reads'} * 100;
        $percentage_orphan_hairpins_reads = $bed_sizes->{'orphan_hairpins'}{'reads'} / $bed_sizes->{$basename_bed}{'reads'} * 100;

        # Known miRNAs
        $nb_reads_known_miRNAs = miRkwood::Results::count_reads_in_basic_yaml_file( $basic_known_yaml );
        $percentage_known_miRNAs_reads = $nb_reads_known_miRNAs / $bed_sizes->{$basename_bed}{'reads'} * 100;

        # New miRNAs
        $nb_reads_new_miRNAs = miRkwood::Results::count_reads_in_basic_yaml_file( $basic_yaml );
        $percentage_new_miRNAs_reads = $nb_reads_new_miRNAs / $bed_sizes->{$basename_bed}{'reads'} * 100;

        # Orphan reads
        $nb_orphan_reads = $bed_sizes->{$basename_bed}{'reads'}
                            - $nb_reads_known_miRNAs
                            - $nb_reads_new_miRNAs
                            - $bed_sizes->{'CDS'}{'reads'}
                            - $bed_sizes->{'tRNA_rRNA_snoRNA'}{'reads'}
                            - $bed_sizes->{'multimapped'}{'reads'}
                            - $bed_sizes->{'orphan_clusters'}{'reads'}
                            - $bed_sizes->{'orphan_hairpins'}{'reads'};

        ### Create HTML
        my $arguments = '?jobID=' . $id_job;
        my $known_url = miRkwood::WebTemplate::get_cgi_url('BAMresults_for_mirnas.pl') . $arguments . '&type=known_miRNA';
        my $new_url = miRkwood::WebTemplate::get_cgi_url('BAMresults_for_mirnas.pl') . $arguments . '&type=novel_miRNA';
        my $exportFileLink = miRkwood::WebTemplate::get_cgi_url('getBEDFile.pl') . '?jobId=' . $id_job;

        # Create reads barchart
        my $total_width = 635;
        my $barchart = miRkwood::Results->make_reads_barchart( $total_width,
                                                               $percentage_CDS_reads,
                                                               $percentage_other_reads,
                                                               $percentage_multi_reads,
                                                               $percentage_known_miRNAs_reads,
                                                               $percentage_new_miRNAs_reads,
                                                               $percentage_orphan_clusters_reads,
                                                               $percentage_orphan_hairpins_reads, );

        my $nb_total_reads = miRkwood::Utils::make_numbers_more_readable( $bed_sizes->{$basename_bed}{'reads'} );
        my $nb_total_unq_reads = miRkwood::Utils::make_numbers_more_readable( $bed_sizes->{$basename_bed}{'unique_reads'} );
        $HTML_results .= "<div class='results_summary'><ul>";
        $HTML_results .= '<h2>Results summary</h2>';
        $HTML_results .= '<br />';
        $HTML_results .= "<li><em>Total number of reads:</em> $nb_total_reads ($nb_total_unq_reads unique reads)</li>";
        $HTML_results .= '<br />';
        $HTML_results .= "$barchart<br />";


        if ( $bed_sizes->{'CDS'}{'reads'} > 0 ){
            my $nb_reads_CDS = miRkwood::Utils::make_numbers_more_readable( $bed_sizes->{'CDS'}{'reads'} );
            $HTML_results .= "<li id='li_CDS'><span id='normal'><em>CoDing Sequences:</em> $nb_reads_CDS reads (<a href='$exportFileLink&type=_CDS'>download</a>)</span></li>";
        }
        else {
            $HTML_results .= "<li id='li_CDS'><span id='normal'><em>CoDing Sequences:</em> 0 reads</span></li>";
        }

        if ( $bed_sizes->{'tRNA_rRNA_snoRNA'}{'reads'} > 0 ){
            my $nb_reads_tRNA_rRNA_snoRNA = miRkwood::Utils::make_numbers_more_readable( $bed_sizes->{'tRNA_rRNA_snoRNA'}{'reads'} );
            $HTML_results .= "<li id='li_other'><span id='normal'><em>tRNA/rRNA/snoRNA:</em> $nb_reads_tRNA_rRNA_snoRNA reads (<a href='$exportFileLink&type=_tRNA_rRNA_snoRNA'>download</a>)</span></li>";
        }
        else {
            $HTML_results .= "<li id='li_other'><span id='normal'><em>tRNA/rRNA/snoRNA:</em> 0 reads</span></li>";
        }

        if ( $bed_sizes->{'multimapped'}{'reads'} > 0 ){
            my $nb_reads_multimapped = miRkwood::Utils::make_numbers_more_readable( $bed_sizes->{'multimapped'}{'reads'} );
            $HTML_results .= "<li id='li_multimapped'><span id='normal'><em>Multiply mapped reads:</em> $nb_reads_multimapped reads (<a href='$exportFileLink&type=_multimapped'>download</a>)</span></li>";
        }
        else {
            $HTML_results .= "<li id='li_multimapped'><span id='normal'><em>Multiply mapped reads:</em> 0 reads</span></li>";
        }

        if ( $bed_sizes->{'orphan_clusters'}{'reads'} > 0 ){
            my $nb_reads_orphan_clusters = miRkwood::Utils::make_numbers_more_readable( $bed_sizes->{'orphan_clusters'}{'reads'} );
            $HTML_results .= "<li id='li_orphan_clusters'><span id='normal'><em>Orphan clusters of reads:</em> $nb_reads_orphan_clusters reads (<a href='$exportFileLink&type=_orphan_clusters'>download</a>)</span></li>";
        }
        else {
            $HTML_results .= "<li id='li_orphan_clusters'><span id='normal'><em>Orphan clusters of reads:</em> 0 reads</span></li>";
        }

        if ( $bed_sizes->{'orphan_hairpins'}{'reads'} > 0 ){
            my $nb_reads_orphan_hairpins = miRkwood::Utils::make_numbers_more_readable( $bed_sizes->{'orphan_hairpins'}{'reads'} );
            $HTML_results .= "<li id='li_orphan_hairpins'><span id='normal'><em>Orphan hairpins:</em> $nb_reads_orphan_hairpins reads (<a href='$exportFileLink&type=_orphan_hairpins'>download</a>)</span></li>";
        }
        else {
            $HTML_results .= "<li id='li_orphan_hairpins'><span id='normal'><em>Orphan hairpins:</em> 0 reads</span></li>";
        }

        $nb_orphan_reads = miRkwood::Utils::make_numbers_more_readable( $nb_orphan_reads );
        $HTML_results .= "<li id='li_unclassified_reads'><span id='normal'><em>Unclassified reads:</em> $nb_orphan_reads reads</span></li>";

        $nb_reads_known_miRNAs = miRkwood::Utils::make_numbers_more_readable( $nb_reads_known_miRNAs );
        if ( $nb_known_results > 0 ){
            $nb_known_results = miRkwood::Utils::make_numbers_more_readable( $nb_known_results );
            $HTML_results .= "<li id='li_known_miRNAs'><span id='normal'><em>Known miRNAs:</em> $nb_known_results sequence(s) - $nb_reads_known_miRNAs reads (<a href=$known_url>see results</a>)</span></li>";
        }
        else {
            $HTML_results .= "<li id='li_known_miRNAs'><span id='normal'><em>Known miRNAs:</em> $nb_known_results sequence(s) - $nb_reads_known_miRNAs reads</span></li>";
        }

        $nb_reads_new_miRNAs = miRkwood::Utils::make_numbers_more_readable( $nb_reads_new_miRNAs );
        if ( $nb_new_results > 0 ){
            $nb_new_results = miRkwood::Utils::make_numbers_more_readable( $nb_new_results );
            $HTML_results .= "<li id='li_new_miRNAs'><span id='normal'><em>Novel miRNAs:</em> $nb_new_results sequence(s) - $nb_reads_new_miRNAs reads (<a href=$new_url>see results</a>)</span></li>";
        }
        else {
            $HTML_results .= "<li id='li_new_miRNAs'><span id='normal'><em>Novel miRNAs:</em> $nb_new_results sequence(s) - $nb_reads_new_miRNAs reads</span></li>";
        }
        $HTML_results .= "</ul>";
        $HTML_results .= '<p style="margin-left: 625px;"><a href="' . $help_page . '#overview">help</a></p>';
        $HTML_results .= "</div>";

    }

    my $arborescence = "<a href='/'>home</a> &gt; ";
    $arborescence .= "<a href='/theme_page/rna.html'>RNA</a> &gt; ";
    $arborescence .= "<a href='/mirkwood/mirkwood.php'>mirkwood</a>";

    $page = <<"END_TXT";
<body onload="main('all');">
    <div class="frametitle">
        <h1 id="title">miRkwood small RNA-seq</h1>
    </div>

    <div id="center_sup">
        <div id="link_home" style="display:inline-block"><a href="/mirkwood/smallRNAseq/index.php" class="text_onglet"><img src="/Style/icon/home_w.png" alt="home_general"/></a></div>
        <div class="tabs" id="menu_central" style="display:inline-block"> 
            $header_menu
        </div>
        <div id="arborescence">$arborescence</div>
    </div>

    <div id="main">
        $HTML_additional

        $HTML_results

    </div><!-- main -->

    $footer
</body>
END_TXT

    my $title = 'miRkwood small RNA-seq - Summary of results';
    $html = miRkwood::WebTemplate::get_HTML_page_for_body($page, \@css, \@js, $title);

}   # end if valid
else{   # job id is not a valid ID
    $HTML_additional .= '</div>';
	$page = <<"END_TXT";
<div id="main">
    $HTML_additional
    <p>No results available for the given job identifier $id_job.</p>
</div><!-- main -->
END_TXT

    my $title = 'miRkwood - No results';
	$html = miRkwood::WebTemplate::get_HTML_page_for_content( 'smallrnaseq', $page, \@css, \@js, $title );
}


print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
