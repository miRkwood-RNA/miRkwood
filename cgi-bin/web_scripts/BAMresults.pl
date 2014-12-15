#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Spec;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood;
use miRkwood::Pipeline;
use miRkwood::WebTemplate;
use miRkwood::Results;


##### Page settings
my $header_menu = miRkwood::WebTemplate::get_header_menu();
my $footer      = miRkwood::WebTemplate::get_footer();

my @css = (miRkwood::WebTemplate->get_server_css_file(), miRkwood::WebTemplate->get_css_file());
my @js  = (miRkwood::WebTemplate->get_js_file());


##### Parameters
my $cgi = CGI->new;
my $id_job = $cgi->param('run_id');    # get id job

my $genome_file;
my $mirbase_file;
#~ my $cfg = miRkwood->CONFIG();



##### Create page
my $valid = miRkwood::Results->is_valid_jobID($id_job);
my $html = '';
my $HTML_additional = '';
my $page = '';

$HTML_additional .= '<p class="header-results" id="job_id"><b>Job ID:</b> ' . $id_job . '</p>';

if ( $valid ){

    my $absolute_job_dir = miRkwood::Results->jobId_to_jobPath($id_job);
    #~ my $cfg = File::Spec->catfile( $absolute_job_dir, 'run_options.cfg' );
    #~ my $cfg = miRkwood->CONFIG();
    #~ 
    #~ my $species = $cfg->param('job.plant');
    my $species = 'Arabidopsis_thaliana';
    if ( $species ne '' ){
        $genome_file = File::Spec->catfile( miRkwood::Paths->get_data_path(), "genomes/$species.fa");
        $mirbase_file = File::Spec->catfile( miRkwood::Paths->get_data_path(), "miRBase/${species}_miRBase.gff3");
    }
    else{
        $genome_file = "$absolute_job_dir/genome.fa";
    }    

    my $known_mirnas = '';
    my $mirna_bed    = '';
    my $other_bed    = '';
    my $final_bed    = '';
    
    opendir (my $dh, $absolute_job_dir) or die "Cannot open $absolute_job_dir : $!";
    while (readdir $dh) {
        if ( /_miRNAs.bed/ ){
            $mirna_bed = File::Spec->catfile($absolute_job_dir, $_);
        }
        elsif ( /_other.bed/ ){
            $other_bed = File::Spec->catfile($absolute_job_dir, $_);
        }
        elsif ( /_filtered.bed/ ){
            $final_bed = File::Spec->catfile($absolute_job_dir, $_);
        }
    }
    closedir $dh;
    
    if ( $mirna_bed ne '' ){
        $known_mirnas  = '<div id="table">';
        $known_mirnas .= miRkwood::Results->known_mirnas_for_jobID($mirna_bed, $mirbase_file, $genome_file, "$absolute_job_dir/reads/known" );
        $known_mirnas .= '</div>';
    }


    $page = <<"END_TXT";
<body>
    <div class="theme-border"></div>
    <div class="logo"></div>
    <div class="bloc_droit">
        $header_menu
        <div class="main main-full">
            $HTML_additional
            <br />
            <h2>Known miRNAs :</h2>
            $known_mirnas
        </div><!-- main -->
    </div><!-- bloc droit--> 
    $footer  
</body>    
    
END_TXT

    $html = miRkwood::WebTemplate::get_HTML_page_for_body($page, \@css, \@js);
}
else{
	$page = <<"END_TXT";
<div class="main">
    $HTML_additional
    <p>No results available for the given job identifier $id_job.</p>
</div><!-- main -->
END_TXT

	$html = miRkwood::WebTemplate::get_HTML_page_for_content( $page, \@css, \@js, 1 );
}




print <<"DATA" or die("Error when displaying HTML: $!");
Content-type: text/html

$html
DATA
