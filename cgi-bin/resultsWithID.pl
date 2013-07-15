#!/usr/bin/perl -w
use strict;
use warnings;

use Class::Struct;
use CGI;
my $cgi = CGI->new;
use CGI::Carp qw(fatalsToBrowser);
use Cwd qw( abs_path );
use File::Basename qw(dirname);
use File::Spec;
use Data::Dumper;
use FindBin;                     # locate this script
use lib "$FindBin::Bin/../lib";  # use the parent directory
use PipelineMiRNA::WebFunctions;


my $id_job = $cgi->param('run_id'); # récupération id job
my $name_job = $cgi->param('nameJob'); # récupération id job

my $HTML_header = <<'END_TXT';
Content-type: application/xhtml+xml

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

    <head>
        <link rel="stylesheet" type="text/css" href="/arn/css/script.css" />
        <script type="text/javascript" language="Javascript" src="/arn/js/results.js"> </script>
        <script type="text/javascript" src="/arn/js/graphics.js"></script>
        <script type="text/javascript" src="/arn/js/miARN.js"></script>
    </head>
END_TXT

my $HTML_additional = "";
if ($name_job ne "")
{
    $HTML_additional .= "<div class='titleJob' ><li>Title Job : ".$name_job."</li></div>";
}

my $valid = PipelineMiRNA::WebFunctions->is_valid_jobID($id_job);

if($valid){
    my %myResults = PipelineMiRNA::WebFunctions->get_structure_for_jobID($id_job);
    my $HTML_results = PipelineMiRNA::WebFunctions->resultstruct2pseudoXML(\%myResults);

    print <<"HTML";
$HTML_header    <body onload="main();">
        <div class="titreDiv"> Identification of miRNA/miRNA hairpins results:</div>
$HTML_additional
        <div id="table" ></div>
        <div id="singleCell"> </div>
$HTML_results
    <a href="./resultsAsCSV.pl?run_id=$id_job">Download as CSV</a>
    </body>
</html>
HTML

}else{
    print <<"HTML";
$HTML_header
    <body>
        <p>No results available for the given job identifier $id_job: $valid </p>
    </body>
</html>
HTML
}
