#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename qw(dirname);
use File::Spec;
use CGI::Carp qw(fatalsToBrowser);
use CGI;
my $cgi = CGI->new;
use Cwd qw( abs_path );
use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory
use PipelineMiRNA::WebFunctions;
use PipelineMiRNA::GFF;

my $chaine = q{};
my $id_job = $cgi->param('run_id');    # récupération id job

my $valid = PipelineMiRNA::WebFunctions->is_valid_jobID($id_job);

if ($valid) {
    my $gff = PipelineMiRNA::GFF->generate_GFF($id_job);

    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/gff
Content-disposition: attachment;filename=Results-$id_job.gff

$gff
DATA
}
else {
    die('Error: the job ID provided is not valid');
}
