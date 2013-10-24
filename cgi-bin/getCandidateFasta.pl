#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;                       # locate this script
use lib "$FindBin::Bin/../lib";    # use the parent directory
use PipelineMiRNA::Candidate;
use PipelineMiRNA::Results;
use PipelineMiRNA::WebTemplate;

my $cgi = CGI->new();

my $jobId    = $cgi->param('jobId');
my $name     = $cgi->param('name');
my $position = $cgi->param('position');

my $candidate_name = $name . '__' . $position;
my $job            = PipelineMiRNA::Results->jobId_to_jobPath($jobId);

my %candidate;

if (
    !eval {
        %candidate =
          PipelineMiRNA::Candidate->retrieve_candidate_information( $job,
            $name, $candidate_name );
    }
  )
{

    # Catching exception
    print PipelineMiRNA::WebTemplate::get_error_page(
        'No results for the given identifiers');
}
else {
    my $fasta = PipelineMiRNA::Candidate->candidateAsFasta(\%candidate);
    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/txt
Content-disposition: attachment;filename=$candidate_name.fa

$fasta
DATA
}
