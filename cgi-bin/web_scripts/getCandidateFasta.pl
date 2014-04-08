#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::Candidate;
use PipelineMiRNA::Results;
use PipelineMiRNA::WebTemplate;

my $cgi = CGI->new();

my $jobId    = $cgi->param('jobId');
my $candidate_id   = $cgi->param('id');

my $job            = PipelineMiRNA::Results->jobId_to_jobPath($jobId);

my %candidate;

if (
    !eval {
        %candidate =
          PipelineMiRNA::Candidate->retrieve_candidate_information( $job,
            $candidate_id );
    }
  )
{

    # Catching exception
    print PipelineMiRNA::WebTemplate::get_error_page(
        'No results for the given identifiers');
}
else {
    my $candidate_name = PipelineMiRNA::Candidate->get_shortened_name(\%candidate);
    my $fasta = PipelineMiRNA::Candidate->candidateAsFasta(\%candidate);
    print <<"DATA" or die "Error when printing content: $!";
Content-type: text/txt
Content-disposition: attachment;filename=$candidate_name.fa

$fasta
DATA
}
