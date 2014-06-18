#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use feature 'switch';

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use PipelineMiRNA::Candidate;
use PipelineMiRNA::Results;
use PipelineMiRNA::WebTemplate;

my $cgi = CGI->new();

my $jobId        = $cgi->param('jobId');
my $candidate_id = $cgi->param('id');
my $export_type  = $cgi->param('type');
my $optimal      = $cgi->param('optimal');

my $job = PipelineMiRNA::Results->jobId_to_jobPath($jobId);

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
    my ( $filename, $contents, $disposition );
    given ($export_type) {
        when (/fas/) {
            ( $filename, $contents, $disposition ) =
              exportAsFasta( \%candidate )
        }
        when (/dot/) {
            ( $filename, $contents, $disposition ) =
              exportAsDotBracket( \%candidate, $optimal )
        }
        when (/alt/) {
            ( $filename, $contents, $disposition ) =
              exportAlternatives( \%candidate )
        }
        when (/yml/) {
            ( $filename, $contents, $disposition ) =
              exportAsYAML( \%candidate )
        }
        default {
            PipelineMiRNA::WebTemplate::web_die("Error: the export type '$export_type' is not supported");
        }
    }
    print <<"DATA" or PipelineMiRNA::WebTemplate::web_die("Error when printing content: $!");
Content-type: text/plain
Content-disposition: $disposition;filename=$filename

$contents
DATA
}

sub exportAsFasta {
    my @args      = @_;
    my $candidate = shift @args;
    my $candidate_name =
      PipelineMiRNA::Candidate->get_shortened_name($candidate);
    my $fasta = PipelineMiRNA::Candidate->candidateAsFasta($candidate);
    return ( "$candidate_name.fa", $fasta, 'inline' );
}

sub exportAsDotBracket {
    my @args      = @_;
    my $candidate = shift @args;
    my $optimal   = shift @args;
    my $candidate_name =
      PipelineMiRNA::Candidate->get_shortened_name($candidate);
    my $filename = $candidate_name;
    my $header   = ">$candidate_name";
    my $vienna =
      PipelineMiRNA::Candidate->candidateAsVienna( $candidate, $optimal );
    if ($optimal) {
        $filename .= "_optimal";
    }
    return ( "$filename.txt", $vienna, 'attachment' );
}

sub exportAlternatives {
    my @args      = @_;
    my $candidate = shift @args;
    my $filename  = PipelineMiRNA::Candidate->get_shortened_name($candidate)
      . '_alternatives';
    my $alternatives =
      PipelineMiRNA::Candidate->alternativeCandidatesAsVienna($candidate);
    return ( "$filename.txt", $alternatives, 'attachment' );
}

sub exportAsYAML {
    my @args      = @_;
    my $candidate = shift @args;
    my $filename  = PipelineMiRNA::Candidate->get_shortened_name($candidate);
    use YAML::XS;
    my $yaml = YAML::XS::Dump($candidate);
    return ( "$filename.yml", $yaml, 'attachment' );
}
