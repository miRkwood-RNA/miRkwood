#!/usr/bin/perl -w
use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use feature 'switch';

BEGIN { require File::Spec->catfile( $FindBin::Bin, 'requireLibrary.pl' ); }
use miRkwood::CandidateHandler;
use miRkwood::Results;
use miRkwood::WebTemplate;

my $cgi = CGI->new();

my $jobId        = $cgi->param('jobId');
my $candidate_id = $cgi->param('id');
my $export_type  = $cgi->param('type');
my $optimal      = $cgi->param('optimal');

my $job = miRkwood::Results->jobId_to_jobPath($jobId);

my $candidate;

if (
    !eval {
        $candidate =
          miRkwood::CandidateHandler->retrieve_candidate_information( $job,
            $candidate_id );
    }
  )
{

    # Catching exception
    print miRkwood::WebTemplate::get_error_page(
        'No results for the given identifiers');
}
else {
    my ( $filename, $contents, $disposition );
    given ($export_type) {
        when (/fas/) {
            ( $filename, $contents, $disposition ) =
              exportAsFasta( $candidate )
        }
        when (/dot/) {
            ( $filename, $contents, $disposition ) =
              exportAsDotBracket( $candidate, $optimal )
        }
        when (/alt/) {
            ( $filename, $contents, $disposition ) =
              exportAlternatives( $candidate )
        }
        when (/yml/) {
            ( $filename, $contents, $disposition ) =
              exportAsYAML( $candidate )
        }
        when (/reads/) {
            ( $filename, $contents, $disposition ) =
              exportReadsClouds( $candidate_id, $jobId )
        }
        default {
            miRkwood::WebTemplate::web_die("Error: the export type '$export_type' is not supported");
        }
    }
    print <<"DATA" or miRkwood::WebTemplate::web_die("Error when printing content: $!");
Content-type: text/plain
Content-disposition: $disposition;filename=$filename

$contents
DATA
}

sub exportAsFasta {
    my @args      = @_;
    my $candidate = shift @args;
    my $candidate_name =
      $candidate->get_shortened_name();
    my $fasta = $candidate->candidateAsFasta();
    return ( "$candidate_name.fa", $fasta, 'inline' );
}

sub exportAsDotBracket {
    my @args      = @_;
    my $candidate = shift @args;
    my $optimal   = shift @args;
    my $candidate_name =
      $candidate->get_shortened_name();
    my $filename = $candidate_name;
    my $vienna =
      $candidate->candidateAsVienna( $optimal );
    if ($optimal) {
        $filename .= "_optimal";
    }
    return ( "$filename.txt", $vienna, 'attachment' );
}

sub exportAlternatives {
    my @args      = @_;
    my $candidate = shift @args;
    my $filename  = $candidate->get_shortened_name() . '_alternatives';
    my $alternatives = $candidate->alternativeCandidatesAsVienna();
    return ( "$filename.txt", $alternatives, 'attachment' );
}

sub exportAsYAML {
    my @args      = @_;
    my $candidate = shift @args;
    my $filename  = $candidate->get_shortened_name();
    use YAML::XS;
    my $yaml = YAML::XS::Dump($candidate);
    return ( "$filename.yml", $yaml, 'attachment' );
}

sub exportReadsClouds {
    my @args           = @_;
    my $candidate_id   = shift @args;
    my $job_id         = shift @args;
    my $linkReadsCloud = miRkwood::CandidateHandler::get_candidate_reads_cloud_file( $job_id, $candidate_id );
    my $contents = '';
    open (my $IN, '<', $linkReadsCloud) or warn "ERROR when opening $linkReadsCloud : $!";
    while ( <$IN> ){
        $contents .= $_;
    }
    close $IN;
    return ( "$candidate->{'identifier'}.txt", $contents, 'attachment' );
}
