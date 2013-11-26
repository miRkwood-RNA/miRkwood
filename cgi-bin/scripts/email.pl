#!/usr/bin/perl -w

my ($results_baseurl, $jobId, $mail, $title) = @ARGV;

my $res_arguments = '?run_id=' . $jobId;
my $results_url   = $results_baseurl . $res_arguments;


$from= 'miarn@lifl.fr';
$subject='MiARN web job';

open(MAIL, "|/usr/sbin/sendmail -t"); 

## Mail Header
print MAIL "To: $mail\n"; 
print MAIL "From: $from\n";
print MAIL "Subject: $subject $title\n\n";
## Mail Body
print MAIL <<"DATA";
Dear MicroARN program user,

Your MicroARN web job is completed.
Results are available at $results_url
Thank you for using program.\n";
DATA

close(MAIL);




