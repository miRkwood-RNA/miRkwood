#!/usr/bin/perl -w

my ($title, $job) = @ARGV;
if ($title eq "noTitle")
{
	$title= "";	
}
$to='benmounah.mohcen@gmail.com';
$from= 'miarn@lifl.fr';
$subject='MiARN web job';

open(MAIL, "|/usr/sbin/sendmail -t"); 

## Mail Header
print MAIL "To: $to\n"; 
print MAIL "From: $from\n";
print MAIL "Subject: $subject $title\n\n";
## Mail Body
print MAIL "Dear MicroARN program user,\n\nYour MicroARN web job is completed.\nResults are available at http://monprojet.com/cgi-bin/resultsWithID.pl?run_id=".$job." \n\nThank you for using program.\n";

close(MAIL);





