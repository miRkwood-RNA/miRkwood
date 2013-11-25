#!/usr/bin/perl -w

my ($job, $title) = @ARGV;
if ($title){
    if ($title eq "noTitle")
    {
        $title= "";
    }
}else{
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
print MAIL <<"DATA";
Dear MicroARN program user,

Your MicroARN web job is completed.
Results are available at http://monprojet.com/cgi-bin/resultsWithID.pl?run_id=".$job."
Thank you for using program.\n";
DATA

close(MAIL);





