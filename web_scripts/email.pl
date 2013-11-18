#!/usr/bin/perl
print "Content-type: text/html\n\n";

$title='Perl Mail demo';
$to='benmounah.mohcen@gmail.com';
$from= 'benmounah.mohcen@gmail.com';
$subject='YOUR SUBJECT';

open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: $subject\n\n";
## Mail Body
print MAIL "This is a test message from Cyberciti.biz! You can write your
mail body text here\n";

close(MAIL);

print "<html><head><title>$title</title>
</head>\n<body>\n\n";

## HTML content let use know we sent an email
print "<h1>$title</h1>\n";
print "<p>A message has been sent from $from to $to";
print "\n\n</body></html>";
