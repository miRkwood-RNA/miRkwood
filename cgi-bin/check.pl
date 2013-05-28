#!/usr/bin/perl -w
print "Content-type: text/html\n\n"; 
use CGI; 
my $cgi = new CGI; 
my $seq = $cgi->param('seqArea');
	
if ($seq !~ /^( *>.+[\r\n]+([-\. atcgu0-9]+[\r\n]+)+){2,}$/)
{
	print "error";
}
print <<DATA
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html>
	<head>
		
		<title>Identification des microARN</title>
	</head>
	<body>
		<h2> Identification des microARN </h2>
		<h2>  </h2>
		
	</body>
</html>

DATA
