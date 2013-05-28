#!/usr/bin/perl -w

@tab=();

push(@tab,74);
push(@tab,8);
push(@tab,77);
push(@tab,1);

print <<DATA;
Content-type: application/xhtml+xml

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

	<head>
		
	</head>
	<body>
		<resutls id="toto">
DATA

for ($i=0;$i<scalar (@tab); $i++)
{
	print "<Sequence id='".$i."' adn='".$tab[$ii]."' crit1='9.8'  crit2='1.0'></Sequence>
	";
}
print <<DATA;		
	</resutls>
	</body>

</html>

DATA
###End###
