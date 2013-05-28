#!/usr/bin/perl -w
print "Content-type: text/html\n\n"; 
use CGI; 
my $cgi = new CGI; 

$dirData = "/var/www/data/"; # chemin séquence
$dirBlast = "/var/www/ncbi-blast-2.2.26+/bin/"; # chemin Blast
my $upload = $cgi->upload('seq') || die "$!";
open (INPUT, '>> '.$dirData.'sequence.fas') || die "$!";
while (my $ligne = <$upload>) 
{
	print INPUT $ligne;
}
close INPUT;
system('chmod 777 '.$dirData.'sequence.fas');
$plante =  $cgi->param('db'); # récupération paramètre menu déroulant
print $plante;
print system($dirBlast.'blastx -query '.$dirData.'sequence.fas -db '.$dirData.$plante.'.fas -outfmt 6 -max_target_seqs 1 -evalue 1E-5 -out '.$dirData.'outBlast.txt');
system('chmod 777 '.$dirData.'outBlast.txt');
##Filtrage du ficher : Elimination des régions codantes###
@SeqNamesOUT=(); # tableau contenant la liste des séquences issues du Blast
@SeqNames=(); # tableau contenant la liste des séquences à traiter
@SeqDiff=(); # tableau contenant la liste des séquences non codantes
$i=0;

# Ouverture du fichier outBlast.TXT et construction du tableau contenant la liste des séquences issues du blast
open(FOut, $dirData.'outBlast.txt') || die "Problème à l\'ouverture : $!";
while(my $line = <FOut>)
{
	my @name = split('\t', $line);
	push(@SeqNamesOUT,@name[0]);
}
close FOut || die "Problème à la fermeture : $!";

# Ouverture du fichier sequence.fas et construction du tableau contenant la liste des séquences issues du blast
open(FSeq, $dirData.'sequence.fas') || die "Problème à l\'ouverture : $!";
while(my $line = <FSeq>)
{
	if ( grep( /^>/,$line) )
	{
		$lineSeq = substr $line, 1,-1;	
		push(@SeqNames,$lineSeq);
	}
	
}
close FSeq || die "Problème à la fermeture : $!";

#Remplissage du tableau @SeqDiff avec les séquences non codantes
foreach my $SeqName (@SeqNames) {
	$exists = true;
	foreach my $SeqNameOUT (@SeqNamesOUT) {
		if ($SeqName eq $SeqNameOUT)
		{
			$exists = false;
		}
	}
	if ($exists eq true)
	{
		push(@SeqDiff,$SeqName);
		#print "\nUne seq ".$SeqName."\n";
	}	
	$i = $i++;	
}

# Création du fichier FASTA filtrée des régions codantes
$record = false;
$i = 0;
open (RES, '>>'.$dirData.'SequencesFiltreesCDS.txt');
open(FSeq, $dirData.'sequence.fas') || die "Problème à l\'ouverture : $!";
while(my $line = <FSeq>)
{
	if ( grep( /^>/,$line) )
	{
		$lineSeq = substr $line, 1,-1;
		if (@SeqDiff[$i] eq $lineSeq)
		{
			print @SeqDiff[$i];
			$record = true;		
			$i++;
		}
		else
		{
			$record = false;
		}		
	}
		if ($record eq true)
		{			
			printf RES $line;
		}
}
close FSeq || die "Problème à la fermeture : $!";
close (RES);
system('chmod 777 '.$dirData.'SequencesFiltreesCDS.txt');

print <<DATA

<html>
	<head>
		<LINK rel="stylesheet" type="text/css" href="../css/script.css">
		<script src="../js/miARN.js" type="text/javascript" LANGUAGE="JavaScript"></script>
		<title>Identification des microARN</title>
	</head>
	<body>
		<h2> Identification des microARN </h2>
		<h2>  </h2>
		
	</body>
</html>

DATA
