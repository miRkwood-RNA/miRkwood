#!/usr/bin/perl -w

use CGI; 
my $cgi = new CGI; 
$dirScript = "/var/www/scripts/"; # chemin script
$dirData = "/var/www/data/"; # chemin séquence
$dirBlast = "/var/www/scripts/ncbi-blast-2.2.26+/bin/"; # chemin Blast
$dirRandfold = "/var/www/scripts/randfold-2.0/"; # chemin Randfold
#appel script upload.pl qui charge le fichier dans le serveur
#system('perl '.$dirScript.'upload.pl '.$dirData.' '.$upload); 

my $upload = $cgi->upload('seq') || die "$!";
open (INPUT, '>> '.$dirData.'sequence.fas') || die "$!";
while (my $ligne = <$upload>) 
{
	print INPUT $ligne;
}
close INPUT;
system('chmod 777 '.$dirData.'sequence.fas');

$plant =  $cgi->param('db'); # récupération paramètre menu déroulant
#appel script filterCDS.pl qui permet de filter les CDS
system('perl '.$dirScript.'filterCDS.pl '.$dirBlast.' '.$dirData.' '.$plant);


##Passage du multifasta -> fasta et appel au script Stemloop
open ENTREE, $dirData.'noCdsSequence.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
%tab=();
while  ($line=<ENTREE>)
{
	if ( grep( /^>/,$line) )
	{
		$nameSeq = substr $line,0;	
	}
	else
	{
		$tab{$nameSeq}= $line;
	}
}
close ENTREE;	
foreach $name (keys %tab)
{ 	open tempFile,  '>'.$dirData.'tempFile.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
	system('chmod 777 '.$dirData.'tempFile.txt');
	print tempFile $name.@tab{$name};
	system ($dirScript.'Stemloop/stemloop -f '.$dirData.'tempFile.txt -p75  >>'.$dirData.'outStemLoop.txt');
	unlink $dirData."tempFile.txt";
}
system('chmod 777 '.$dirData.'outStemLoop.txt');

##Traitement fichier de sortie outStemloop
open outStem, $dirData.'outStemLoop.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
open outStemTraited, '>'.$dirData.'outStemTraited.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
while  ($line=<outStem>)
{
	if (($line=~/^[^M\n]/)) 
	{
		$line=~s/\s\[/\_/;
		$line=~s/\]//;
		print outStemTraited $line;
	} 	
}
system('chmod 777 '.$dirData.'outStemTraited.txt');

####Traitement fichier de sortie outStemloop
system($dirScript.'ViennaRNA-1.8.5/Progs/RNAfold -noPS < '.$dirData.'outStemTraited.txt > '.$dirData.'outRNAFold.txt');
system('chmod 777 '.$dirData.'outRNAFold.txt');
	
####conversion en format CT 
system($dirScript.'ViennaRNA-1.8.5/Utils/b2ct < '.$dirData.'outRNAFold.txt > '.$dirData.'outB2ct.ct');
system('chmod 777 '.$dirData.'outB2ct.ct');

####détection des séquences représentants une structure tige boucle (appel script tigeBoucle.pl) et génération d'un nouveau fichier ct filtrés
system('perl '.$dirScript.'tigeBoucle.pl '.$dirData.'outB2ct.ct '.$dirData);
	
####calcul MFEI (appel script energie.pl) 
system('perl '.$dirScript.'energie.pl '.$dirData.'outTigeBoucle.ct '.$dirData);

	
####Filtrage fichier MULTIFASTA des séquences non tiges boucles  (appel script filterFasta.pl) 
system('perl '.$dirScript.'filterFasta.pl '.$dirData.'outTigeBoucle.ct '.$dirData.'outStemTraited.txt '.$dirData);

####calcul p-value randfold  

system($dirRandfold.'randfold -d '.$dirData.'tigeBouclesSequence.txt 7 > '.$dirData.'pvalue.txt');
system('chmod 777 '.$dirData.'pvalue.txt');

