#!/usr/bin/perl -w



$dirScript = "/home/mohcen/arn/scripts/"; # chemin script
$dirImages = "/home/mohcen/arn/images/"; # chemin images
$dirData = "/home/mohcen/arn/data/"; # chemin séquence
$dirBlast = "/home/mohcen/arn/scripts/ncbi-blast-2.2.26+/bin/"; # chemin Blast
$dirRandfold = "/home/mohcen/arn/scripts/randfold-2.0/"; # chemin Randfold
#appel script upload.pl qui charge le fichier dans le serveur
#system('perl '.$dirScript.'upload.pl '.$dirData.' '.$upload); 
my ($check,$mfei,$randfold,$dirJob,$plant) = @ARGV;




if ($check eq "checked")
{
	 # récupération paramètre menu déroulant
	#appel script filterCDS.pl qui permet de filter les CDS
	system('perl '.$dirScript.'filterCDS.pl '.$dirBlast.' '.$dirData.' '.$dirJob.' '.$plant);
}
else
{
	system('mv '.$dirJob.'sequenceUpload.fas '.$dirJob.'Sequences.fas');
}

##Passage du multifasta -> fasta et appel au script Stemloop
open ENTREE, $dirJob.'Sequences.fas' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
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
{ 	open tempFile,  '>'.$dirJob.'tempFile.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
open test,  '>>'.$dirJob.'test.txt';
	system('chmod 777 '.$dirJob.'tempFile.txt');
	print tempFile $name.@tab{$name};
	$name = substr $name, 1,-1 ;
	#mkdir ($dirJob.$name);
	system('mkdir '.$dirJob.$name);

	#open outStemTraited,  '>'.$dirJob.$name.'/outStemTraited.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
	print test '-p75 -o '.$dirJob.$name."/out.txt"	;
	system ($dirScript.'Stemloop/stemloop -f '.$dirJob.'tempFile.txt -p75 -o '.$dirJob.$name."/outStemLoop.txt");
	unlink $dirJob."tempFile.txt";
	
	system('chmod 777 '.$dirJob.$name.'/outStemLoop.txt');
}

##Traitement fichier de sortie outStemloop

opendir DIR, $dirJob;
my @dirs;
@dirs=readdir DIR;
closedir DIR;
foreach $dir(@dirs)
{
	if ($dir ne "." && $dir ne ".." && -d $dirJob.$dir) 
	{
		print $dir."\n";
		opendir DIR, $dirJob.$dir;
		my @files;
		@files=readdir DIR;
		closedir DIR;
		foreach $file(@files)
		{
			
			if ($dirJob.$dir."/".$file ne "." && $dirJob.$dir."/".$file ne ".." && -f $dirJob.$dir."/".$file) 
			{
				print $file."\n";
				open outStem, $dirJob.$dir.'/outStemLoop.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
				open outStemTraited, '>'.$dirJob.$dir.'/outStemTraited.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
				while  ($line=<outStem>)
				{
					if (($line=~/^>/)) 
					{
						$line=~s/\s+\[/\_/;
						$line=~s/\]//;
						print outStemTraited $line;
					} 
					else
					{
						
						print outStemTraited $line;
					}	
				}
				system('chmod 777 '.$dirJob.$dir.'outStemTraited.txt');



			}


		}



	}


}

=header
####Traitement fichier de sortie outStemloop
system($dirScript.'ViennaRNA-1.8.5/Progs/RNAfold -noPS < '.$dirData.'outStemTraited.txt > '.$dirData.'outRNAFold.txt');
system('chmod 777 '.$dirData.'outRNAFold.txt');
	
####conversion en format CT 
system($dirScript.'ViennaRNA-1.8.5/Utils/b2ct < '.$dirData.'outRNAFold.txt > '.$dirData.'outB2ct.ct');
system('chmod 777 '.$dirData.'outB2ct.ct');

####détection des séquences représentants une structure tige boucle (appel script tigeBoucle.pl) et génération d'un nouveau fichier ct filtrés
system('perl '.$dirScript.'tigeBoucle.pl '.$dirData.'outB2ct.ct '.$dirData);
	
####calcul MFEI (appel script energie.pl) 

if ($mfei eq "mfeiChecked") 
{
	system('perl '.$dirScript.'energie.pl '.$dirData.'outTigeBoucle.ct '.$dirData);
}
	
####Filtrage fichier MULTIFASTA des séquences non tiges boucles  (appel script filterFasta.pl) 
system('perl '.$dirScript.'filterFasta.pl '.$dirData.'outTigeBoucle.ct '.$dirData.'outStemTraited.txt '.$dirData);

#### Génération structure sequence en post scritp et le fichier contenant format Vienna
system($dirScript.'ViennaRNA-1.8.5/Progs/RNAfold < '.$dirData.'tigeBouclesSequence.txt > '.$dirData.'outVienna.txt');
system('chmod 777 '.$dirData.'outVienna.txt');
## traitement du fichier OutVienna pour la récupération des données(Format Vienna, sequence ADN) 
open (TRAITED, '>>'.$dirData.'outViennaTraited.txt') || die "$!";

open (INPUT, $dirData.'outVienna.txt') || die "$!";
	while (my $line = <INPUT>) 
	{
		if (($line=~/^> (.*) \(/))	{$nameSeq = $1;}#nom sequence
		elsif (($line=~/^[a-zA-Z]/))	{$dna= substr $line,0,-1;}#récupération de la sequence adn
		elsif (($line=~/(.*) /))	{ 
			$Vienna= $1;print TRAITED $nameSeq."\t".$dna."\t".$Vienna."\n";#récupération du format Vienna
		}
	}
	close INPUT;
	close TRAITED;
system('mv ./*ss.ps '.$dirImages);
system('chmod 777 '.$dirData.'outViennaTraited.txt');

####calcul p-value randfold  

if ($randfold eq "randfoldChecked") 
{
	system($dirRandfold.'randfold -d '.$dirData.'tigeBouclesSequence.txt 7 > '.$dirData.'pvalue.txt');
	system('chmod 777 '.$dirData.'pvalue.txt');
 }
=cut
