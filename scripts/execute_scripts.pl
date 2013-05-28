#!/usr/bin/perl -w

use File::Path 'rmtree';

$dirScript = "/var/www/arn/scripts/"; # chemin script
$dirImages = "/var/www/arn/images/"; # chemin images
$dirData = "/var/www/arn/data/"; # chemin séquence
$dirBlast = "/var/www/arn/scripts/ncbi-blast-2.2.26+/bin/"; # chemin Blast
$dirRandfold = "/var/www/arn/scripts/randfold-2.0/"; # chemin Randfold

my ($check,$mfei,$randfold, $SC,$align,$dirJob,$plant) = @ARGV;

if ($check eq "checked")
{
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
{ 	
	
	open tempFile,  '>'.$dirJob.'tempFile.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
	system('chmod 777 '.$dirJob.'tempFile.txt');
	print tempFile $name.@tab{$name};
	$name = substr $name, 1,-1 ;
	system('mkdir '.$dirJob.$name);
		
	system('chmod 777 '.$dirJob.$name.'/fichierRNAfold.txt');
	system($dirScript.'ViennaRNA-1.8.5/Progs/RNAfold -noPS < '.$dirJob.'/tempFile.txt > '.$dirJob.$name.'/fichierRNAfold.txt');
	system('chmod 777 '.$dirJob.$name.'/fichierRNAfold.txt');
	####conversion en format CT 
	system($dirScript.'ViennaRNA-1.8.5/Utils/b2ct < '.$dirJob.$name.'/fichierRNAfold.txt > '.$dirJob.$name.'/fichierOutB2ct.ct');
	system('chmod 777 '.$dirJob.$name.'/fichierOutB2ct.ct');
	system('perl '.$dirScript.'Tigeboucle3bc.pl  '.$dirJob.$name.'/fichierOutB2ct.ct '.$dirJob.$name);
	print tempFile $dirScript.'Tigeboucle3bc.pl  '.$dirJob.$name.'/fichierOutB2ct.ct '.$dirJob.$name;
	unlink $dirJob."tempFile.txt";
	
}

##Traitement fichier de sortie outStemloop
opendir DIR, $dirJob; #ouverture répertoire job
my @dirs;
@dirs=readdir DIR;
closedir DIR;
foreach $dir(@dirs) # parcours du contenu
{
	if ($dir ne "." && $dir ne ".." && -d $dirJob.$dir) #si fichier est un répertoire
	{
		
		opendir DIR, $dirJob.$dir; # ouverture du sous répertoire 
		my @files;
		@files=readdir DIR;
		closedir DIR;
		foreach $file(@files)
		{
			
				if ($file ne "." && $file ne ".." && -d $dirJob.$dir."/".$file) # si le fichier est de type repertoire
				{
					####Traitement fichier de sortie outStemloop
					chmod 0777, $dirJob.$dir."/".$file;
					#system('perl /var/www/arn/scripts/image.pl '.$dirScript.' '.$dirJob.$dir."/".$file.'/');
					system($dirScript.'ViennaRNA-1.8.5/Progs/RNAfold -noPS < '.$dirJob.$dir."/".$file.'/seq.txt > '.$dirJob.$dir."/".$file.'/outRNAFold.txt');
					system('chmod 777 '.$dirJob.$dir."/".$file.'/outRNAFold.txt');
					#system ('/usr/bin/java -cp '.$dirScript.'VARNAv3-8.jar fr.orsay.lri.varna.applications.VARNAcmd -i '.$dirJob.$dir.'/'.$file.'/outB2ct.ct -o '.$dirJob.$dir.'/'.$file.'/image.png'); 
					
					
					#open coucou, '>>'.$dirJob.$dir.'/'.$file.'/coucou.txt';
					#print coucou 'java -cp '.$dirScript.'VARNAv3-8.jar fr.orsay.lri.varna.applications.VARNAcmd -i '.$dirJob.$dir."/".$file.'/outB2ct.ct -o '.$dirJob.$dir."/".$file.'/image.png';
					#print coucou '\nperl /var/www/arn/scripts/image.pl '.$dirScript.' '.$dirJob.$dir."/".$file.'/';
					####conversion en format CT 
					system($dirScript.'ViennaRNA-1.8.5/Utils/b2ct < '.$dirJob.$dir."/".$file.'/outRNAFold.txt > '.$dirJob.$dir."/".$file.'/outB2ct.ct');
					system('chmod 777 '.$dirJob.$dir."/".$file.'/outB2ct.ct');
					
					#system('./test.sh '.$dirScript.' '.$dirJob.$dir.'/'.$file.'/ '.' > trash 2> errors');
					system ('/usr/bin/java -cp '.$dirScript.'VARNAv3-8.jar fr.orsay.lri.varna.applications.VARNAcmd -i '.$dirJob.$dir.'/'.$file.'/outB2ct.ct -o '.$dirJob.$dir.'/'.$file.'/image.png > trash');
					## traitement du fichier OutVienna pour la récupération des données(Format Vienna, séquence ADN) 
					$dirSeq = $1;
					open (TRAITED, '>'.$dirJob.$dir."/".$file.'/outViennaTraited.txt') || die "$!";
					open (INPUT, $dirJob.$dir."/".$file.'/outRNAFold.txt') || die "$!";
					while (my $line = <INPUT>) 
					{
						if (($line=~/^>(.*)/))  {$nameSeq = $1;}#nom sequence
						elsif (($line=~/^[a-zA-Z]/))    {$dna= substr $line,0,-1;}#récupération de la sequence adn
						elsif (($line=~/(.*) /))        
						{ 
							$Vienna= $1;print TRAITED $nameSeq."\t".$dna."\t".$Vienna."\n";#récupération du format Vienna
						}
					 }
					close INPUT;
					close TRAITED;
					system('chmod 777 '.$dirJob.$dir."/".$file.'/outViennaTraited.txt');
					####calcul MFEI (appel script energie.pl)
					if ($mfei eq "mfeiChecked") 
					{
						system('perl '.$dirScript.'energie.pl '.$dirJob.$dir."/".$file.'/outB2ct.ct '.$dirJob.$dir.' '.$file);
					}     
					####calcul p-value randfold  
					if ($randfold eq "randfoldChecked") 
					{
						system($dirRandfold.'randfold -d '.$dirJob.$dir."/".$file.'/seq.txt 7 > '.$dirJob.$dir."/".$file.'/pvalue.txt');
						system('chmod 777 '.$dirJob.$dir."/".$file.'/pvalue.txt');
					}  
					####calcul self-contain   
					if ($SC eq "SCChecked") 
					{
						system('python '.$dirScript.'selfcontain_unix/selfcontain.py -i '.$dirJob.$dir."/".$file.'/seq.txt -n 100  > '.$dirJob.$dir."/".$file.'/selfContain.txt');
						system('chmod 777 '.$dirJob.$dir."/".$file.'/selfContain.txt');
					}
					####### creation sequence boucle terminale masquee avec des N pour chaque sequence (repertoire ) et resultat alignement mirBASE
					if ($align eq "alignChecked")
					{
						open seqN, '>>'.$dirJob.$dir.'/'.$file.'/seqWithN.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
						open seq, $dirJob.$dir.'/'.$file.'/seq.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
						while  ($line=<seq>)
						{
							if ($line=~/^>/)
							{
								print seqN $line;
							}
							else 
							{
								$line =~ s/U/T/g;
								print seqN $line;
							}
						}
						close seqN;
			
						system($dirScript.'exonerate-2.2.0-i386/bin/exonerate -E --model affine:bestfit '.$dirData.'MirbaseFile.txt  '.$dirJob.$dir.'/'.$file.'/seqWithN.txt -d '.$dirData.'matrix --bestn 1 --score -3 -e -1 -o -1 > '.$dirJob.$dir.'/'.$file.'/alignement.txt');
					}
				}
		}
	}
}



