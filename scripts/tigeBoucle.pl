#!/usr/local/bin/perl
######################################################################
# Script permettant de tester si éventuellement la séquence peut     #
#  représenter une structure tige boucle ou pas. Permet également de #
#  générer un nouveau fichier ct contenant que les séquences         #
#   représentant une structure tige boucle                           #
#                                                                    #
#	 						             							 #
#  @author : Mohcen BENMOUNAH                                        #
#                                                                    #
#                                                                    #
######################################################################
@tab=();
my ($CT,$dirData) = @ARGV;
@tigeBoucle = ();
open(FOut, $CT) || die "Problème à l\'ouverture : $!"; #ouverture du fichier de sortie de RNAFOLD
$ordre = false;
while(my $ligne = <FOut>)
{
	chomp $ligne;
	$ligne=~s/\s+//; #second type de fichier ct, permet d'éviter le décalage des colonnes par rapport aux sorties de RNAstructure.
	@tbl1=split /\s+/, $ligne;
	#Recherche des identifiants
	if ($ligne=~/ENERGY/)
	{#####Seuils à appliquer
		if ($lgtige>22 and $lgtige<100 and $perc>70 and $boucleterm<150 and $rep==0) 
		{
			push(@tigeBoucle, $name);#Stocke les identifiants des séquences qui passent le test
		}
		$name=$tbl1[4];
		$val=0; #numero base de 1 à longueur du premiR
		$app=0; #numero base appariée en face de $val
		$count0=0;#nombre de 0 du début jusque la boucle terminale comprise.
		$count=0;
		$test=1;
		$est=0;
		$count2=0; #nb bases non appariées en 5'
		$test3=1;
		$lgtige=0;
		$perc=0;
		$rep=0;
	}
	else 
	{#Si pas une ligne d'identifiant
		if($test2==1)
		{ #Considérer la tige qu'à partir du premier appariement
			if($test==1)
			{ #Ne passer ici qu'une fois dès la première base de l'autre coté de la tige
				if ((($tbl1[4]<$tbl1[0]) and ($tbl1[4]!=0)) or ($tbl1[4]<($app-7) and ($tbl1[4]!=0)))
				{#cette ligne conditionne la fin de la tige 1) en est passé de l'autre coté de la tige 2) il y a une 					boucle de l'autre coté
					$test=0;#remet la valeur test à 0 pour ne plus y passer
					$lgtige=$val-$count2;#longueur de la tige = position dernière base appariée - nb de 0 avant la tige
					$perc=($lgtige-($count0-$count2+$val-$tbl1[0]+1))/$lgtige*100; #nombre de 0 comptés moins les 0 de la première boucle boucle terminale
					$boucleterm=$app-$val-1; #longueur de la boucle à l'extremité de la tige
				}

			}
		}
		if ($tbl1[4]!=0)
		{ #Si une base est appariée
			if ($test3==1)
			{#Passe ici une seule fois = première base appariée
				$count2=$count0;#compte le nombre de bases non appariées en 5' : $count2;
				$test3=0;
				$prembaseapp=$tbl1[4];#Stocke le numéro de la première base appariée
			}
			$val=$tbl1[0]; #numéro Base
			$app=$tbl1[4]; #numéro Base appariée
			$test2=1;#Active le début de la tige : A partir de là on est dans la tige
		}
		else
		{
			$count0++;#compte les 0 : également ceux de la boucle terminale et ceux non appariés au début
		}
		#Recherche d'une boucle en 3' Observée car base appariée > première base appariée	
		if ($tbl1[4]>$prembaseapp)
		{
			$rep=1;
		}
	}
}	
close FOut || die "Problème à la fermeture : $!";
# Création du fichier CT en gardant que les séquences à structure tige boucle 
$record = false;
$i = 0;
open(FOut, $CT) || die "Problème à l\'ouverture : $!"; #ouverture du fichier de sortie de Unafold
open (RES, '>> '.$dirData.'/outTigeBoucle.ct');
while(my $line = <FOut>)
{
	my @variable = split('\s+', $line);
	if (@variable[2]=~ /^ENERGY/)
	{	
		$nomSeq = substr @variable[5], 0; # Récupération du nom de la séquence
		if (@tigeBoucle[$i] eq $nomSeq)
		{	
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
			print RES $line;
		}
}
system('chmod 777 '.$dirData.'/outTigeBoucle.ct');
close FOut || die "Problème à la fermeture : $!";
close (RES)|| die "Problème à la fermeture : $!";
