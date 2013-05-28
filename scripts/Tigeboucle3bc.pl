#!/usr/local/bin/perl
use File::Path 'rmtree';

$eraseDir = True;

my ($CT,$dirData) = @ARGV;
$seq = "";

open(FOut, $CT) || die "Problème à l\'ouverture : $!"; #ouverture du fichier de sortie de RNAFOLD
while($ligne = <FOut>)
{	
	@tbl1=split /\s+/, $ligne;
	if (! ($tbl1[2]=~/ENERGY/))
	{
		$seq.=$tbl1[2];
	}
}
close FOut;

open(FOut, $CT) || die "Problème à l\'ouverture : $!"; #ouverture du fichier de sortie de RNAFOLD
while($ligne = <FOut>)
{	
	@tbl1=split /\s+/, $ligne;
    #recherche des identifiants/début d'une nouvelle séquence à traiter
    if ($ligne=~/ENERGY/)
    {
		$ID=$tbl1[5];#Mise en mémoire de l'identifiant de la séquence
		$newtige=1; #on va démarré une nouvelle tige forcemment
		$test1=0;#On ne vient pas d'un cas 3 (tige dans une tige qui conditionne de mémoriser la base précédente et d'enlever un appariement)
		$test2=1;#pour calcul longueur boucle cas2
		$nbapp=0;#mise à zero du nombre d'appariements
	}#if
	#si on n'est pas dans une ligne d'identifiant, nous sommes donc dans le tableau des appariements
	else 
	{
		#Démarrage nouvelle tige activé,pas 0 et on est sur un coté 5' d'une tige
		if (($newtige==1) and ($tbl1[5]!=0) and ($tbl1[5]>$tbl1[1]))
		{
			#si on a une tige ds une tige c'est la position précédente que l'on mémorise puisque là on est déjà sur la seconde base de la tige suivante, donc on ne passe pas par la boucle = test<>0
			if($test1==0)
			{			
				$idep=$tbl1[1];
				$BPidep=$tbl1[5];
				$nbapp=0;#remise à zero des appariements
			}
			$newtige=0; #on est plus dans une nouvelle tige on vient de l'initier
			#mise en mémoire des bases i et BPi comme étant les dernières
			$ider=$tbl1[1];
			$BPider=$tbl1[5];		
		}
		#On ne démarre pas une nouvelle tige ; on continue dans la tige
		if ($newtige==0)
		{
			###on est sur le 2ème brin de la tige : cas2 :c'est une boucle terminale simple
			if (($tbl1[5]<$tbl1[1]) and ($tbl1[5]!=0))
			{ 
				#mise en mémoire de la première base après la boucle terminale : permet de calculer la longueur de la tige, la longueur de la boucle terminale et le nombre d'appariements
				if($test2==1)
				{#test2 pour ne mémoriser que la première base BPi<i
					$lgtige=$tbl1[5]-$idep+1;
					$lgboucle=$tbl1[1]-$tbl1[5]-1;
					$test2=0;
									
					#si on vient de cas 3 (tige dans une tige) il faut enlever un appariemment
					if ($test1==1)
					{
						$nbapp=$nbapp-1;
					}
					$long=$BPidep-$idep+1; #longueur de la tige
					$percapp=$nbapp/$lgtige*100; #pourcentage d'appariement
					#paramètres de sélection des tiges
					if (($long>80) and ($long<500) and ($lgboucle<150) and ($lgtige>30) and ($percapp>64))
					{
						system('mkdir '.$dirData."/".$ID."__$idep-$BPidep");
						open tempFile,  '>>'.$dirData."/".$ID."__$idep-$BPidep".'/seq.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
						$mySeq = substr ( $seq,$idep-1,$BPidep-$idep+1);
						print tempFile '>'.$ID."__$idep-$BPidep\n$mySeq\n";
						$eraseDir = False;
					}
					$newtige=1;#On va avoir une nouvelle tige
					$nbapp=0; #remise à zero du nombre d'appariements	
					$test1=0; #Pas passage par un cas3 (tige dans tige)
					$test2=1; #reactivation de la boucle de mise en mémoire de la première base après une boucle terminale simple	
				}
			}
			
			
			#on a une tige dans la tige cas 3 : conditions : (1)saut de BPi par rapport à la dernière base appariée en + ou en - (2) On a pas zero (3) on est coté 5' de la tige 
			if ((($tbl1[5]<$BPider-10) or ($tbl1[5]>$BPider+10)) and ($tbl1[5]!=0) and ($tbl1[5]>$tbl1[1]))
			{
				#si on vient de cas 3 on doit enlever un appariement
				if ($test1==1)
				{
					$nbapp=$nbapp-1;
				}
				$long=$BPidep-$idep+1; #longueur du premiR
				$lgboucle=$BPider-$ider-1; #longueur de la boucle terminale
				$lgtige=$ider-$idep+1; #longueur de la tige
				$percapp=$nbapp/$lgtige*100; #% d'appariement
				#paramètres de sélection des tiges
				if (($long>80) and ($long<500) and ($lgboucle<150) and ($lgtige>30) and ($percapp>64))
				{
					system('mkdir '.$dirData."/".$ID."__$idep-$BPidep");
					open tempFile,  '>>'.$dirData."/".$ID."__$idep-$BPidep".'/seq.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
					$mySeq = substr ( $seq,$idep-1,$BPidep-$idep+1);
					print tempFile '>'.$ID."__$idep-$BPidep\n$mySeq\n";
					$eraseDir = False;
				}
				$newtige=1;	#On va avoir une nouvelle tige
				$idep=$tbl1[1]; #mise en mémoire ici de idep (base de départ) on est déjà sur la première base de la tige suivante, contrairement au cas 2 où c'est la base suivante qui est la base de départ
				$BPidep=$tbl1[5]; #idem pour BPidep (base appariée de départ)
				$nbapp=1; #on vient déjà de voir le premier appariement de la tige suivante
				$test1=1;# Pour après un cas cas3 pour ne pas mémoriser la ligne suivante mais considérer cette position comme départ, et enlever un appariement
			}
		}
		if ($tbl1[5]!=0)
		{# pour chaque base différente de zero, on mémorise la dernière base lu (Position ider et appariement BPider)
			$BPider=$tbl1[5];
			$ider=$tbl1[1];
			$nbapp++;# on incrément le nombre d'appariements		
		}
	}
}	
