#!/usr/local/bin/perl
use File::Path 'rmtree';

$eraseDir = True;

my ($CT,$dirData) = @ARGV;

#fichier d'entrée contenant la sortie de RNAfold (.ct)

open(FOut, $CT) || die "Problème à l\'ouverture : $!"; #ouverture du fichier de sortie de RNAFOLD
while($ligne = <FOut>)
{	
	@tbl1=split /\s+/, $ligne;
    #recherche des identifiants/début d'une nouvelle séquence à traiter
    if ($ligne=~/ENERGY/)
    {
		#Mise en mémoire de l'identifiant de la séquence
		$ID=$tbl1[5];
		$newtige=1; #on démarre une nouvelle tige
		$test=1; #on 
		$nbapp=0;
		$seq="";
	}#if
	#si on n'est pas dans une ligne d'identifiant, nous sommes donc dans le tableau des appariements
	else 
	{
		$seq.=$tbl1[2];
		#Démarrage nouvelle tige et pas 0
		if (($newtige==1) and ($tbl1[5]!=0) and ($tbl1[5]>$tbl1[1]))
		{
			#si on a une tige ds une tige c'est la position précédente que l'on mémorise
			if($test!=0)
			{			
				$idep=$tbl1[1];
				$BPidep=$tbl1[5];
				$nbapp=0;
			}
			$newtige=0; #on est plus dans une nouvelle tige on vient de l'initier
			$ider=$tbl1[1];
			$BPider=$tbl1[5];
		}
		#Démarre pas une nouvelle tige		
		if ($newtige==0)
		{
			#on est sur le 2ème brin de la tige : cas2
			if (($tbl1[5]<$tbl1[1])  and ($tbl1[5]==$idep) and ($tbl1[5]!=0))
			{
				#si on vient de cas 3
				if ($test2==1)
				{
					$nbapp=$nbapp-1;
				}
				$nbapp=$nbapp+1;
				$nbapp=$nbapp/2;
				$long=$BPidep-$idep+1;
				if (($long>50) and ($long<500) and ($nbapp>25))
				{
					system('mkdir '.$dirData."/".$ID."_$idep-$BPidep");
					open tempFile,  '>>'.$dirData."/".$ID."_$idep-$BPidep".'/seq.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
					$mySeq = substr $seq,$idep,$BPidep-$idep;
					print tempFile '>'.$ID."_$idep-$BPidep\n$mySeq";
					$eraseDir = False;
				} 
				$newtige=1;
				$test=1;
				$nbapp=0;	
				$test2=0;
					
			}
			
			
			#on a une tige dans la tige cas 3
			if ((($tbl1[5]<$BPider-6) or ($tbl1[5]>$BPider+6)) and ($tbl1[5]!=0) and ($tbl1[5]>$tbl1[1]))
			{
				#si on vient de cas 3
				if ($test2==1)
				{
					$nbapp=$nbapp-1;
				}
				$long=$BPidep-$idep+1;
				if (($long>50) and ($long<500) and ($nbapp>25)){
				system('mkdir '.$dirData."/".$ID."_$idep-$BPidep");
					open tempFile,  '>>'.$dirData."/".$ID."_$idep-$BPidep".'/seq.txt' or die "Impossible d'ouvrir le fichier d'entree  : $!"; 
					$mySeq = substr $seq,$idep,$BPidep-$idep;
					print tempFile '>'.$ID."_$idep-$BPidep\n$mySeq";
					$eraseDir = False;
				}
				$newtige=1;  	
				$idep=$tbl1[1];
				$BPidep=$tbl1[5];
				$test=0;# pour ne pas mémoriser la ligne suivante mais considérer cette position comme départ
				$nbapp=1;
				$test2=1;# Pour après un cas cas3 
			}
		}
		if ($tbl1[5]!=0){
			$BPider=$tbl1[5];
			$nbapp++;
			#print "$nbapp\t $BPider\n";
		}
	}
}
	
if ($eraseDir eq True)
{
	rmtree($dirData);	# supprimer repertoire 
}		
		
		


