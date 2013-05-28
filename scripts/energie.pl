#!/usr/local/bin/perl
######################################################################
#  Lancement du script : perl energie.pl FichierFasta				 #		                                         #
#  @author : Mohcen BENMOUNAH                                        #                                                        #
#                                                                	 #										         #
#                                                                    #                                                                                           #
######################################################################
my ($CT, $dirData, $seq ) = @ARGV;


open(FOut, $CT) || die "Problème à l\'ouverture : $!"; #ouverture du fichier de sortie de Unafold
$cg = 0;
while(my $line = <FOut>)
{
	my @variable = split('\s+', $line);
	{	
		if (@variable[2]=~ /^ENERGY/)
		{	
			$nameSeq = substr @variable[5], 0; # Récupération du nom de la séquence
			$mfe = substr @variable[4],0;  # Récupération du MFE 
			$longueur = @variable[1]; 
			$cg = 0;
		} 
		else
		{
			if (@variable[2]=~ /[CG]/) # on rencontre un 'g' ou un 'c' 
			{
				$cg++;
			}
			if (($longueur == @variable[1]) && ($nameSeq eq $seq)) # tester si la longueur actuelle est égale à la longueur de la séquence
			{
			
					
					open (RES, '>> '.$dirData.'/'.$seq.'/outMFEI.txt');
					print RES $nameSeq."\t".(($mfe/$longueur)*100)/(($cg/$longueur)*100)."\t".$mfe."\t".($mfe/$longueur)*100;
					system('chmod 777 '.$dirData.'/'.$seq.'/outMFEI.txt');
					close RES || die "Problème à la fermeture : $!";
			
				
			
			
			}
		}
	}
}


