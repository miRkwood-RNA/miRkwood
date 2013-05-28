#!/usr/local/bin/perl
######################################################################
#  Lancement du script : perl filtrage.pl sequenceTest DBArapidopsis #
#  @author : Mohcen BENMOUNAH                                        #
#                                                                    #
#                                                                    #
######################################################################
@SeqNamesCT=(); # tableau contenant la liste des séquences du fichier CT qui reste après passage du script tigeBoucle.pl
 
my ( $CT,$MultiFas,$dirData) = @ARGV;
open(FOut,  $CT) || die "Problème à l\'ouverture : $!"; 
while(my $line = <FOut>)
{	
	my @variable = split('\s+', $line);
	if (@variable[2]=~ /^ENERGY/)
	{	
		$nomSeq = substr @variable[5], 0; # Récupération du nom de la séquence
		push(@SeqNamesCT,$nomSeq);
	} 
	
}

## Création du fichier FASTA filtrée des séquences tiges boucle
$record = false;
$i = 0;
open (RES, '>>'.$dirData.'tigeBouclesSequence.txt');
open(FSeq, $MultiFas) || die "Problème à l\'ouverture : $!";
while(my $line = <FSeq>)
{
	if ( $line =~/^>\s(.*) \(/ )
	{
		$lineSeq = $1;   
		
		if (@SeqNamesCT[$i] eq $lineSeq)
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
			printf RES $line;
		}
}
close FSeq || die "Problème à la fermeture : $!";
close (RES);
system('chmod 777 '.$dirData.'tigeBouclesSequence.txt');
