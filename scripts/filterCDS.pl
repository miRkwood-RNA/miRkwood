#!/usr/local/bin/perl
######################################################################
#  Lancement du script : perl filtrage.pl sequenceTest DBArapidopsis #
#  @author : Mohcen BENMOUNAH                                        #
#                                                                    #
#                                                                    #
######################################################################
 
my ($dirBlast, $dirData, $dirJob,$plant) = @ARGV;
system($dirBlast.'blastx -query '.$dirJob.'sequenceUpload.fas -db '.$dirData.$plant.'.fas -outfmt 6 -max_target_seqs 1 -evalue 1E-5 -out '.$dirJob.'outBlast.txt');
system('chmod 777 '.$dirJob.'outBlast.txt');
##Filtrage du ficher : Elimination des régions codantes###
@SeqNamesOUT=(); # tableau contenant la liste des séquences issues du Blast
@SeqNames=(); # tableau contenant la liste des séquences à traiter
@SeqDiff=(); # tableau contenant la liste des séquences non codantes
$i=0;

# Ouverture du fichier outBlast.TXT et construction du tableau contenant la liste des séquences issues du blast
open(FOut, $dirJob.'outBlast.txt') || die "Problème à l\'ouverture : $!";
while(my $line = <FOut>)
{
	my @name = split('\t', $line);
	push(@SeqNamesOUT,@name[0]);
}
close FOut || die "Problème à la fermeture : $!";

# Ouverture du fichier sequenceUpload.fas et construction du tableau contenant la liste des séquences issues du blast
open(FSeq, $dirJob.'sequenceUpload.fas') || die "Problème à l\'ouverture : $!";
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
open (RES, '>>'.$dirJob.'Sequences.fas');
open(FSeq, $dirJob.'sequenceUpload.fas') || die "Problème à l\'ouverture : $!";
while(my $line = <FSeq>)
{
	if ( grep( /^>/,$line) )
	{
		$lineSeq = substr $line, 1,-1;
		if (@SeqDiff[$i] eq $lineSeq)
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
system('chmod 777 '.$dirJob.'Sequences.fas');

