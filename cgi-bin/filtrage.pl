#!/usr/local/bin/perl
######################################################################
#  Lancement du script : perl filtrage.pl sequenceTest DBArapidopsis #
#  @author : Mohcen BENMOUNAH                                        #
#                                                                    #
#                                                                    #
######################################################################
 
print "Blast en cours .......\n";
my ($SEQ, $BD) = @ARGV;

system ('./blastx -query '.$SEQ.' -db '.$BD.' -outfmt 6 -max_target_seqs 1 -evalue 1E-5 -out OUTcgi.txt');
print "Fin du traitement Blast.\n";

print "Filtrage des CDS en cours .......\n";

@SeqNamesOUT=(); # tableau contenant la liste des séquences issues du Blast
@SeqNames=(); # tableau contenant la liste des séquences à traiter
@SeqDiff=(); # tableau contenant la liste des séquences non codantes
$i=0;

# Ouverture du fichier OUT.TXT et construction du tableau contenant la liste des séquences issues du blast
open(FOut, 'OUT.txt') || die "Problème à l\'ouverture : $!";
while(my $line = <FOut>)
{
	my @name = split('\t', $line);
	push(@SeqNamesOUT,@name[0]);
}
close FOut || die "Problème à la fermeture : $!";

# Ouverture du fichier SEQUENCE.TXT et construction du tableau contenant la liste des séquences issues du blast
open(FSeq, $SEQ) || die "Problème à l\'ouverture : $!";
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
open (RES, '>>SequencesFiltreesCDS.txt');
open(FSeq, $SEQ) || die "Problème à l\'ouverture : $!";
while(my $line = <FSeq>)
{
	if ( grep( /^>/,$line) )
	{
		$lineSeq = substr $line, 1,-1;
		if (@SeqDiff[$i] eq $lineSeq)
		{
			print @SeqDiff[$i];
			$record = true;		
			$i+;
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
close FSeq || die "Problème à la fermeture : $!";
close (RES);
print "\nFiltrage terminé et génération des séquences filtrées!!!.\n";
