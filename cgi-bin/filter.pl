#!/usr/local/bin/perl
######################################################################
#  Lancement du script : perl filtrage.pl sequenceTest DBArapidopsis #
#  @author : Mohcen BENMOUNAH                                        #
#                                                                    #
#                                                                    #
######################################################################
 
my ($dirBlast, $dirData,$plant) = @ARGV;
print $dirBlast;
system($dirBlast.'blastx -query '.$dirData.'sequence.fas -db '.$dirData.$plant.'.fas -outfmt 6 -max_target_seqs 1 -evalue 1E-5 -out '.$dirData.'outBlast.txt');
system('chmod 777 '.$dirData.'outBlast.txt');
##Filtrage du ficher : Elimination des régions codantes###
