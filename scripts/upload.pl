#!/usr/local/bin/perl -w
######################################################################
#  Lancement du script : perl filtrage.pl sequenceTest DBArapidopsis #
#  @author : Mohcen BENMOUNAH                                        #
#                                                                    #
#                                                                    #
######################################################################

use CGI; 
my $cgi = new CGI; 

 $dirScript = "/var/www/scripts/"; # chemin script
$dirImages = "/var/www/images/"; # chemin images
$dirData = "/var/www/data/"; # chemin séquence
$dirBlast = "/var/www/scripts/ncbi-blast-2.2.26+/bin/"; # chemin Blast
$dirRandfold = "/var/www/scripts/randfold-2.0/"; # chemin Randfold
my $upload = $cgi->upload("seq");
	open (INPUT, '>> '.$dirData.'sequenceUpload.fas') || die "$!";
	while (my $ligne = <$upload>) 
	{
		print INPUT $ligne;
	}
	close INPUT;
	system('chmod 777 '.$dirData.'sequenceUpload.fas');
