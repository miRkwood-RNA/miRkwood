#!/usr/local/bin/perl
######################################################################
# Script permettant de détecter les boucles terminales               #
#                                                                    #
#	Lancement du script : perl boucleTerm.pl Fichier.txt.ct  	     #
#  @author : Mohcen BENMOUNAH                                        #                                                                	 #										         #
#                                                                    #
######################################################################
use strict;
use warnings;

use PipelineMiRNA::Components;


my ( $CT, $boucleTermWithN_out ) = @ARGV;
PipelineMiRNA::Components::mask_CT_file( $CT, $boucleTermWithN_out );

