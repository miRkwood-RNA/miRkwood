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

use FindBin;

BEGIN {
    use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
    use miRkwood;
    use miRkwood::Components;
}

my ( $CT, $boucleTermWithN_out ) = @ARGV;
miRkwood::Components::mask_CT_file( $CT, $boucleTermWithN_out );
