#!/usr/local/bin/perl
######################################################################
#  Lancement du script : perl energie.pl FichierFasta				 #
#  @author : Mohcen BENMOUNAH                                        #
#                                                                	 #
#                                                                    #
######################################################################

use strict;
use warnings;

use File::Spec;

my ( $CT, $dirData, $seq ) = @ARGV;

# ouverture du fichier de sortie de Unafold
open( my $FOut, '<', $CT ) || die "Unable to open: $!";

my $cg = 0;
my ( $nameSeq, $mfe, $longueur );
while ( my $line = <$FOut> ) {
	my @variable = split( '\s+', $line );

	if ( $variable[2] =~ /^ENERGY/ ) {
		$nameSeq = substr $variable[5],
		  0;    # Récupération du nom de la séquence
		$mfe      = substr $variable[4], 0;    # Récupération du MFE
		$longueur = $variable[1];
		$cg       = 0;
	}
	else {                                     # Not the ENERGY line
		if ( $variable[2] =~ /[CGcg]/ )        # on rencontre un 'g' ou un 'c'
		{
			$cg++;
		}
		if (
			( $longueur == $variable[1] )

			&& ( $nameSeq eq $seq )
		  ) # tester si la longueur actuelle est égale à la longueur de la séquence
		{
			my $MFEI_output =
			  File::Spec->catfile( $dirData, $seq, 'outMFEI.txt' );
			my $num = ( $mfe / $longueur ) * 100;
			my $other = $num / ( ( $cg / $longueur ) * 100 );
			open( my $RES, '>>', $MFEI_output ) || die "Unable to open: $!";
			print $RES $nameSeq . "\t" . $other . "\t" . $mfe . "\t" . $num;
			close $RES || die "Unable to close: $!";
			system( 'chmod 777 ' . $MFEI_output );
		}

	}    # /Else− not ENERGY line
}    # /While
close $FOut || die "Unable to close: $!";
