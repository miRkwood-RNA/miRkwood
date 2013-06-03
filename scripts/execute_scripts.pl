#!/usr/bin/perl -w

use warnings;
use strict;

use File::Path 'rmtree';
use File::Basename;
use Cwd qw( abs_path );

my $local_dir = dirname(abs_path($0));
my $rootdir = File::Spec->catdir($local_dir, "..");

my $dirScript = File::Spec->catdir( $rootdir, 'scripts' );  # chemin script
my $dirImages =  File::Spec->catdir( $rootdir, 'images');  # chemin images
my $dirData =  File::Spec->catdir( $rootdir, 'data');  # chemin séquence
my $dirBlast =  File::Spec->catdir( $dirScript, 'ncbi-blast-2.2.28+-src', 'c++', 'GCC460-Debug', 'bin'); # chemin Blast

## Programs ##
my $rnafold_bin = File::Spec->catfile($dirScript, 'ViennaRNA-1.8.5', 'Progs', 'RNAfold');
my $randfold_bin = File::Spec->catfile($dirScript, 'randfold-2.0', 'randfold');
my $b2ct_bin = File::Spec->catfile($dirScript, 'ViennaRNA-1.8.5', 'Utils', 'b2ct');
my $selfcontain_bin = File::Spec->catfile($dirScript, 'selfcontain_unix', 'selfcontain.py');
my $exonerate_bin = File::Spec->catfile($dirScript, 'exonerate-2.2.0-i386', 'bin', 'exonerate');
my $varna_bin = File::Spec->catfile($dirScript, 'VARNAv3-9.jar');

## Data ##
my $mirbase_file = File::Spec->catfile($dirData, 'MirbaseFile.txt');
my $matrix_file = File::Spec->catfile($dirData, 'matrix');


my ($check,$mfei,$randfold, $SC,$align,$dirJob,$plant) = @ARGV;

if ($check eq "checked")
{
	#appel script filterCDS.pl qui permet de filter les CDS
	my $filter_script = File::Spec->catfile($dirScript, 'filterCDS.pl');
	system("perl $filter_script $dirBlast $dirData $dirJob $plant");
}
else
{
	system('mv '.$dirJob.'sequenceUpload.fas '.$dirJob.'Sequences.fas');
}

##Passage du multifasta -> fasta et appel au script Stemloop
my %tab=();
my $nameSeq;
open my $ENTREE_FH, '<', File::Spec->catfile($dirJob, 'Sequences.fas') or die "Impossible d'ouvrir le fichier d'entree  : $!";
while  (my $line=<$ENTREE_FH>)
{	
	if ( grep {/^>/} $line )
	{
		$nameSeq = substr $line,0;	
	}else{
		$tab{$nameSeq}= $line;
	}
}
close $ENTREE_FH;	
foreach my $name (keys %tab)
{ 	
	my $temp_file = File::Spec->catfile($dirJob, 'tempFile.txt');
	open (my $TEMPFILE_FH, '>', $temp_file) or die "Impossible d'ouvrir le fichier d'entree  : $!";
	system("chmod 777 $temp_file");
	print $TEMPFILE_FH $name.$tab{$name};
	my $name = substr $name, 1,-1 ;

	my $sequence_dir = File::Spec->catdir($dirJob, $name);
	system("mkdir $sequence_dir");
		
	my $rnafold_output = File::Spec->catfile($sequence_dir, 'fichierRNAfold.txt');
	system("$rnafold_bin -noPS < $temp_file > $rnafold_output");
	system("chmod 777 $rnafold_output");
	####conversion en format CT 
	my $ct_file = File::Spec->catfile($dirJob, $name, 'fichierOutB2ct.ct');
	system("$b2ct_bin < $rnafold_output > $ct_file");
	system("chmod 777 $ct_file");
	my $tigeboucle_script = File::Spec->catfile($dirScript, 'Tigeboucle3bc.pl');
	system("perl $tigeboucle_script $ct_file $sequence_dir");
	print $TEMPFILE_FH "perl $tigeboucle_script $ct_file $sequence_dir";
	unlink $temp_file;
	
}

##Traitement fichier de sortie outStemloop
opendir DIR, $dirJob; #ouverture répertoire job
my @dirs;
@dirs=readdir DIR;
closedir DIR;
foreach my $dir(@dirs) # parcours du contenu
{
	if ($dir ne "." && $dir ne ".." && -d $dirJob.$dir) #si fichier est un répertoire
	{

		my $sequence_dir = File::Spec->catdir($dirJob, $dir);
		opendir DIR, $sequence_dir; # ouverture du sous répertoire
		my @files;
		@files=readdir DIR;
		closedir DIR;
		foreach my $file(@files)
		{
	
				if ($file ne "." && $file ne ".." && -d File::Spec->catdir($sequence_dir, $file)) # si le fichier est de type repertoire
				{
					####Traitement fichier de sortie outStemloop
					my $candidate_dir = File::Spec->catdir($sequence_dir, $file);
					chmod 0777, $candidate_dir;
					#system('perl /var/www/arn/scripts/image.pl '.$dirScript.' '.$dirJob.$dir."/".$file.'/');
					my $candidate_rnafold_out = File::Spec->catfile($candidate_dir, 'outRNAFold.txt');
					my $seq_file = File::Spec->catfile($candidate_dir, 'seq.txt');
					system("$rnafold_bin -noPS < $seq_file > $candidate_rnafold_out");
					system("chmod 777 $candidate_rnafold_out");

					####conversion en format CT
					my $candidate_ct_file = File::Spec->catfile($candidate_dir, 'outB2ct.ct');
					system("$b2ct_bin < $candidate_rnafold_out > $candidate_ct_file");
					system("chmod 777 $candidate_ct_file");
					
					#system('./test.sh '.$dirScript.' '.$dirJob.$dir.'/'.$file.'/ '.' > trash 2> errors');
					my $varna_image = File::Spec->catfile($candidate_dir, 'image.png');
					system("/usr/bin/java -cp $varna_bin fr.orsay.lri.varna.applications.VARNAcmd -i $candidate_ct_file -o $varna_image > trash");

					## traitement du fichier OutVienna pour la récupération des données(Format Vienna, séquence ADN)
					my $out_Vienna = File::Spec->catfile($candidate_dir, 'outViennaTraited.txt');
					open (my $TRAITED_FH, '>', $out_Vienna) || die "$!";
					open (my $INPUT_FH, '<', $candidate_rnafold_out) || die "$!";
					my ($nameSeq, $dna, $Vienna);
					while (my $line = <$INPUT_FH>)
					{
						if (($line=~/^>(.*)/)){  # nom sequence
							$nameSeq = $1;
						}
						elsif (($line=~/^[a-zA-Z]/)){  # récupération de la sequence adn
							$dna= substr $line,0,-1;
						}
						elsif (($line=~/(.*) /))        
						{ 
							$Vienna= $1;
							print $TRAITED_FH $nameSeq."\t".$dna."\t".$Vienna."\n";#récupération du format Vienna
						}
					 }
					close $INPUT_FH;
					close $TRAITED_FH;
					system("chmod 777 $out_Vienna");
					####calcul MFEI (appel script energie.pl)
					if ($mfei eq "mfeiChecked") 
					{
						my $energie_script = File::Spec->catfile($dirScript, 'energie.pl');
						system("perl $energie_script $candidate_ct_file $sequence_dir $file");
					}
					####calcul p-value randfold  
					if ($randfold eq "randfoldChecked") 
					{
						my $randfold_out = File::Spec->catfile($candidate_dir, 'pvalue.txt');
						system("$randfold_bin -d $seq_file 7 > $randfold_out");
						system("chmod 777 $randfold_out");
					}  
					####calcul self-contain   
					if ($SC eq "SCChecked") 
					{
						my $selfcontain_out = File::Spec->catfile($candidate_dir, 'selfContain.txt');
						system("python $selfcontain_bin -i $seq_file -n 100  > $selfcontain_out");
						system("chmod 777 selfcontain_out");
					}
					####### creation sequence boucle terminale masquee avec des N pour chaque sequence (repertoire ) et resultat alignement mirBASE
					if ($align eq "alignChecked")
					{
						my $seqN = File::Spec->catfile($candidate_dir, 'seqWithN.txt');
						open (my $SEQN_FH, '>>', $seqN) or die "Impossible d'ouvrir le fichier d'entree  : $!";
						open (my $SEQUENCE_FH, '<', $seq_file) or die "Impossible d'ouvrir le fichier d'entree  : $!";
						while  (my $line=<$SEQUENCE_FH>)
						{
							if ($line=~/^>/)
							{
								print $SEQN_FH $line;
							}
							else 
							{
								$line =~ s/U/T/g;
								print $SEQN_FH $line;
							}
						}
						close $SEQN_FH;
						close $SEQUENCE_FH;
						my $exonerate_out = File::Spec->catfile($candidate_dir, 'alignement.txt');
						system("$exonerate_bin -E --model affine:bestfit $mirbase_file $seqN -d $matrix_file --bestn 1 --score -3 -e -1 -o -1 > $exonerate_out");
					}
				}
		}
	}
}



