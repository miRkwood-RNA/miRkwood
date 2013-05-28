#!/usr/bin/perl -w

use Class::Struct;

$dirJob = "/home/mohcen/arn/data/job/"; # chemin séquence de base
$names =[]; 
$pvalues =[]; 
$positions =[]; 
$mfeis =[]; 
$namesPositions =[]; 
$DNASequence =[];
$Vienna =[];  
struct results => {        # déclaration de la structure de données 
    names  => '@',         
    pvalues => '@',
    positions => '@',
    mfeis => '@',
	namesPositions=> '@',
	DNASequence=> '@',
	Vienna=> '@', 
};
	
	
opendir DIR, $dirJob; #ouverture répertoire job
my @dirs;
@dirs=readdir DIR;
closedir DIR;
foreach $dir(@dirs) # parcours du contenu
{
	if ($dir ne "." && $dir ne ".." && -d $dirJob.$dir) #si fichier est un répertoire
	{
		#print $dir."\n";
		opendir DIR, $dirJob.$dir; # ouverture du sous répertoire 
		my @files;
		@files=readdir DIR;
		closedir DIR;
		foreach $subDir(@files)
		{
			if (($subDir ne ".") && ($subDir ne "..") && -d $dirJob.$dir."/".$subDir) # si le fichier est de type repertoire
			{
				print $dirJob.$dir."/".$subDir.'/pvalue.txt'."\n";
				#print $dirJob.$dir."/".$subDir."/".$subSubDir."\n";
				push(@$names,$dir); #récupération nom séquence
				@position = split(/_/,$subDir);
			
				push(@$positions,$position[1]); # récupération position
				push(@$namesPositions,$subDir); # récupération noms + positions
				#print $subSubDir."\n";
				open (PVALUE, $dirJob.$dir."/".$subDir.'/pvalue.txt') || die "$!";
				while (my $line = <PVALUE>) 
				{
					if ( $line =~/(.*)\t(.*)\t(.*)/ )
					{    
						push(@$pvalues,$3); # récupération pvalues
					}
				}
				close PVALUE;
			
				open (MFEI, $dirJob.$dir."/".$subDir.'/outMFEI.txt') || die "$!";
				while (my $line = <MFEI>) 
				{
					if ( $line =~/(.*)\t(.*)/ )
					{    
						push(@$mfeis,$2); # récupération mfei
					}
				}
				close PVALUE;
				#Récupération séquence et format Vienna
		
				open (Vienna, $dirJob.$dir."/".$subDir.'/outViennaTraited.txt') || die "$!";
				while (my $line = <Vienna>) 
				{
					if ( $line =~/(.*)\t(.*)\t(.*)/ )
					{    
						push(@$DNASequence,$2); # récupération sequence
						push(@$Vienna,$3); # récupération Vienna
					
					}
				}
				
			}
		}
	}
}
	
	
$myResults=results->new(); # initialisation du struct 
$myResults->names($names);
$myResults->pvalues($pvalues);
$myResults->namesPositions($namesPositions);
$myResults->mfeis($mfeis);
$myResults->positions($positions);
$myResults->DNASequence($DNASequence);
$myResults->Vienna($Vienna);	
	
	
	
	
for ($i=0;$i<scalar(@$names); $i++)
{
	#print "<Sequence name='".$myResults->names($i)."' position='".$myResults->positions($i)."' mfei='".$myResults->mfeis($i)."'  p_value='".$myResults->pvalues($i)."' image='".$myResults->names($i)."_ss.ps' Vienna='".$myResults->Vienna($i)."' DNASequence='".$myResults->DNASequence($i)."'></Sequence>\n";	 

#print "<Sequence name='".$myResults->names($i)."' position='".$myResults->positions($i)."' mfei='".$myResults->mfeis($i)."'  p_value='".$myResults->pvalues($i)."' image='".$dirImages.$myResults->names($i)."_ss.ps' ' DNA='".$myResults->DNASequence($i)."' Vienna='".$myResults->Vienna($i)."></Sequence>\n";	

}

