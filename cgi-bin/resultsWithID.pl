#!/usr/bin/perl -w
use Class::Struct;
use CGI; 
my $cgi = new CGI; 
use Cwd qw( abs_path );
use File::Basename qw(dirname);
use File::Spec;

my $local_dir = dirname( abs_path($0) );
my $rootdir = File::Spec->catdir( $local_dir, ".." );

$id_job = $cgi->param('run_id'); # récupération id job

$dirJob = abs_path(File::Spec->catdir( $rootdir, 'data', 'job'.$id_job)).'/';
#TODO Remove the trailing slash...

$names =[]; 
$pvalues =[]; 
$positions =[]; 
$mfeis =[]; 
$mfes =[];
$amfes =[];
$namesPositions =[]; 
$DNASequence =[];
$Vienna =[];  
$alignement =[];
$selfContain =[];
struct results => {        # déclaration de la structure de données 
	names  => '@',         
	pvalues => '@',
	positions => '@',
	mfeis => '@',
	mfes => '@',   
	amfes => '@',
	namesPositions=> '@',
	DNASequence=> '@',
	Vienna=> '@', 
	alignement=> '@',
	selfContain=> '@',
};
$name_job = $cgi->param('nameJob'); # récupération id job	

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
				
				#print $dirJob.$dir."/".$subDir."/".$subSubDir."\n";
				push(@$names,$dir); #récupération nom séquence
				@position = split(/__/,$subDir);
			
				push(@$positions,$position[1]); # récupération position
				push(@$namesPositions,$subDir); # récupération noms + positions
				$MfeiExist = false;
				$RandfoldExist = false;
				$alignExist = false;
				$SCExist = false;
				if( -e $dirJob.$dir."/".$subDir.'/pvalue.txt' ) # si fichier existe
				{
					$RandfoldExist = true;
					open (PVALUE, $dirJob.$dir."/".$subDir.'/pvalue.txt') || die "$!";
					while (my $line = <PVALUE>) 
					{
						if ( $line =~/(.*)\t(.*)\t(.*)/ )
						{    
							push(@$pvalues,$3); # récupération pvalues
						}
					}
					close PVALUE;
					
				}	
				#Récupération valeur MFEI
				if( -e $dirJob.$dir."/".$subDir.'/outMFEI.txt' ) # si fichier existe
				{
					$MfeiExist = true;
					open (MFEI, $dirJob.$dir."/".$subDir.'/outMFEI.txt') || die "$!";
					while (my $line = <MFEI>) 
					{
						if ( $line =~/(.*)\t(.*)\t(.*)\t(.*)/ )
						{    
							push(@$mfeis,$2); # récupération mfei
							push(@$mfes,$3); # récupération mfei
							push(@$amfes,$4); # récupération mfei
						}
					}
					close MFEI;
				}
				#Récupération valeur self contain
				if( -e $dirJob.$dir."/".$subDir.'/selfContain.txt' ) # si fichier existe
				{
					$SCExist = true;
					open (SC, $dirJob.$dir."/".$subDir.'/selfContain.txt') || die "$!";
					while (my $line = <SC>) 
					{
						if ( $line =~/(.*) (.*)/ )
						{    
							push(@$selfContain,$2); # récupération mfei
						}
					}
					close SC;
				}
					
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
				#Récupération alignement avec mirBase 
				if( -e $dirJob.$dir."/".$subDir.'/alignement.txt' ) # si fichier existe
				{
					$alignExist = true;
					open (align, $dirJob.$dir."/".$subDir.'/alignement.txt') || die "$!";
					$align = "none";
					while (my $line = <align>) 
					{
						if ( $line =~/^C4/ )
						{    
							$align = "/arn/data/job".$id_job.'/'.$dir."/".$subDir.'/alignement.txt';
							last;
						}
					
					}
					push(@$alignement,$align); 
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
$myResults->mfes($mfes);
$myResults->amfes($amfes);
$myResults->positions($positions);
$myResults->DNASequence($DNASequence);
$myResults->Vienna($Vienna);
$myResults->alignement($alignement);
$myResults->selfContain($selfContain);
print <<DATA;
Content-type: application/xhtml+xml

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

	<head>
		<link rel="stylesheet" type="text/css" href="/arn/css/script.css" />
		<script type="text/javascript" language="Javascript" src="/arn/js/results.js"> </script>
		<script type="text/javascript" src="/arn/js/graphics.js"></script>
	<script type="text/javascript" src="/arn/js/miARN.js"></script>	
	</head>
	<body onload='main();'>
		<div class="titreDiv"> Identification of miRNA/miRNA hairpins results:</div>
DATA

if ($name_job ne "")
{
	print "<div class='titleJob' ><li>Title Job : ".$name_job."</li></div>";
	
}
print <<DATA;
		<div id="table" ></div>
		<div id="singleCell"> </div>
		<results id="all">
DATA

for ($i=0;$i<scalar(@$names); $i++)
{
	$baliseSeq="<Sequence name='";
	$baliseSeq=$baliseSeq.$myResults->names($i);
	$baliseSeq=$baliseSeq."' position='".$myResults->positions($i);
	if ($MfeiExist eq "true")
	{
		$baliseSeq=$baliseSeq."' mfei='".$myResults->mfeis($i)."' mfe='".$myResults->mfes($i)."' amfe='".$myResults->amfes($i);
	}
	if ($RandfoldExist eq "true")
	{
		$baliseSeq=$baliseSeq."'  p_value='".$myResults->pvalues($i);
	}
	if ($SCExist eq "true")
	{
		$baliseSeq=$baliseSeq."'  self_contain='".$myResults->selfContain($i);
	}
	if ($alignExist eq "true")
	{
		$baliseSeq=$baliseSeq."'  alignment='".$myResults->alignement($i);
	}
	
		
	$baliseSeq=$baliseSeq."' image='"."/arn/data/job".$id_job.'/'.$myResults->names($i)."/".$myResults->names($i)."__".$myResults->positions($i)."/image.png' Vienna='".$myResults->Vienna($i)."' DNASequence='".$myResults->DNASequence($i)."'></Sequence>\n";
	print $baliseSeq;
	#print "<Sequence name='".$myResults->names($i)."' position='".$myResults->positions($i)."' mfei='".$myResults->mfeis($i)."'  p_value='".$myResults->pvalues($i)."' image='".$myResults->names($i)."_ss.ps' Vienna='".$myResults->Vienna($i)."' DNASequence='".$myResults->DNASequence($i)."'></Sequence>\n";	 
	#print "<Sequence name='".$myResults->names($i)."' position='".$myResults->positions($i)."' mfei='".$myResults->mfeis($i)."'  p_value='".$myResults->pvalues($i)."' image='".$dirImages.$myResults->names($i)."_ss.ps' ' DNA='".$myResults->DNASequence($i)."' Vienna='".$myResults->Vienna($i)."></Sequence>\n";	
}
print <<DATA;		
		</results>
		
	</body>

</html>

DATA

