#!/usr/bin/perl -w
use Class::Struct;
use CGI; 
use Time::gmtime;

my $cgi = new CGI; 
$now = gmctime();

$now =~ s/[: ]//g;
$now = substr($now, 3);

$dirScript = "/var/www/arn/scripts/"; # chemin script
$dirData = "/var/www/arn/data/"; # chemin séquence de base
$dirImages = "/var/www/arn/images/"; # chemin images
$dirJob = $dirData."job".$now."/"; # chemin séquence de base
#mkdir($dirJob,0777);
system('mkdir '.$dirJob);



my $seqArea = $cgi->param('seqArea');
if ($seqArea eq "")  # cas upload fichier
{	
	$seq = "";
	my $upload = $cgi->upload('seqFile') || die "$!";
	open (INPUT, '>> '.$dirJob.'sequence.fas') || die "$!";
	while (my $ligne = <$upload>) 
	{
		$seq.= $ligne;
	}
	print INPUT lc($seq)."\n"; ##mise en minuscules
	close INPUT;
	system('chmod 777 '.$dirJob.'sequence.fas');
	$seq = lc($seq)."\n";
	system('sh '.$dirScript.'formatFasta.sh '.$dirJob.'sequence.fas '.$dirJob.'sequenceLoad.fas' ) ; # script qui elimine les retour a la ligne 
	if ($seq !~ /^( *>.+[\r\n]+([-\. atcgunwkmsydr0-9]+[\r\n]+)+){1,}$/)  
	{# erreur de syntaxe
		print $cgi->redirect('http://'.$ENV{SERVER_NAME}.'/cgi-bin/error.pl');       exit; 
	}
}
else #cas textArea
{
	open (INPUT, '>> '.$dirJob.'sequence.fas') || die "$!";
	print INPUT lc($seqArea);
	system ("sed -i 'N;s/\r//g' ".$dirJob."sequence.fas"); #commande sed permettant d'eliminer les caracteres du (copier-coller)
	system('sh '.$dirScript.'formatFasta.sh '.$dirJob.'sequence.fas '.$dirJob.'sequenceLoad.fas' ) ; # script qui elimine les retour a la ligne 
	while (my $ligne = <INPUT>) 
	{
		$seqArea.= $ligne;
	}
	$seqArea = lc($seqArea)."\n";#mise en miniscule
	
	if ($seqArea !~ /^( *>.+[\r\n]+([-\. atcgunwkmsydr0-9]+[\r\n]+)+){1,}$/)  
	{# erreur de syntaxe
		print $cgi->redirect('http://'.$ENV{SERVER_NAME}.'/cgi-bin/error.pl');       exit; 
	}
	close INPUT;
}

	
	open (LOAD, $dirJob.'sequenceLoad.fas') || die "$!";
	open (OUTPUT, '>> '.$dirJob.'sequenceUpload.fas') || die "$!";
	while (my $line = <LOAD>) 
	{
		if ($line=~/^>/)
		{
			$name = $line ; 
		}
		else
		{
			$line=~s/w/n/g;
			$line=~s/r/n/g;
			$line=~s/k/n/g;
			$line=~s/m/n/g;
			$line=~s/s/n/g;
			$line=~s/y/n/g;
			$line=~s/d/n/g;
			print OUTPUT $name.$line;
			$name = "";
		}
	
	}
	unlink ($dirJob.'sequenceLoad.fas');

# redirection vers la page wait en attendant le calcul
my $mail = $cgi->param('mail');

my $nameJob = $cgi->param('job');

print $cgi->redirect(-uri=>'http://'.$ENV{SERVER_NAME}.'/cgi-bin/wait.pl?jobId='.$now.'&mail='.$mail.'&nameJob='.$nameJob);



print $cgi->redirect(-uri=>'http://'.$ENV{SERVER_NAME}.'/cgi-bin/wait.pl?jobId='.$now.'&mail='.$mail);
print "Location: http://".$ENV{SERVER_NAME}."wait.pl?jobId=".$now."&mail=".$mail."\n\n";

my $check = $cgi->param('check');
my $mfei = $cgi->param('mfei');
my $randfold = $cgi->param('randfold');
my $SC = $cgi->param('selfContain');
my $align = $cgi->param('align');
$plant =  $cgi->param('db');
if ($mfei eq "") { $mfei = "notChecked" };
if ($randfold eq "") { $randfold = "notChecked" };
if ($SC eq "") { $SC = "notChecked" };
if ($align eq "") { $align = "notChecked" };
if ($check eq "") { $check = "notChecked" };
#execution de tous les scripts de traitements
system('perl '.$dirScript.'execute_scripts.pl '.$check.' '.$mfei.' '.$randfold.' '.$SC.' '. $align.' '.$dirJob.' '.$plant);

open (finish,'>', $dirJob.'finished') || die "$!"; 

