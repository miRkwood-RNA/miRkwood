#!/usr/bin/perl -w
use Class::Struct;
use CGI; 
use Time::gmtime;
use File::Spec;
use FindBin qw($Bin);
use Cwd qw( abs_path );
use File::Basename qw(dirname);

my $cgi = new CGI; 
$now = gmctime();

$now =~ s/[: ]//g;
$now = substr($now, 3);
my $dirJob_name = 'job'.$now;

my $local_dir = dirname( abs_path($0) );
my $rootdir = File::Spec->catdir( $local_dir, ".." );

$dirScript = File::Spec->catdir( $rootdir, 'scripts');
$dirData = File::Spec->catdir( $rootdir, 'data');
$dirLib = File::Spec->catdir( $rootdir, 'lib');

my $formatFasta_bin = File::Spec->catfile($dirScript, 'formatFasta.sh');

$dirJob = File::Spec->catdir( $dirData, $dirJob_name); # chemin séquence de base
mkdir $dirJob;

my $sequence = File::Spec->catfile($dirJob, 'sequence.fas');
my $sequence_load = File::Spec->catfile($dirJob, 'sequenceLoad.fas');
my $sequence_upload = File::Spec->catfile($dirJob, 'sequenceUpload.fas');

my $seqArea = $cgi->param('seqArea');
if ($seqArea eq "")  # cas upload fichier
{	
	$seq = "";
	my $upload = $cgi->upload('seqFile') || die "$!";
	open (INPUT, '>>', $sequence) || die "$!";
	while (my $ligne = <$upload>) 
	{
		$seq.= $ligne;
	}
	print INPUT lc($seq)."\n"; ##mise en minuscules
	close INPUT;
	chmod 777, $sequences;
	$seq = lc($seq)."\n";
	my $formatFasta_cmd = "sh $formatFasta_bin $sequence $sequence_load";
	system($formatFasta_cmd) ; # script qui elimine les retour a la ligne 
	if ($seq !~ /^( *>.+[\r\n]+([-\. atcgunwkmsydr0-9]+[\r\n]+)+){1,}$/)  
	{# erreur de syntaxe
		print $cgi->redirect('http://'.$ENV{SERVER_NAME}.'/cgi-bin/error.pl');       exit; 
	}
}
else #cas textArea
{
	open (INPUT, '>>', $sequences) || die "$!";
	print INPUT lc($seqArea);
	system ("sed -i 'N;s/\r//g' ".$sequences); #commande sed permettant d'eliminer les caracteres du (copier-coller)
	system("sh $formatFasta_bin $sequence $sequence_load" ) ; # script qui elimine les retour a la ligne 
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

	open (LOAD, $sequence_load) || die "$!";
	open (OUTPUT, '>>', $sequence_upload) || die "$!";
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
	unlink ($sequence_load);

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
my $perl_script = File::Spec->catfile($dirScript, 'execute_scripts.pl');
my $cmd = "perl -I$dirLib $perl_script $check $mfei $randfold $SC $align $dirJob $plant";
system($cmd);

my $finish_file = File::Spec->catfile($dirJob, 'finished');
open (my $finish, '>', $finish_file) || die "$!";
close $finish_file;

