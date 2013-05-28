#!/usr/bin/perl -w

my $rep = $ARGV[0];
print $rep, "\n";
opendir DIR, $rep;
my @fichiers;
@fichiers=readdir DIR;
closedir DIR;
foreach $fichier(@fichiers)
{
    print $fichier, "\n";
}
mkdir($rep."/hai");
