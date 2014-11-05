use strict;
use warnings;

my $filename = $ARGV[0];
my $threshold = $ARGV[1]; 
my $filename_out = $ARGV[2]; # given without the .bam suffix
if (!$filename_out) {
	die("No output filename given");
}

my $samtools_cmd = "samtools view $filename";
open( my $SAM, '-|', $samtools_cmd ) or die "Can't open '$samtools_cmd': $!";
my %read_ids = ();
while (my $line = <$SAM>) {
	chomp $line;
	my @fields = split( "\t", $line);
	if (defined $read_ids{$fields[0]}) {
		$read_ids{$fields[0]}++;
	}
	else {
		$read_ids{$fields[0]} = 1;
	}
}

$samtools_cmd = "samtools view -h $filename";
open(my $SAM2, '-|', $samtools_cmd ) or die "Can't open '$samtools_cmd': $!";
open(my $fOut, '>', $filename_out. '.sam');
while (my $line = <$SAM2>) {
	chomp $line;
	if ($line =~ /^@/) {
		print $fOut $line, "\n";
		next;
	}
	my @fields = split("\t", $line);
	if ($read_ids{$fields[0]} < $threshold) {
		print $fOut $line, "\n";
	}
	else {
		print "Read $fields[0] evicted.\n";
	}
}

$samtools_cmd = "samtools view -Sb $filename_out.sam > $filename_out.bam";
system($samtools_cmd);
system("rm $filename_out.sam");