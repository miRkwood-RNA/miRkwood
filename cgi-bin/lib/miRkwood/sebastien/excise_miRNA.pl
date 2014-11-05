use strict;
use warnings;

my $filename = $ARGV[0];
my $filename_out = $ARGV[1];
if (!$filename_out) {
	die("No output filename given");
}

open(my $annotated_dna, '<', $filename) or die "Can't open '$filename': $!";
open(my $annotated_miRNA, '>', $filename_out) or die "Can't open '$filename_out': $!";
my $counter = 0;
while (my $line = <$annotated_dna>) {
	chomp $line;
	if ($line =~ /^#/) {
		print $annotated_miRNA $line. "\n";
		next;
	}
	my @fields = split("\t", $line);
	if ($fields[2] eq 'miRNA') {
		print $annotated_miRNA $line. "\n";
		$counter++;
	}
}
print "Done ! ($counter miRNA extracted).\n"