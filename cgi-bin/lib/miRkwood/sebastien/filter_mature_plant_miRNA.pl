use strict;
use warnings;

my $filename_miRNA = $ARGV[0];
my $filename_ids = $ARGV[1];
my $filename_out = $ARGV[2];
if (!$filename_out) {
	die("No output filename given");
}

open( my $FH_ids, '<', $filename_ids ) or die "Can't open '$filename_ids': $!";
my %read_ids = ();
while (my $line = <$FH_ids>) {
	chomp $line;
	if ($line ne '') {
		$read_ids{$line} = 1;
	}
}

open( my $FH_mature, '<', $filename_miRNA ) or die "Can't open '$filename_miRNA': $!";
open(my $fOut, '>', $filename_out);
my $evicted = 0;
while (my $line = <$FH_mature>) {
	chomp $line;
	if ($line =~ /^>/) {
		my @fields = split(' ', $line);
		my $name = '';
		if (scalar @fields == 5) {
			$name = $fields[2]. ' ' .$fields[3];
		}
		elsif (scalar @fields == 4) {
			$name = $fields[2];
		}
		if (defined $read_ids{$name}) {
			print $fOut $line, "\n";
			$line = <$FH_mature>;
			chomp $line;
			print $fOut $line, "\n";
		}
		else {
			$evicted++;
		}
	}
}

print "$evicted miRNA evicted.\n";