# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use KMean;

my $classifier = KMean->new();

my @values = (4, 6, 0, 5);

# my @array = ();
# foreach my $v (@values) {
# 	KMean::__binary_insert(\@array, $v);
# 	print "Printing array:\n";
# 	foreach my $v (@array) {
# 		print "\t$v\n";
# 	}
# }
# 
# foreach my $v (@array) {
# 	print "$v\n";
# }

foreach my $v (@values) {
	print "Added $v: ", $classifier->add_point($v), "\n";
}