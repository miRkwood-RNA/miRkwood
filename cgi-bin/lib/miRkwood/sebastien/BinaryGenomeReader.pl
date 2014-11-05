# ABSTRACT: BAM Pipeline object

use strict;
use warnings;

use Time::HiRes ();

use BinaryGenomeReader;
use Bio::DB::Fasta;

# my $genome = BinaryGenomeReader->new("../data/Athaliana_167.2bit");

my $db = Bio::DB::Fasta->new("../data/Athaliana_167.fa");

# foreach my $chr (sort keys %{$genome->{index}}) {
foreach my $chr ($db->ids) {
# 	$genome->set_current_chr($chr);
# 	print "$chr:99-119\t", $genome->read_chr(99, 119), "\n";
	print "$chr:99-119\t", $db->seq($chr, 99, 119), "\n";
}

# $genome->close_genome();