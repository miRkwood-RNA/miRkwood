package BinaryGenomeReader;

use strict;
use warnings;

use Inline CPP => <<'END_OF_CPP_CODE';

// A function that returns the 4 bases contained in the byte
SV* __getDnaQuadruplet(I32 byte) {
	char seq[4];
	char* s = seq+3;
	for (int i = 0; i < 4; i++, s--) {
		int r = byte & 3;
		*s = r == 0 ? 'T' : r == 1 ? 'C' : r == 2 ? 'A' : 'G';
		byte >>= 2;
	}
	return newSVpv(seq, 4);
}

END_OF_CPP_CODE

=method new

# Constructor

=cut
sub new {
    my ( $class, $genome_file ) = @_;
    my $this = bless {
        genome_file => $genome_file
    }, $class;
    $this->open_genome();
    return $this;
}

=method open_genome

# Opens the genome file. No input argument other than $self. This function is automatically called within the constructor. You don't have to call this function.

=cut
sub open_genome {
    my $this = shift;
	open (my $filehandle, '<', $this->{genome_file}) or die "Can't open '$this->{genome_file}': $!";
	binmode($filehandle);
	# Reading the header
	read($filehandle, my $signature, 4);
	read($filehandle, my $version, 4);
	read($filehandle, my $sequenceCount, 4);
	$sequenceCount = unpack('l', $sequenceCount);
	read($filehandle, my $reserved, 4);
	# Reading the index
	my %index = ();
	for (my $i = 0; $i < $sequenceCount; $i++) {
		read($filehandle, my $nameSize, 1);
		$nameSize = unpack('c', $nameSize);
		read($filehandle, my $name, $nameSize);
		read($filehandle, my $offset, 4);
		$index{$name} = {offset => unpack('l', $offset)};
	}

	foreach my $name (keys %index) {
		seek($filehandle, $index{$name}{offset}, 0);
		read($filehandle, my $dnaBaseCount, 4);
		read($filehandle, my $nBlockCount, 4);
		$dnaBaseCount = unpack('l', $dnaBaseCount);
		$nBlockCount = unpack('l', $nBlockCount);
		for (my $i = 0; $i < $nBlockCount; $i++) {
			read($filehandle, my $nBlockStarts, 4);
		}
		for (my $i = 0; $i < $nBlockCount; $i++) {
			read($filehandle, my $nBlockSizes, 4);
		}
		read($filehandle, my $maskBlockCount, 4);
		$maskBlockCount = unpack('l', $maskBlockCount);
		for (my $i = 0; $i < $maskBlockCount; $i++) {
			read($filehandle, my $maskBlockStarts, 4);
		}
		for (my $i = 0; $i < $maskBlockCount; $i++) {
			read($filehandle, my $maskBlockSizes, 4);
		}
		read($filehandle, $reserved, 4);
		$index{$name}{offset} = tell($filehandle); # We update the offset to the begining of the sequence
		$index{$name}{size} = $dnaBaseCount;
	}
	$this->{filehandle} = $filehandle;
	$this->{index} = \%index;
}


=method __getDnaBase

# A private function that reads a base.
# Input arguments:
# 	pos: An integer in the range 0 <= pos < 4 which determines which nucleotide to look for within the byte
# 	byte: The byte that codes for the nucleotides

=cut
sub __getDnaBase {
	my ($pos, $byte) = @_;
	my $s = 6 - 2*$pos;
	my $r = ((3 << $s) & $byte) >> $s;
	return $r == 0 ? 'T' : $r == 1 ? 'C' : $r == 2 ? 'A' : 'G';
}


=method set_current_chr

# A function that sets the current chromosome to read. You must then use read_chr to read the chr
# Input arguments:
# 	$self
# 	chr: The name of the chr to read

=cut
sub set_current_chr {
	my ($this, $chr) = @_;
	$this->{chr_offset} = $this->{index}{$chr}{offset};
}


=method read_sequence

# A function that reads a sequence from a given chr and position.
# Input arguments:
# 	$self
# 	chr: The name of the chr to read
# 	from: The begining position to read from. The first position of a chr starts at 1 (and not 0).
# 	to: The last position to read. The last position will be read and will be the last base of the sequence. Hence $self->read_sequence("chr", 1, 1); will return the first nucleotide and not an empty string.

=cut
sub read_sequence {
	my ($this, $chr, $from, $to) = @_;
	return __priv_readSequence($this->{filehandle}, $this->{index}{$chr}{offset}, $from-1, $to);
}


=method read_chr

# A function that reads a sequence from the current chr and a position.
# Input arguments:
# 	$self
# 	from: The begining position to read from. The first position of a chr starts at 1 (and not 0).
# 	to: The last position to read. The last position will be read and will be the last base of the sequence. Hence $self->read_chr(1, 1); will return the first nucleotide and not an empty string.

=cut
sub read_chr {
	my ($this, $from, $to) = @_;
	return __priv_readSequence($this->{filehandle}, $this->{chr_offset}, $from-1, $to);
}


=method get_chr_length

# A function that returns the number of nucleotide from a chr.
# Input arguments:
# 	$self
# 	chr: The name of the chr.

=cut
sub get_chr_length {
	my ($this, $chr) = @_;
	return $this->{index}{$chr}{size};
}


=method __priv_readSequence

# A private function that reads and return the sequence.
# Input arguments:
# 	filehandle: The genomic file handle
# 	chr_offet: A file offset that points to the begining of the chromosomal sequence
#	from: The begining position to read from. from must be 0 to read the first nucleotide.
#	to: The excluded end pos of the sequence

=cut
sub __priv_readSequence {
	my ($filehandle, $chr_offet, $from, $to) = @_;
	my $byte_begin_offset = $chr_offet + int($from/4);
	my $byte_end_offset = $chr_offet + int($to/4);
	my $bit_begin_offet = $from % 4;
	my $bit_end_offet = $to % 4;

	my $seq = '';
	
	seek($filehandle, $byte_begin_offset, 0);

	read($filehandle, my $byte, 1);
	$byte = unpack('c', $byte);
	if ($byte_begin_offset == $byte_end_offset) {
		for (my $i = $bit_begin_offet; $i < $bit_end_offet; $i++) {
			$seq .= __getDnaBase($i, $byte);
		}
	}
	else {
		# Reading the first byte
		for (my $i = $bit_begin_offet; $i < 4; $i++) {
			$seq .= __getDnaBase($i, $byte);
		}
		# Reading bytes
		for (my $i = $byte_begin_offset, my $e = $byte_end_offset-1; $i < $e; $i++) {
			read($filehandle, $byte, 1);
			$byte = unpack('c', $byte);
			$seq .= __getDnaQuadruplet($byte);
		}
		# Reading the last byte
		read($filehandle, $byte, 1);
		$byte = unpack('c', $byte);
		for (my $i = 0; $i < $bit_end_offet; $i++) {
			$seq .= __getDnaBase($i, $byte);
		}
	}
	return $seq;
}


sub close_genome {
	my $this = shift;
	close $this->{filehandle};
}


1;
