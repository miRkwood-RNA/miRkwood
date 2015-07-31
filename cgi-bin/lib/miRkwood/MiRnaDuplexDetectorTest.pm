package MiRnaDuplexDetectorTest;

use strict;
use warnings;
use Inline (Config => DIRECTORY => '/tmp/',);
use Inline CPP => <<'END_OF_C_CODE';
// Here starts C++ code

#include <vector>
#include <utility>
#include <cmath>


// This function returns true iff the nucleotide 'a' can bind with 'b'
// If either nucleotide is 'N', it returns false
// 'U' and 'T' are treated equally
// You may call this function from Perl. Simply write '$result = equiv_forward_strand("C", "T");' and it will work.
bool equiv_forward_strand(char a, char b) {
	if (a == 'G')
		return b == 'C' || b == 'U' || b == 'T';
	else if (a == 'A')
		return b == 'U' || b == 'T';
	else if (a == 'C')
		return b == 'G';
	else if (a == 'U' || a == 'T')
		return b == 'A' || b == 'G';
	return false; // a == 'N'
}

// This function returns true iff the nucleotide 'a' can bind with 'b' on the reverse strand. (i.e. if the complement of 'a' can bind to the complement of 'b')
// If either nucleotide is 'N', it returns false
// 'U' and 'T' are treated equally
// You may call this function from Perl. Simply write '$result = equiv_reverse_strand("C", "T");' and it will work.
bool equiv_reverse_strand(char a, char b) {
	if (a == 'A')
		return b == 'U' || b == 'T' || b == 'C';
	else if (a == 'C')
		return b == 'G' || b == 'A';
	else if (a == 'G')
		return b == 'C';
	else if (a == 'U' || a == 'T')
		return b == 'A';
	return false; // a == 'N'
}

// This function just create a perl array reference. It returns: [r, s] (a reference on (r, s)).
// SV stands for Scalar value
// You shouldn't call this function from Perl
SV* __get_pos(int r, int s) {
	AV* array; // AV stands for Array Value
	array = newAV();
	av_push(array, newSViv(r)); // iv stands for Integer value
	av_push(array, newSViv(s));
	return newRV_noinc((SV*)array); // rv stands for Reference value
}

struct EditDistanceScore {
		int value, from_j;

		EditDistanceScore(int v = 0, int j = 0) : value(v), from_j(j) {}

		EditDistanceScore& increased() { value++; return *this; }
		EditDistanceScore increasedBy(int val) const { EditDistanceScore score(*this); score.value+= val; return score; }
};

class matrix {
	public:
		//typedef EditDistanceScore value_type;

	private:
		std::size_t myRows, myColumns;
		std::vector<EditDistanceScore> myData;

	public:
		matrix(std::size_t rows, std::size_t columns) : myRows(rows), myColumns(columns), myData(rows*columns) {}
		matrix(std::size_t rows, std::size_t columns, EditDistanceScore const& value) : myRows(rows), myColumns(columns), myData(rows*columns, value) {}

		std::size_t rows() const { return myRows; }
		std::size_t columns() const { return myColumns; }

		std::size_t size() const { return myData.size(); }

		//EditDistanceScore* data() { return myData.data(); }

		EditDistanceScore* scanRow(std::size_t row) { return myData.data() + myColumns*row; }
		EditDistanceScore const* scanRow(std::size_t row) const { return myData.data() + myColumns*row; }

		EditDistanceScore& operator()(std::size_t i_row, std::size_t j_column) { return myData[myColumns*i_row+j_column]; }
		EditDistanceScore const& operator()(std::size_t i_row, std::size_t j_column) const { return myData[myColumns*i_row+j_column]; }
};

// This function just create a perl array reference. It returns: [r, s] (a reference on (r, s)).
// SV stands for Scalar value
// You shouldn't call this function from Perl
SV* __get_region(int begin_read, int end_read, int begin_seq, int end_seq) {
	AV* array; // AV stands for Array Value
	array = newAV();
	av_push(array, __get_pos(begin_read, end_read)); // iv stands for Integer value
	av_push(array, __get_pos(begin_seq, end_seq));
	return newRV_noinc((SV*)array); // rv stands for Reference value
}

template <bool (*_equiv_function)(char, char)>
EditDistanceScore process_score(matrix& score_matrix, std::size_t i, std::size_t j, char read, char sequence, bool allowBinding) {
	EditDistanceScore s;
	if (_equiv_function(read, sequence) && allowBinding)
		s = (score_matrix(i, j) = score_matrix(i-1, j-1));
	else {
		EditDistanceScore subs = score_matrix(i-1,j-1).increasedBy(1);
		EditDistanceScore insertion = score_matrix(i-1,j).increasedBy(2);
		EditDistanceScore deletion = score_matrix(i,j-1).increasedBy(2);
		if (subs.value <= insertion.value and subs.value <= deletion.value)
			s = (score_matrix(i, j) = subs);
		else if (insertion.value <= subs.value and insertion.value <= deletion.value)
			s = (score_matrix(i, j) = insertion);
		else
			s = (score_matrix(i, j) = deletion);
	}
	return s;
}

// EDIT DISTANCE EXTENDED

#include <limits>

struct Locus {
		unsigned int begin, end;
		Locus(unsigned int b, unsigned int e) : begin(b), end(e) {}

		SV* toPerl() const { return __get_pos(begin, end); }
};

enum ExtensionResult {
	DistanceAboveThreshold, TestPassed, TooShortToTest, ReadAndSeqOverlap
};

#include <vector>

SV* vectorToPerl(std::vector<Locus> const& loci) {
	AV* results = newAV();
	for (std::vector<Locus>::const_iterator it = loci.begin(), e = loci.end(); it != e; ++it) {
		av_push(results, it->toPerl());
	}
	return newRV_noinc((SV*)results);
}


typedef int TextBounds;

class MiRnaDetector;

template <bool (*_equiv_function)(char, char)>
void computeDistance(char* text, Locus read, Locus sequence, TextBounds const& bounds, MiRnaDetector& self);
template <bool (*_equiv_function)(char, char)>
ExtensionResult checkInwardExtent(char* text, Locus read, Locus sequence, Locus locus, TextBounds const& bounds, MiRnaDetector& self);
template <bool (*_equiv_function)(char, char)>
ExtensionResult checkOutwardExtent(char* text, Locus read, Locus sequence, Locus locus, TextBounds const& bounds, MiRnaDetector& self);
template <bool (*_equiv_function)(char, char)>
std::vector<Locus> tp_detect_on_strand(char* text, int text_size, int read_beg, int read_end, int seq_beg, int seq_end, MiRnaDetector& self);

class MiRnaDetector {

	public:

		matrix m_matrix;
		int m_miRnaMinErrorsThreshold, m_miRnaMaxErrorsThreshold;
		int m_miRnaMinLengthToExtend;
		int m_increaseMiRnaEachNt;
		double m_miRnaExtendedRatioErrorThreshold;

		void initRegion(int from_i, int to_i, int from_j, int to_j, bool realDistance) {
			for (int i = 0, e = to_i-from_i; i <= e; i++) {
				EditDistanceScore* data = m_matrix.scanRow(i+from_i) + from_j;
				*data = EditDistanceScore(i, 0);
			}
			if (realDistance) {
				EditDistanceScore* data = m_matrix.scanRow(from_i) + from_j+1;
				for (int j = 1, e = to_j-from_j; j <= e; j++, data++)
					*data = EditDistanceScore(j, j);
			}
			else {
				EditDistanceScore* data = m_matrix.scanRow(from_i) + from_j+1;
				for (int j = 1, e = to_j-from_j; j <= e; j++, data++)
					*data = EditDistanceScore(0, j);
			}
		}

		std::vector<Locus> lookForMiRnaResult(int to_i, int from_j, int to_j, bool& meetMinThreshold) {
			meetMinThreshold = false;
			EditDistanceScore* data = m_matrix.scanRow(to_i) + from_j+1;
			std::vector<Locus> result;
			for (int j = 1, e = to_j-from_j; j <= e; j++, data++) {
				if (data->value <= m_miRnaMinErrorsThreshold) {
					if (!meetMinThreshold)
						result.clear();
					meetMinThreshold = true;
					result.push_back(Locus(e-j, e-data->from_j));
				}
				else if (!meetMinThreshold && data->value <= m_miRnaMaxErrorsThreshold)
					result.push_back(Locus(e-j, e-data->from_j));
			}
			return result;
		}
		bool lookForResult(Locus read, Locus sequence) const {
			// Bottom side
			for (uint i = read.begin+m_miRnaMinLengthToExtend; i <= read.end; i++) {
				EditDistanceScore const* data = m_matrix.scanRow(i) + sequence.begin+1;
				for (uint j = 1, e = sequence.end-sequence.begin; j <= e; j++, data++) {
					if (data->value <= m_miRnaExtendedRatioErrorThreshold*std::max(i-read.begin, j))
						return true;
				}
			}
			// Top right
			for (uint i = read.begin+1, e = read.begin+1+m_miRnaMinLengthToExtend; i <= e; i++) {
				EditDistanceScore const* data = m_matrix.scanRow(i) + sequence.begin+m_miRnaMinLengthToExtend;
				for (uint j = sequence.begin+m_miRnaMinLengthToExtend; j <= sequence.end; j++, data++) {
					if (data->value <= m_miRnaExtendedRatioErrorThreshold*std::max(i-read.begin, j-sequence.begin))
						return true;
				}
			}
			return false;
		}

	public:
		MiRnaDetector(int textSize) : m_matrix(textSize+1, textSize+1, EditDistanceScore(0, 0)) {
			setMiRnaMinErrorsThreshold(4);
			setMiRnaMaxErrorsThreshold(7);
			setMiRnaMinLengthToExtend(15);
			setIncreaseMiRnaEachNt(5);
			setMiRnaExtendedRatioErrorThreshold(0.34);
		}

		int admissibleTextLength() const { return m_matrix.rows()-1; }

		SV* detect_forward_strand(char* text, int text_size, int read_beg, int read_end, int seq_beg, int seq_end) {
			return vectorToPerl(tp_detect_on_strand<equiv_forward_strand>(text, text_size, read_beg, read_end, seq_beg, seq_end, *this));
		}
		SV* detect_reverse_strand(char* text, int text_size, int read_beg, int read_end, int seq_beg, int seq_end) {
			return vectorToPerl(tp_detect_on_strand<equiv_reverse_strand>(text, text_size, read_beg, read_end, seq_beg, seq_end, *this));
		}
		SV* detect_on_strand(char strand, int text_size, char* text, int read_beg, int read_end, int seq_beg, int seq_end) {
			return strand == '+' ? detect_forward_strand(text, text_size, read_beg, read_end, seq_beg, seq_end) :
								   detect_reverse_strand(text, text_size, read_beg, read_end, seq_beg, seq_end);
		}

		int miRnaMinErrorsThreshold() const {
			return m_miRnaMinErrorsThreshold;
		}
		void setMiRnaMinErrorsThreshold(int miRnaMinErrorsThreshold) {
			m_miRnaMinErrorsThreshold = miRnaMinErrorsThreshold;
		}

		int miRnaMaxErrorsThreshold() const {
			return m_miRnaMaxErrorsThreshold;
		}
		void setMiRnaMaxErrorsThreshold(int miRnaMaxErrorsThreshold) {
			m_miRnaMaxErrorsThreshold = miRnaMaxErrorsThreshold;
		}

		int miRnaMinLengthToExtend() const {
			return m_miRnaMinLengthToExtend;
		}
		void setMiRnaMinLengthToExtend(int miRnaMinLengthToExtend) {
			m_miRnaMinLengthToExtend = miRnaMinLengthToExtend;
		}

		int increaseMiRnaEachNt() const {
			return m_increaseMiRnaEachNt;
		}
		void setIncreaseMiRnaEachNt(int increaseMiRnaEachNt) {
			m_increaseMiRnaEachNt = increaseMiRnaEachNt;
		}

		double miRnaExtendedRatioErrorThreshold() const {
			return m_miRnaExtendedRatioErrorThreshold;
		}
		void setMiRnaExtendedRatioErrorThreshold(double miRnaExtendedRatioErrorThreshold) {
			m_miRnaExtendedRatioErrorThreshold = miRnaExtendedRatioErrorThreshold;
		}
};


template <bool (*_equiv_function)(char, char)>
//void computeDistance(char* text, int from_i, int to_i, int from_j, int to_j, MiRnaDetector& self) {
void computeDistance(char* text, Locus read, Locus sequence, TextBounds const& bounds, MiRnaDetector& self) {
	char* read_str = text + read.begin;
	for (int i = read.begin+1; i <= read.end; i++, read_str++) {
		char previousRead = i <= 1 ? 'N' : *(read_str-1);
		char nextRead = i >= bounds ? 'N' : *(read_str+1);
		char* sequence_str = text + sequence.end-1;
		for (std::size_t j = sequence.begin+1; j <= sequence.end; j++, sequence_str--) {
			char previousSeq = j <= 1 ? 'N' : *(sequence_str+1);
			char nextSeq = j >= bounds ? 'N' : *(sequence_str-1);
			process_score<_equiv_function>(self.m_matrix, i, j, *read_str, *sequence_str, 
			_equiv_function(previousRead, previousSeq) || _equiv_function(nextRead, nextSeq));
		}
	}
}

template <bool (*_equiv_function)(char, char)>
//bool checkInwardExtent(char* text, int read_beg, int read_end, int sequence_begin, int seq_end, Locus locus, MiRnaDetector& self) {
ExtensionResult checkInwardExtent(char* text, Locus read, Locus sequence, Locus locus, TextBounds const& bounds, MiRnaDetector& self) {
	Locus seq_locus(sequence.begin + locus.begin, sequence.begin + locus.end);
	Locus newReadLocus(0, 0), newSeqLocus(0, 0);
	int length;
	if (seq_locus.end <= read.begin) { // Read is the 3p
		length = std::ceil(static_cast<float>(read.begin - seq_locus.end)/self.m_increaseMiRnaEachNt);
		length = std::min((uint)length, std::min(read.begin, std::min((uint)self.m_matrix.rows()-1 - seq_locus.end, (seq_locus.end - read.begin)/2)));
		if (length < self.m_miRnaMinLengthToExtend)
			return TooShortToTest;
		newReadLocus.begin = read.begin-length;
		newSeqLocus.begin = seq_locus.end;
	}
	else if (seq_locus.begin >= read.end) { // Read is the 5p
		length = std::ceil(static_cast<float>(seq_locus.begin - read.end)/self.m_increaseMiRnaEachNt);
		length = std::min((uint)length, std::min(seq_locus.begin, std::min((uint)(self.m_matrix.rows()-1 - read.end), (seq_locus.begin - read.end)/2)));
		if (length < self.m_miRnaMinLengthToExtend)
			return TooShortToTest;
		newReadLocus.begin = read.end;
		newSeqLocus.begin = seq_locus.begin - length;
	}
	else
		return ReadAndSeqOverlap;

	newReadLocus.end = newReadLocus.begin + length;
	newSeqLocus.end = newSeqLocus.begin + length;

	self.initRegion(newReadLocus.begin, newReadLocus.end, newSeqLocus.begin, newSeqLocus.end, true);
	computeDistance<_equiv_function>(text, newReadLocus, newSeqLocus, bounds, self);
//	std::cout << "checkInwardExtent: Read: " << newReadLocus << " | Seq: " << newSeqLocus << std::endl;
//	self.m_matrix.print(newReadLocus.begin, newReadLocus.end, newSeqLocus.begin, newSeqLocus.end, text+newReadLocus.begin, text+newSeqLocus.begin);
	return self.lookForResult(newReadLocus, newSeqLocus) ?
				TestPassed : DistanceAboveThreshold;
}

template <bool (*_equiv_function)(char, char)>
ExtensionResult checkOutwardExtent(char* text, Locus read, Locus sequence, Locus locus, TextBounds const& bounds, MiRnaDetector& self) {
	Locus seq_locus(sequence.begin + locus.begin, sequence.begin + locus.end);
	Locus newReadLocus(0, 0), newSeqLocus(0, 0);
	int length;
	if (seq_locus.end <= read.begin) { // Read is the 3p
		length = std::ceil(static_cast<float>(read.begin - seq_locus.end)/self.m_increaseMiRnaEachNt);
		length = std::min((uint)length, std::min(seq_locus.begin, std::min((uint)(self.m_matrix.rows()-1 - read.end), (seq_locus.begin - read.end)/2)));
		if (length < self.m_miRnaMinLengthToExtend)
			return TooShortToTest;
		newReadLocus.begin = read.end;
		newSeqLocus.begin = seq_locus.begin - length;
	}
	else if (seq_locus.begin >= read.end) { // Read is the 5p
		length = std::ceil(static_cast<float>(seq_locus.begin - read.end)/self.m_increaseMiRnaEachNt);
		length = std::min((uint)length, std::min(read.begin, std::min((uint)self.m_matrix.rows()-1 - seq_locus.end, (seq_locus.end - read.begin)/2)));
		if (length < self.m_miRnaMinLengthToExtend)
			return TooShortToTest;
		newReadLocus.begin = read.begin-length;
		newSeqLocus.begin = seq_locus.end;
	}
	else
		return ReadAndSeqOverlap;;

	newReadLocus.end = newReadLocus.begin + length;
	newSeqLocus.end = newSeqLocus.begin + length;

	self.initRegion(newReadLocus.begin, newReadLocus.end, newSeqLocus.begin, newSeqLocus.end, true);
	computeDistance<_equiv_function>(text, newReadLocus, newSeqLocus, bounds, self);
//	std::cout << "checkOutwardExtent: Read: " << newReadLocus << " | Seq: " << newSeqLocus << std::endl;
//	self.m_matrix.print(newReadLocus.begin, newReadLocus.end, newSeqLocus.begin, newSeqLocus.end, text+newReadLocus.begin, text+newSeqLocus.begin);
	return self.lookForResult(newReadLocus, newSeqLocus) ?
				TestPassed : DistanceAboveThreshold;
}
template <bool (*_equiv_function)(char, char)>
std::vector<Locus> tp_detect_on_strand(char* text, int text_size, int read_beg, int read_end, int seq_beg, int seq_end, MiRnaDetector& self) {
	TextBounds bounds(text_size);
	Locus readLocus(read_beg, read_end), seqLocus(seq_beg, seq_end);
	self.initRegion(read_beg, read_end, seq_beg, seq_end, false);
	computeDistance<_equiv_function>(text, readLocus, seqLocus, bounds, self);
//	std::cout << "tp_detect_on_strand" << std::endl;
//	self.m_matrix.print(read_beg, read_end, seq_beg, seq_end, text+read_beg, text+seq_beg);
	bool hasMinThreshold = false;
	std::vector<Locus> candidates = self.lookForMiRnaResult(read_end, seq_beg, seq_end, hasMinThreshold);
//	std::cout << hasMinThreshold << std::endl;
	if (hasMinThreshold)
		return candidates;
	std::vector<Locus> result;
	for (std::vector<Locus>::const_iterator it = candidates.begin(), e = candidates.end(); it != e; ++it) {
//		std::cout << "Candidate: " << *it << std::endl;
		ExtensionResult inwardResult = checkInwardExtent<_equiv_function>(text, readLocus, seqLocus, *it, bounds, self);
		if (inwardResult == TestPassed) {
			result.push_back(*it);
			continue;
		}
		ExtensionResult outwardResult = checkOutwardExtent<_equiv_function>(text, readLocus, seqLocus, *it, bounds, self);
		if (outwardResult == TestPassed ||
				(inwardResult == TooShortToTest && outwardResult == TooShortToTest)) {
			result.push_back(*it);
		}
	}
	return result;
}
// Here ends C++ code
END_OF_C_CODE

sub generate_string {
	my ($length) = @_;
	my $str = '';
	my @alph = qw{A U C G};
	for (my $i = 0; $i < $length; $i++) {
		$str .= $alph[int(rand(4))];
	}
	return $str;
}

sub reverse_complement {
	my $str = shift;
	$str = reverse $str;
	return complement($str);
}

sub complement {
	my $str = shift;
	my %comp = ('A' => 'U', 'U' => 'A', 'T' => 'A', 'C' => 'G', 'G' => 'C');
	for (my $i = 0, my $e = length $str; $i < $e; $i++) {
		substr($str, $i, 1) = $comp{substr($str, $i, 1)};
	}
	return $str;
}

sub generate_string_with_mismatch {
	my ($ref, $mismatch) = @_;
	my $length = length $ref;
	my $str = '';
	my %equiv = ('A' => 'U', 'C' => 'G', 'G' => 'C', 'U' => 'A');
	my %mis = ('A' => 'C', 'C' => 'U', 'G' => 'A', 'U' => 'C');
	for (my $i = 0; $i < $length; $i++) {
		$str .= $equiv{substr($ref, $i, 1)};
	}
	for (my $i = 0, my $e = $mismatch; $i < $e; $i++) {
		substr($str, int(rand($length)), 1) = $mis{substr($ref, $i, 1)};
	}
	return $str;
}

sub generate_matching_candidates {
	my ($lengthA, $lengthB, $w_len, $mismatch) = @_;
	my $read = generate_string($lengthA);
	my $posA = rand($lengthA - $w_len);
	my $posB = rand($lengthB - $w_len);
	my $seq = generate_string($posB);
	$seq .= generate_string_with_mismatch(substr($read, $posA, $w_len), $mismatch);
	$seq .= generate_string($lengthB - $posB - $w_len);
	return ($read, $seq);
}

sub generate_candidates {
	my ($lengthA, $lengthB) = @_;
	return (generate_string($lengthA), generate_string($lengthB));
}

1;
