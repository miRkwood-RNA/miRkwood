package MiRnaDuplexDetector;

use strict;
use warnings;
use Inline (Config => DIRECTORY => '/tmp/',);
use Inline CPP => <<'END_OF_C_CODE';
// Here starts C++ code

#include <vector>
#include <utility>
#include <cmath>
#include <iostream>

int ntToIndex(char nt) {
	switch (nt) {
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
		case 'U':
			return 3;
		default:
			return 4;
	}
}

struct forward_strand {};
struct reverse_strand {};

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

template <class>
struct cost_function;

#define INDEL -3
#define SUBS -2

template<>
struct cost_function<forward_strand> {

		static int match[5][5];

		static int matchCost(char a, char b) { return match[ntToIndex(a)][ntToIndex(b)]; }

		static int const insertion = INDEL;
		static int const deletion = INDEL;

		static bool equiv(char a, char b) {
			return equiv_forward_strand(a, b);
		}

};
template<>
struct cost_function<reverse_strand> {

		static int match[5][5];

		static int matchCost(char a, char b) { return match[ntToIndex(a)][ntToIndex(b)]; }

		static int const insertion = INDEL;
		static int const deletion = INDEL;

		static bool equiv(char a, char b) {
			return equiv_reverse_strand(a, b);
		}

};

int cost_function<forward_strand>::match[5][5] = {
//		A		C		G		U		N
/*A*/	{SUBS,	SUBS,	SUBS,	2,		SUBS},
/*C*/	{SUBS,	SUBS,	3,		SUBS,	SUBS},
/*G*/	{SUBS,	3,		SUBS,	1,		SUBS},
/*U*/	{2,		SUBS,	1,		SUBS,	SUBS},
/*N*/	{SUBS,	SUBS,	SUBS,	SUBS,	SUBS}
};

int cost_function<reverse_strand>::match[5][5] = {
//		A		C		G		U		N
/*A*/	{SUBS,	1,		SUBS,	2,		SUBS},
/*C*/	{1,		SUBS,	3,		SUBS,	SUBS},
/*G*/	{SUBS,	3,		SUBS,	SUBS,	SUBS},
/*U*/	{2,		SUBS,	SUBS,	SUBS,	SUBS},
/*N*/	{SUBS,	SUBS,	SUBS,	SUBS,	SUBS}
};

#undef SUBS
#undef INDEL

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


struct Locus {
		unsigned int begin, end;
		Locus(unsigned int b, unsigned int e) : begin(b), end(e) {}

		SV* toPerl() const { return __get_pos(begin, end); }
};

enum ExtensionResult {
	DistanceAboveThreshold, TestPassed, TooShortToTest, ReadAndSeqOverlap
};

inline std::ostream& operator<<(std::ostream& os, Locus const& pos) {
	return os << "Locus: (" << pos.begin << " -> " << pos.end << ")";
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

		void print(std::size_t from_i, std::size_t to_i, std::size_t from_j, std::size_t to_j, char* read, char* seq) const {
			std::cout << "x\te\t";
			seq += to_j-from_j-1;
			for (std::size_t j = from_j+1; j <= to_j; ++j, --seq) {
				std::cout << *seq << "\t";
			}
			std::cout << std::endl;
			for (std::size_t i = from_i; i <= to_i; ++i) {
				if (i == from_i)
					std::cout << "e\t";
				else
					std::cout << *read++ << "\t";
				for (std::size_t j = from_j; j <= to_j; ++j) {
					std::cout << (*this)(i, j).value << "\t";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
};

template <class _Strand>
EditDistanceScore process_score(matrix& score_matrix, std::size_t i, std::size_t j, char read, char sequence, bool allowBinding) {
	EditDistanceScore s;
//	if (_equiv_function(read, sequence) && allowBinding)
//		s = (score_matrix(i, j) = score_matrix(i-1, j-1));
//	else {
		EditDistanceScore subs = score_matrix(i-1,j-1).increasedBy(cost_function<_Strand>::matchCost(allowBinding ? read : 'N', sequence));
		EditDistanceScore insertion = score_matrix(i-1,j).increasedBy(cost_function<_Strand>::insertion);
		EditDistanceScore deletion = score_matrix(i,j-1).increasedBy(cost_function<_Strand>::deletion);
		if (subs.value >= insertion.value and subs.value >= deletion.value)
			s = (score_matrix(i, j) = subs);
		else if (insertion.value >= subs.value and insertion.value >= deletion.value)
			s = (score_matrix(i, j) = insertion);
		else
			s = (score_matrix(i, j) = deletion);
//	}
	return s;
}

// EDIT DISTANCE EXTENDED

#include <limits>
#include <algorithm>
#include <vector>


SV* vectorToPerl(std::vector<Locus> const& loci) {
	AV* results = newAV();
	for (std::vector<Locus>::const_iterator it = loci.begin(), e = loci.end(); it != e; ++it) {
		av_push(results, it->toPerl());
	}
	return newRV_noinc((SV*)results);
}

//struct TextBounds {
//		Locus read, sequence;

//		TextBounds(Locus r, Locus s) : read(r), sequence(s) {}
//};

typedef int TextBounds;

class MiRnaDetector;

template <class _Strand>
void computeDistance(char* text, Locus read, Locus sequence, TextBounds const& bounds, MiRnaDetector& self);
template <class _Strand>
ExtensionResult checkInwardExtent(char* text, Locus read, Locus sequence, Locus locus, TextBounds const& bounds, MiRnaDetector& self);
template <class _Strand>
ExtensionResult checkOutwardExtent(char* text, Locus read, Locus sequence, Locus locus, TextBounds const& bounds, MiRnaDetector& self);
template <class _Strand>
std::vector<Locus> tp_detect_on_strand(char* text, int text_size, Locus read, Locus sequence, MiRnaDetector& self);

class MiRnaDetector {

	public:

		matrix m_matrix;
		int m_miRnaHigherScoreThreshold, m_miRnaLowerScoreThreshold;
		int m_miRnaMinLengthToExtend;
		int m_increaseMiRnaEachNt;
		double m_miRnaExtendedRatioScoreThreshold;

		void initRegion(Locus read, Locus sequence, bool realDistance) {
			for (int i = 0, e = read.end-read.begin; i <= e; i++) {
				EditDistanceScore* data = m_matrix.scanRow(i+read.begin) + sequence.begin;
				*data = EditDistanceScore(i*cost_function<forward_strand>::insertion, 0);
			}
			if (realDistance) {
				EditDistanceScore* data = m_matrix.scanRow(read.begin) + sequence.begin+1;
				for (int j = 1, e = sequence.end-sequence.begin; j <= e; j++, data++)
					*data = EditDistanceScore(j*cost_function<forward_strand>::deletion, j);
			}
			else {
				EditDistanceScore* data = m_matrix.scanRow(read.begin) + sequence.begin+1;
				for (int j = 1, e = sequence.end-sequence.begin; j <= e; j++, data++)
					*data = EditDistanceScore(0, j);
			}
		}

		std::vector<Locus> lookForMiRnaResult(int to_i, Locus sequence, bool& meetHigherThreshold) {
			meetHigherThreshold = false;
			EditDistanceScore* data = m_matrix.scanRow(to_i) + sequence.begin+1;
			std::vector<Locus> result;
			for (int j = 1, e = sequence.end-sequence.begin; j <= e; j++, data++) {
				if (data->value >= m_miRnaHigherScoreThreshold) {
					if (!meetHigherThreshold)
						result.clear();
					meetHigherThreshold = true;
					result.push_back(Locus(e-j, e-data->from_j));
				}
				else if (!meetHigherThreshold && data->value >= m_miRnaLowerScoreThreshold)
					result.push_back(Locus(e-j, e-data->from_j));
			}
			return result;
		}
		bool lookForResult(Locus read, Locus sequence) const {
			// Last row
			EditDistanceScore const* data = m_matrix.scanRow(read.end) + sequence.begin+1;
			for (uint j = sequence.begin+1; j <= sequence.end; j++, data++) {
				if (data->value >= m_miRnaExtendedRatioScoreThreshold*std::max(read.end-read.begin, j-sequence.begin))
					return true;
			}
			// Last column
			for (uint i = read.begin+1; i <= read.end; i++) {
				if (m_matrix(i, sequence.end).value >= m_miRnaExtendedRatioScoreThreshold*std::max(i-read.begin, sequence.end-sequence.begin))
					return true;
			}
			return false;
		}

	public:
		MiRnaDetector(int textSize) : m_matrix(textSize+1, textSize+1, EditDistanceScore(0, 0)) {
			setMiRnaHigherScoreThreshold(32);
			setMiRnaLowerScoreThreshold(14);
			setMiRnaMinLengthToExtend(24);
			setIncreaseMiRnaEachNt(1);
			setMiRnaExtendedRatioScoreThreshold(0.85);
		}

		int admissibleTextLength() const { return m_matrix.rows()-1; }

		SV* detect_forward_strand(char* text, int text_size, int read_beg, int read_end, int seq_beg, int seq_end) {
			return vectorToPerl(tp_detect_on_strand<forward_strand>(text, text_size, Locus(read_beg, read_end), Locus(seq_beg, seq_end), *this));
		}
		SV* detect_reverse_strand(char* text, int text_size, int read_beg, int read_end, int seq_beg, int seq_end) {
			return vectorToPerl(tp_detect_on_strand<reverse_strand>(text, text_size, Locus(read_beg, read_end), Locus(seq_beg, seq_end), *this));
		}
		SV* detect_on_strand(char strand, char* text, int text_size, int read_beg, int read_end, int seq_beg, int seq_end) {
			return strand == '+' ? detect_forward_strand(text, text_size, read_beg, read_end, seq_beg, seq_end) :
								   detect_reverse_strand(text, text_size, read_beg, read_end, seq_beg, seq_end);
		}

		int miRnaHigherScoreThreshold() const {
			return m_miRnaHigherScoreThreshold;
		}
		void setMiRnaHigherScoreThreshold(int miRnaHigherScoreThreshold) {
			m_miRnaHigherScoreThreshold = miRnaHigherScoreThreshold;
		}

		int miRnaLowerScoreThreshold() const {
			return m_miRnaLowerScoreThreshold;
		}
		void setMiRnaLowerScoreThreshold(int miRnaLowerScoreThreshold) {
			m_miRnaLowerScoreThreshold = miRnaLowerScoreThreshold;
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

		double miRnaExtendedRatioScoreThreshold() const {
			return m_miRnaExtendedRatioScoreThreshold;
		}
		void setMiRnaExtendedRatioScoreThreshold(double miRnaExtendedRatioScoreThreshold) {
			m_miRnaExtendedRatioScoreThreshold = miRnaExtendedRatioScoreThreshold;
		}
};


template <class _Strand>
void computeDistance(char* text, Locus read, Locus sequence, TextBounds const& bounds, MiRnaDetector& self) {
//	std::cout << "Compute distance between read: " << read << " and sequence: " << sequence << std::endl;
	char* read_str = text + read.begin;
	for (int i = read.begin+1; i <= read.end; i++, read_str++) {
		char previousRead = i <= 1 ? 'N' : *(read_str-1);
		char nextRead = i >= bounds ? 'N' : *(read_str+1);
		char* sequence_str = text + sequence.end-1;
		for (std::size_t j = sequence.begin+1; j <= sequence.end; j++, sequence_str--) {
			char previousSeq = j <= 1 ? 'N' : *(sequence_str+1);
			char nextSeq = j >= bounds ? 'N' : *(sequence_str-1);
			process_score<_Strand>(self.m_matrix, i, j, *read_str, *sequence_str,
			 cost_function<_Strand>::equiv(previousRead, previousSeq) || cost_function<_Strand>::equiv(nextRead, nextSeq));
		}
	}
}

template <class _Strand>
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

	self.initRegion(newReadLocus, newSeqLocus, true);
	computeDistance<_Strand>(text, newReadLocus, newSeqLocus, bounds, self);
	return self.lookForResult(newReadLocus, newSeqLocus) ?
				TestPassed : DistanceAboveThreshold;
}

template <class _Strand>
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

	self.initRegion(newReadLocus, newSeqLocus, true);
	computeDistance<_Strand>(text, newReadLocus, newSeqLocus, bounds, self);
	return self.lookForResult(newReadLocus, newSeqLocus) ?
				TestPassed : DistanceAboveThreshold;
}

template <class _Strand>
std::vector<Locus> tp_detect_on_strand(char* text, int text_size, Locus readLocus, Locus seqLocus, MiRnaDetector& self) {
	TextBounds bounds(text_size);
	self.initRegion(readLocus, seqLocus, false);
	computeDistance<_Strand>(text, readLocus, seqLocus, bounds, self);
	bool hasMinThreshold = false;
	std::vector<Locus> candidates = self.lookForMiRnaResult(readLocus.end, seqLocus, hasMinThreshold);
	if (hasMinThreshold)
		return candidates;
	std::vector<Locus> result;
	for (std::vector<Locus>::const_iterator it = candidates.begin(), e = candidates.end(); it != e; ++it) {
		ExtensionResult inwardResult = checkInwardExtent<_Strand>(text, readLocus, seqLocus, *it, bounds, self);
		if (inwardResult == TestPassed) {
			result.push_back(*it);
			continue;
		}
		ExtensionResult outwardResult = checkOutwardExtent<_Strand>(text, readLocus, seqLocus, *it, bounds, self);
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
	my @alph = ('A', 'U', 'C', 'G');
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
