package MiRnaDuplexDetector;

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

class matrix_int {
	public:
		//typedef unsigned int value_type;

	private:
		std::size_t myRows, myColumns;
		std::vector<unsigned int> myData;

	public:
		matrix_int(std::size_t rows, std::size_t columns) : myRows(rows), myColumns(columns), myData(rows*columns) {}
		matrix_int(std::size_t rows, std::size_t columns, unsigned int value) : myRows(rows), myColumns(columns), myData(rows*columns, value) {}

		std::size_t rows() const { return myRows; }
		std::size_t columns() const { return myColumns; }

		std::size_t size() const { return myData.size(); }

		//EditDistanceScore* data() { return myData.data(); }

		unsigned int& operator()(std::size_t i_row, std::size_t j_column) { return myData[myColumns*i_row+j_column]; }
		unsigned int const& operator()(std::size_t i_row, std::size_t j_column) const { return myData[myColumns*i_row+j_column]; }
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
EditDistanceScore process_score(matrix& score_matrix, std::size_t i, std::size_t j, char read, char sequence) {
	EditDistanceScore s;
	if (_equiv_function(read, sequence))
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

template <bool (*_equiv_function)(char, char)>
unsigned int process_score(matrix_int& score_matrix, std::size_t i, std::size_t j, char read, char sequence) {
	unsigned int s;
	if (_equiv_function(read, sequence))
		s = (score_matrix(i, j) = score_matrix(i-1, j-1));
	else {
		unsigned int subs = score_matrix(i-1,j-1)+1;
		unsigned int insertion = score_matrix(i-1,j)+2;
		unsigned int deletion = score_matrix(i,j-1)+2;
		if (subs <= insertion and subs <= deletion)
			s = (score_matrix(i, j) = subs);
		else if (insertion <= subs and insertion <= deletion)
			s = (score_matrix(i, j) = insertion);
		else
			s = (score_matrix(i, j) = deletion);
	}
	return s;
}

EditDistanceScore process_score_forward_strand(matrix& score_matrix, std::size_t i, std::size_t j, char read, char sequence) {
	return process_score<equiv_forward_strand>(score_matrix, i, j, read, sequence);
}
EditDistanceScore process_score_reverse_strand(matrix& score_matrix, std::size_t i, std::size_t j, char read, char sequence) {
	return process_score<equiv_reverse_strand>(score_matrix, i, j, read, sequence);
}

unsigned int process_score_forward_strand(matrix_int& score_matrix, std::size_t i, std::size_t j, char read, char sequence) {
	return process_score<equiv_forward_strand>(score_matrix, i, j, read, sequence);
}
unsigned int process_score_reverse_strand(matrix_int& score_matrix, std::size_t i, std::size_t j, char read, char sequence) {
	return process_score<equiv_reverse_strand>(score_matrix, i, j, read, sequence);
}

// EDIT DISTANCE EXTENDED

#include <limits>

struct Locus {
		unsigned int begin, end;
		Locus(unsigned int b, unsigned int e) : begin(b), end(e) {}

		SV* toPerl() const { return __get_pos(begin, end); }
};

#include <deque>

SV* dequeToPerl(std::deque<Locus> const& loci) {
	AV* results = newAV();
	for (std::deque<Locus>::const_iterator it = loci.begin(), e = loci.end(); it != e; ++it) {
		av_push(results, it->toPerl());
	}
	return newRV_noinc((SV*)results);
}

class MiRnaDetector;

template <bool (*_equiv_function)(char, char)>
void computeDistance(char* text, int from_i, int to_i, int from_j, int to_j, MiRnaDetector& self);
template <bool (*_equiv_function)(char, char)>
bool checkInwardExtent(char* text, int read_beg, int read_end, int seq_beg, int seq_end, Locus locus, MiRnaDetector& self);
template <bool (*_equiv_function)(char, char)>
bool checkOutwardExtent(char* text, int read_beg, int read_end, int seq_beg, int seq_end, Locus locus, MiRnaDetector& self);
template <bool (*_equiv_function)(char, char)>
std::deque<Locus> tp_detect_on_strand(char* text, int read_beg, int read_end, int seq_beg, int seq_end, MiRnaDetector& self);

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

		std::deque<Locus> lookForMiRnaResult(int to_i, int from_j, int to_j, bool& meetMinThreshold) {
			meetMinThreshold = false;
			EditDistanceScore* data = m_matrix.scanRow(to_i) + from_j+1;
			std::deque<Locus> result;
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
		bool lookForResult(int to_i, int from_j, int to_j, int threshold) const {
			EditDistanceScore const* data = m_matrix.scanRow(to_i) + from_j+1;
			for (int j = 1, e = to_j-from_j; j <= e; j++, data++) {
				if (data->value <= threshold)
					return true;
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

		SV* detect_forward_strand(char* text, int read_beg, int read_end, int seq_beg, int seq_end) {
			return dequeToPerl(tp_detect_on_strand<equiv_forward_strand>(text, read_beg, read_end, seq_beg, seq_end, *this));
		}
		SV* detect_reverse_strand(char* text, int read_beg, int read_end, int seq_beg, int seq_end) {
			return dequeToPerl(tp_detect_on_strand<equiv_reverse_strand>(text, read_beg, read_end, seq_beg, seq_end, *this));
		}
		SV* detect_on_strand(char strand, char* text, int read_beg, int read_end, int seq_beg, int seq_end) {
			return strand == '+' ? detect_forward_strand(text, read_beg, read_end, seq_beg, seq_end) :
								   detect_reverse_strand(text, read_beg, read_end, seq_beg, seq_end);
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
void computeDistance(char* text, int from_i, int to_i, int from_j, int to_j, MiRnaDetector& self) {
	char* read = text + from_i;
	for (int i = from_i+1; i <= to_i; i++, read++) {
		char* sequence = text + to_j-1;
		for (std::size_t j = from_j+1; j <= to_j; j++, sequence--)
			process_score<_equiv_function>(self.m_matrix, i, j, *read, *sequence);
	}
}

template <bool (*_equiv_function)(char, char)>
bool checkInwardExtent(char* text, int read_beg, int read_end, int seq_beg, int seq_end, Locus locus, MiRnaDetector& self) {
	Locus seq_locus(seq_beg + locus.begin, seq_beg + locus.end);
	int length, new_seq_beg, new_read_beg;
	if (seq_locus.end <= read_beg) { // Read is the 3p
		length = std::ceil(static_cast<float>(read_beg - seq_locus.end)/self.m_increaseMiRnaEachNt);
		if (read_beg < length)
			length = read_beg;
		if (seq_locus.end+length > self.m_matrix.rows()-1)
			length = self.m_matrix.rows()-1 - seq_locus.end;
		if (length < self.m_miRnaMinLengthToExtend)
			return true;
		new_read_beg = read_beg-length;
		new_seq_beg = seq_locus.end;
	}
	else if (seq_locus.begin >= read_end) { // Read is the 5p
		length = std::ceil(static_cast<float>(seq_locus.begin - read_end)/self.m_increaseMiRnaEachNt);
		if (seq_locus.begin < length)
			length = seq_locus.begin;
		if (read_end+length > self.m_matrix.rows()-1)
			length = self.m_matrix.rows()-1 - read_end;
		if (length < self.m_miRnaMinLengthToExtend)
			return true;
		new_read_beg = read_end;
		new_seq_beg = seq_locus.begin - length;
	}
	else
		return false;

	self.initRegion(new_read_beg, new_read_beg+length, new_seq_beg, new_seq_beg+length, true);
	computeDistance<_equiv_function>(text, new_read_beg, new_read_beg+length, new_seq_beg, new_seq_beg+length, self);
	return self.lookForResult(new_read_beg+length, new_seq_beg, new_seq_beg+length, length*self.m_miRnaExtendedRatioErrorThreshold);
}

template <bool (*_equiv_function)(char, char)>
bool checkOutwardExtent(char* text, int read_beg, int read_end, int seq_beg, int seq_end, Locus locus, MiRnaDetector& self) {
	Locus seq_locus(seq_beg + locus.begin, seq_beg + locus.end);
	int length, new_seq_beg, new_read_beg;
	if (seq_locus.end <= read_beg) { // Read is the 3p
		length = std::ceil(static_cast<float>(read_beg - seq_locus.end)/self.m_increaseMiRnaEachNt);
		if (seq_locus.begin < length)
			length = seq_locus.begin;
		if (read_end+length > self.m_matrix.rows()-1)
			length = self.m_matrix.rows()-1 - read_end;
		if (length < self.m_miRnaMinLengthToExtend)
			return true;
		new_read_beg = read_end;
		new_seq_beg = seq_locus.begin - length;
	}
	else if (seq_locus.begin >= read_end) { // Read is the 5p
		length = std::ceil(static_cast<float>(seq_locus.begin - read_end)/self.m_increaseMiRnaEachNt);
		if (read_beg < length)
			length = read_beg;
		if (seq_locus.end+length > self.m_matrix.rows()-1)
			length = self.m_matrix.rows()-1 - seq_locus.end;
		if (length < self.m_miRnaMinLengthToExtend)
			return true;
		new_read_beg = read_beg-length;
		new_seq_beg = seq_locus.end;
	}
	else
		return false;

	self.initRegion(new_read_beg, new_read_beg+length, new_seq_beg, new_seq_beg+length, true);
	computeDistance<_equiv_function>(text, new_read_beg, new_read_beg+length, new_seq_beg, new_seq_beg+length, self);
	return self.lookForResult(new_read_beg+length, new_seq_beg, new_seq_beg+length, length*self.m_miRnaExtendedRatioErrorThreshold);
}

template <bool (*_equiv_function)(char, char)>
std::deque<Locus> tp_detect_on_strand(char* text, int read_beg, int read_end, int seq_beg, int seq_end, MiRnaDetector& self) {
	self.initRegion(read_beg, read_end, seq_beg, seq_end, false);
	computeDistance<_equiv_function>(text, read_beg, read_end, seq_beg, seq_end, self);
	bool hasMinThreshold = false;
	std::deque<Locus> candidates = self.lookForMiRnaResult(read_end, seq_beg, seq_end, hasMinThreshold);
	if (hasMinThreshold)
		return candidates;
	std::deque<Locus> result;
	for (std::deque<Locus>::const_iterator it = candidates.begin(), e = candidates.end(); it != e; ++it) {
		if (checkInwardExtent<_equiv_function>(text, read_beg, read_end, seq_beg, seq_end, *it, self)) {
			result.push_back(*it);
		}
		if (checkOutwardExtent<_equiv_function>(text, read_beg, read_end, seq_beg, seq_end, *it, self)) {
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
