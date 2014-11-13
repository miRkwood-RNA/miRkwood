package MiRnaDuplexDetector;

use strict;
use warnings;

use Inline CPP => <<'END_OF_C_CODE';
// Here starts C++ code

#include <cstring>
#include <deque>

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

// This function tries to map the reversed 'read' dna sequence to the 'sequence' dna sequence.
// A match is found if the sequences can bind on a window of length 'w_len', allowing up to 'allowedMismatches' mismatches. (Insertions and deletions are not considered, only substitutions are).
// 'read' is automatically reversed by this function, so you must not do it manually.
// Returns an array reference of array references:
//		Returns [] if there are no results
//		If there is a single match found, returns [[$read_id, $sequence_id]]; This means that '$reversed_read = reverse substr($read, $read_id, $w_len)' can bind to 'substr($sequence, $sequence_id, $w_len)'
//		If N matches are found, returns [[$read_id_1, $sequence_id_1], ..., [$read_id_N, $sequence_id_N]]
//		
// Best value for w_len = 12
// Best value for allowedMismatches = 2
// You may call this function from Perl. Simply write '$result = find_match_forward_strand($read, $sequence, $w_len, $allowedMismatches);' and it will work.
SV* find_match_forward_strand(char* read, char* sequence, int w_len, int allowedMismatches) {
	int read_length, seq_length, i, e_seq, seq_id, read_id, reversed_read_id, counter, counter_end;
	AV* results;
	read_length = std::strlen(read);
	seq_length = std::strlen(sequence);
	results = newAV();
	for (i = -read_length+w_len, e_seq = seq_length - w_len + 1; i < e_seq; i++) {
		std::deque<int> substitution;
		seq_id = std::max(i, 0);
		read_id = std::max(-i, 0);
		reversed_read_id = read_length-1-read_id;
		int returned_seq_id = seq_id+1-w_len;
		counter = 0;
		counter_end = std::min(seq_length-seq_id, read_length-read_id);
		for (; counter < w_len; counter++) {
			if (!equiv_forward_strand(read[reversed_read_id-counter], sequence[seq_id+counter]))
				substitution.push_back(counter);
		}
		if (substitution.size() <= allowedMismatches)
			av_push(results, __get_pos(reversed_read_id-w_len+1, seq_id));
		for (; counter < counter_end - w_len; counter++) {
			if (substitution.size() && substitution.front() == counter - w_len)
				substitution.pop_front();
			if (!equiv_forward_strand(read[reversed_read_id-counter], sequence[seq_id+counter]))
				substitution.push_back(counter);
			if (substitution.size() <= allowedMismatches)
				av_push(results, __get_pos(reversed_read_id-counter, returned_seq_id+counter));
		}
		int errors_to_remove = substitution.size();
		for (; counter < counter_end; counter++) {
			if (substitution.size() && substitution.front() == counter - w_len) {
				substitution.pop_front();
				errors_to_remove--;
			}
			if (!equiv_forward_strand(read[reversed_read_id-counter], sequence[seq_id+counter]))
				substitution.push_back(counter);
			if (substitution.size() <= allowedMismatches)
				av_push(results, __get_pos(reversed_read_id-counter, returned_seq_id+counter));
			else if (substitution.size() - errors_to_remove > allowedMismatches)
				break;
		}
	}
	return newRV_noinc((SV*)results);
}

// This function tries to map the complement of 'read' dna sequence to the reverse complement of 'sequence' dna sequence. Both sequences are complemented by the function so you mustn't complement them.
// 'sequence' is also reversed by the function so mustn't reverse it beforehand as well.
// A match is found if the sequences can bind on a window of length 'w_len', allowing up to 'allowedMismatches' mismatches. (Insertions and deletions are not considered, only substitutions are).
// Returns an array reference of array references:
//		Returns [] if there are no results
//		If there is a single match found, returns [[$read_id, $sequence_id]]; This means that 'complement substr($read, $read_id, $w_len)' can
//			bind to 'reverse_complement substr($sequence, $sequence_id, $w_len)'
//		If N matches are found, returns [[$read_id_1, $sequence_id_1], ..., [$read_id_N, $sequence_id_N]]
//
// Best value for w_len = 12
// Best value for allowedMismatches = 2
// You may call this function from Perl. Simply write '$result = find_match_reverse_strand($read, $sequence, $w_len, $allowedMismatches);' and it will work.
SV* find_match_reverse_strand(char* read, char* sequence, int w_len, int allowedMismatches) {
	int read_length, seq_length, i, e_seq, seq_id, read_id, reversed_seq_id, counter, counter_end;
	AV* results;
	read_length = std::strlen(read);
	seq_length = std::strlen(sequence);
	results = newAV();
	for (i = -read_length+w_len, e_seq = seq_length - w_len + 1; i < e_seq; i++) {
		std::deque<int> substitution;
		seq_id = std::max(i, 0);
		read_id = std::max(-i, 0);
		reversed_seq_id = seq_length-1-seq_id;
		int returned_read_id = read_id+1-w_len;
		counter = 0;
		counter_end = std::min(seq_length-seq_id, read_length-read_id);
		for (; counter < w_len; counter++) {
			if (!equiv_reverse_strand(read[read_id+counter], sequence[reversed_seq_id-counter]))
				substitution.push_back(counter);
		}
		if (substitution.size() <= allowedMismatches)
			av_push(results, __get_pos(read_id, reversed_seq_id-w_len+1));
		for (; counter < counter_end - w_len; counter++) {
			if (substitution.size() && substitution.front() == counter - w_len)
				substitution.pop_front();
			if (!equiv_reverse_strand(read[read_id+counter], sequence[reversed_seq_id-counter]))
				substitution.push_back(counter);
			if (substitution.size() <= allowedMismatches)
				av_push(results, __get_pos(returned_read_id+counter, reversed_seq_id-counter));
		}
		int errors_to_remove = substitution.size();
		for (; counter < counter_end; counter++) {
			if (substitution.size() && substitution.front() == counter - w_len) {
				substitution.pop_front();
				errors_to_remove--;
			}
			if (!equiv_reverse_strand(read[read_id+counter], sequence[reversed_seq_id-counter]))
				substitution.push_back(counter);
			if (substitution.size() <= allowedMismatches)
				av_push(results, __get_pos(returned_read_id+counter, reversed_seq_id-counter));
			else if (substitution.size() - errors_to_remove > allowedMismatches)
				break;
		}
	}
	return newRV_noinc((SV*)results);
}

// This function calls find_match_forward_strand or find_match_reverse_strand depending on the strand parameter
// If strand is '+', calls find_match_forward_strand
// Otherwise calls find_match_reverse_strand
SV* find_match(char strand, char* read, char* sequence, int w_len, int allowedMismatches) {
	return strand == '+' ? find_match_forward_strand(read, sequence, w_len, allowedMismatches) : find_match_reverse_strand(read, sequence, w_len, allowedMismatches);
}

// EDIT DISTANCE
#include <vector>
#include <utility>

struct EditDistanceScore {
		unsigned int value, from_i, from_j;

		EditDistanceScore(unsigned int v = 0, unsigned int i = 0, unsigned int j = 0) : value(v), from_i(i), from_j(j) {}

		EditDistanceScore& increased() { value++; return *this; }
		EditDistanceScore& increasedBy(int amount) { value+=amount; return *this; }
};

class matrix {
	public:
		//typedef EditDistanceScore T;

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

		EditDistanceScore& operator()(std::size_t i_row, std::size_t j_column) { return myData[myRows*j_column+i_row]; }
		EditDistanceScore const& operator()(std::size_t i_row, std::size_t j_column) const { return myData[myRows*j_column+i_row]; }
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
		EditDistanceScore subs = score_matrix(i-1,j-1);
		EditDistanceScore insertion = score_matrix(i-1,j);
		EditDistanceScore deletion = score_matrix(i,j-1);
		if (subs.value+1 <= insertion.value+2 and subs.value+1 <= deletion.value+2)
			s = (score_matrix(i, j) = subs.increased());
		else if (insertion.value+2 <= subs.value+1 and insertion.value+2 <= deletion.value+2)
			s = (score_matrix(i, j) = insertion.increasedBy(2));
		else
			s = (score_matrix(i, j) = deletion.increasedBy(2));
	}
	return s;
}

EditDistanceScore process_score_forward_strand(matrix& score_matrix, std::size_t i, std::size_t j, char read, char sequence) {
	return process_score<equiv_forward_strand>(score_matrix, i, j, read, sequence);
}

EditDistanceScore process_score_reverse_strand(matrix& score_matrix, std::size_t i, std::size_t j, char read, char sequence) {
	return process_score<equiv_reverse_strand>(score_matrix, i, j, read, sequence);
}

typedef matrix dummy;

// The extra dummy* parameter is to prevent Inline CPP from binding this function to Perl
// This is because Inline CPP ignores the template parameter and hence fails to bind to Perl, resulting in compiling errors.
template <EditDistanceScore (*process_score_func)(matrix&, std::size_t, std::size_t, char, char)>
SV* __edit_distance_on_strand(char* read, char* sequence, unsigned int max_distance, dummy*) {
	unsigned int read_length = std::strlen(read);
	unsigned int sequence_length = std::strlen(sequence);
	matrix score_matrix(read_length+1, sequence_length+1, EditDistanceScore(0, 0, 0));

	for (std::size_t i = 1; i <= read_length; i++)
		score_matrix(i, 0) = EditDistanceScore(i, read_length-i, 0);
	for (std::size_t j = 1; j <= sequence_length; j++)
		score_matrix(0, j) = EditDistanceScore(0, read_length-1, j-1);

	unsigned int min_score = ~0u;
	AV* result = newAV();
	std::size_t i = 1;
	for (; i < read_length; i++) for (std::size_t j = 1; j <= sequence_length; j++)
		process_score_func(score_matrix, i, j, read[read_length-i], sequence[j-1]);
	EditDistanceScore current_score;
	for (std::size_t j = 1; j <= sequence_length; j++) {
		current_score = process_score_func(score_matrix, i, j, read[read_length-i], sequence[j-1]);
		if (current_score.value < min_score and current_score.value <= max_distance) {
			min_score = current_score.value;
			av_clear(result);
			av_push(result, __get_region(read_length-i, current_score.from_i+1, current_score.from_j, j));
		}
		else if (current_score.value == min_score)
			av_push(result, __get_region(read_length-i, current_score.from_i+1, current_score.from_j, j));
	}
	return newRV_noinc((SV*)result);
}

SV* edit_distance_forward_strand(char* read, char* sequence, unsigned int max_distance) {
	return __edit_distance_on_strand<process_score_forward_strand>(read, sequence, max_distance, 0);
}

SV* edit_distance_reverse_strand(char* read, char* sequence, unsigned int max_distance) {
	return __edit_distance_on_strand<process_score_reverse_strand>(read, sequence, max_distance, 0);
}

SV* edit_distance_on_strand(char strand, char* read, char* sequence, unsigned int max_distance) {
	return strand == '+' ? edit_distance_forward_strand(read, sequence, max_distance) : edit_distance_reverse_strand(read, sequence, max_distance);
}

// Here ends C++ code
END_OF_C_CODE

sub generate_string {
	my ($length) = @_;
	my $str = '';
	my @alph = ('A', 'U', 'C', 'G');
	for (my $i = 0; $i < $length; $i++) {
		$str .= $alph[rand(4)];
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
		substr($str, rand($length), 1) = $mis{substr($ref, $i, 1)};
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