package miRkwood::CppBinarySearch;

use strict;
use warnings;
use Inline (Config => DIRECTORY => '/tmp/',);
use Inline CPP => <<'END_OF_C_CODE';
// Here starts C++ code

#include <algorithm>
#include <iterator>

#define PERL_HASH_KEY_COMP 			"begin"
#define PERL_HASH_KEY_COMP_LENGTH 5

class perl_array_of_hash_ref_iterator : public std::iterator<std::random_access_iterator_tag, IV> {
	
	public:
	
	AV* array;
	SSize_t index;
	
		perl_array_of_hash_ref_iterator(AV* arr, SSize_t s = 0) : array(arr), index(s) {}
		perl_array_of_hash_ref_iterator(perl_array_of_hash_ref_iterator const& other) : array(other.array), index(other.index) {}
		
		IV operator*() const { 
			HV* hash = (HV*)SvRV(*av_fetch(array, index, 0));
			return SvIV(*hv_fetch(hash, PERL_HASH_KEY_COMP, PERL_HASH_KEY_COMP_LENGTH, 0));
		}
		
		perl_array_of_hash_ref_iterator& operator++() {
			++index;
			return *this;
		}
		perl_array_of_hash_ref_iterator operator++(int) {
			perl_array_of_hash_ref_iterator it(array, index++);
			return it;
		}
		
		perl_array_of_hash_ref_iterator& operator--() {
			--index;
			return *this;
		}
		perl_array_of_hash_ref_iterator operator--(int) {
			perl_array_of_hash_ref_iterator it(array, index--);
			return it;
		}
		
		bool operator==(const perl_array_of_hash_ref_iterator &o) const
			{ return index == o.index; }
		bool operator!=(const perl_array_of_hash_ref_iterator &o) const
			{ return index != o.index; }
		bool operator<(const perl_array_of_hash_ref_iterator& other) const
			{ return index < other.index; }
		bool operator<=(const perl_array_of_hash_ref_iterator& other) const
			{ return index <= other.index; }
		bool operator>(const perl_array_of_hash_ref_iterator& other) const
			{ return index > other.index; }
		bool operator>=(const perl_array_of_hash_ref_iterator& other) const
			{ return index >= other.index; }
			
		perl_array_of_hash_ref_iterator &operator+=(int j) { index+=j; return *this; }
		perl_array_of_hash_ref_iterator &operator-=(int j) { index-=j; return *this; }
		perl_array_of_hash_ref_iterator operator+(int j) const { return perl_array_of_hash_ref_iterator(array, index+j); }
		perl_array_of_hash_ref_iterator operator-(int j) const { return perl_array_of_hash_ref_iterator(array, index-j); }
		int operator-(perl_array_of_hash_ref_iterator j) const { return index - j.index; }
	
};

// parsed_bed: an array reference of hash reference. Each hash must have a 'begin' key.
// The array should be sorted in ascending order with regard to the 'begin' key. (i.e. {$a->{'begin'} <=> $b->{'begin'}})
// Returns the index of the first element in the range [begin, end) that is not less than (i.e. greater or equal to) value.
// Returns end if no such element is found
int lower_bound(SV* parsed_bed, int begin, int end, IV value) {
	AV* array = (AV*)SvRV(parsed_bed);
	perl_array_of_hash_ref_iterator first(array, begin);
	perl_array_of_hash_ref_iterator last(array, end);
	return std::lower_bound(first, last, value).index;
}

// parsed_bed: an array reference of hash reference. Each hash must have a 'begin' key.
// The array should be sorted in ascending order with regard to the 'begin' key. (i.e. {$a->{'begin'} <=> $b->{'begin'}})
// Returns the index of the first element in the range [begin, end) that is not less than (i.e. greater or equal to) value.
// Returns end if no such element is found
int upper_bound(SV* parsed_bed, int begin, int end, IV value) {
	AV* array = (AV*)SvRV(parsed_bed);
	perl_array_of_hash_ref_iterator first(array, begin);
	perl_array_of_hash_ref_iterator last(array, end);
	return std::upper_bound(first, last, value).index;
}

// Here ends C++ code
END_OF_C_CODE

1;
