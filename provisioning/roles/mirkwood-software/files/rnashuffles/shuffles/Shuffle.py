#! /usr/bin/env python
# -*- coding: latin-1 -*-

"""Implementation of the Altschul-Erikson dinucleotide shuffle algorithm.

Builds shuffled sequences by permuting a single
sequence while preserving dinucleotide frequencies
<http://clavius.bc.edu/~clotelab/RNAdinucleotideShuffle/ShuffleCodeParts/altschulEriksonDinuclShuffle.txt>

NOTE: One cannot use function "count(s,word)" to count the number
of occurrences of dinucleotide word in string s, since the built-in
function counts only non-overlapping words, presumably in a left to
right fashion.
"""

__author__ = ("Original implementation by Peter Clote (Oct 2003), "
              "organized by Jean-Frédéric (March & November 2012)")


# P. Clote, Oct 2003
# NOTE: One cannot use function "count(s,word)" to count the number
# of occurrences of dinucleotide word in string s, since the built-in
# function counts only nonoverlapping words, presumably in a left to
# right fashion.


import sys
import random
import copy
import itertools
from collections import Counter


class Shuffle:

    """Provide dinucleotides shuffles of a given sequence.

    nucl_count
        A count of each nucleotide, using a dictionary
        Example: {'A': 3, 'C': 1, 'U': 2, 'G': 2}

    _dinucl_count
        A count of each dinucleotide, using a dictionary
        Example: {'A': {'A': 0, 'C': 0, 'U': 1, 'G': 1},
                  'C': {'A': 0, 'C': 0, 'U': 1, 'G': 0},
                  'U': {'A': 1, 'C': 0, 'U': 0, 'G': 1},
                  'G': {'A': 1, 'C': 1, 'U': 0, 'G': 0}}

    edgelists_dict
        A representation of the edge-lists of the multigraph
        (Note: this only serves to be copied in dinuclShuffle())

    """

    def __init__(self):
        self.possible_nucleotides_list = ["A", "C", "G", "U"]
        self.sequence_composition_list = []
        self.nucl_count = None
        self._dinucl_count = None
        self.edgelists_dict = None
        self.firstCh = None
        self.lastCh = None
        self.sequence_length = None
        self.sequence = None

    def init_with_sequence(self, sequence):
        """Initialise the Shuffle with a given sequence."""
        self.sequence = sequence.upper().replace("T", "U")
        self.sequence_length = len(self.sequence)
        self.firstCh = self.sequence[0]  # start with first letter of s
        self.lastCh = self.sequence[-1]
        try:
            self.computeCountAndLists()
        except KeyError:
            msg = ("Source sequence contains unexpected character: "
                   "only acceptable are %s"
                   % ' '.join(self.possible_nucleotides_list))
            raise ShuffleException(msg)
        self.sequence_composition_list = [k for k, v
                                          in self.nucl_count.items()
                                          if v > 0]

    def init_with_data(self, dinucleotides, firstCh, lastCh):
        """Initialise the Shuffle with with a dinucleotides "plain" count."""
        self._dinucl_count = DataStructuresConversion.plain2intricate(dinucleotides)
        self.sequence_length = sum(dinucleotides.values()) + 1
        self.firstCh = firstCh
        self.lastCh = lastCh
        self.edgelists_dict = DataStructuresConversion.intricate2edgelist(self._dinucl_count)
        self.sequence_composition_list = list(set(itertools.chain(*self.edgelists_dict.values())))
        self.check_valid_input(dinucleotides, firstCh, lastCh)

    @staticmethod
    def check_valid_input(dinucleotides, firstCh, lastCh):
        """Check a Shuffle data input (in "plain" format) for validity."""
        outgoing_counter = Counter()
        incoming_counter = Counter()
        for key, value in dinucleotides.items():
            first = key[0]
            second = key[1]
            if first is not second:  # Ignoring loops
                outgoing_counter[first] += value
                incoming_counter[second] += value
        outgoing_counter.subtract(incoming_counter)
        result = (outgoing_counter.values() == [1, -1] or
                  len(list(outgoing_counter.elements())) is 0)
	return result

    def computeCountAndLists(self):
        """Initialise lists and mono- and dinucleotide dictionaries.

        It is called at initialisation time when using a sequence as base.

        """
        #WARNING: Use of function count(s,'UU') returns 1 on word UUU
        #since it apparently counts only non-overlapping words UU
        #For this reason, we work with the indices.
        list_dict = {}  # list_dict is a dictionary of lists
        nucl_count = {}    # empty dictionary
        dinucl_count = {}  # empty dictionary

        for nucleotide in self.possible_nucleotides_list:
            list_dict[nucleotide] = []
            nucl_count[nucleotide] = 0
            dinucl_count[nucleotide] = {}
            for nucleotide_bis in self.possible_nucleotides_list:
                dinucl_count[nucleotide][nucleotide_bis] = 0

        #Compute count and lists
        nucl_count[self.firstCh] = 1
        nuclTotal = 1
        dinuclTotal = 0
        for index in range(self.sequence_length - 1):
            first_nucleotide = self.sequence[index]
            second_nucleotide = self.sequence[index + 1]
            list_dict[first_nucleotide].append(second_nucleotide)
            nucl_count[second_nucleotide] += 1
            nuclTotal += 1
            dinucl_count[first_nucleotide][second_nucleotide] += 1
            dinuclTotal += 1
        assert (nuclTotal == self.sequence_length)
        assert (dinuclTotal == self.sequence_length - 1)
        self.nucl_count = nucl_count
        self._dinucl_count = dinucl_count
        self.edgelists_dict = list_dict

    def build_last_edge_graph(self):
        """Build a random last-edge graph respecting the Eulerian condition.

        Build a last-edge graph Z until it satisfies the Eulerian condition,
        ie if all its vertices are connected in Z to the last vertex.

        """
        ok = False
        while not ok:
            edgeList = self._build_random_last_edge_graph()
            ok = GraphHelper.is_connected_to_last(edgeList,
                                                  self.sequence_composition_list,
                                                  self.lastCh)
        return edgeList

    def _build_random_last_edge_graph(self):
        """Construct randomly a last-edge graph (as a list of last edges)."""
        dinucl_count = copy.deepcopy(self._dinucl_count)
        edgeList = []
        for nucleotide in self.sequence_composition_list:
            if nucleotide != self.lastCh:
                choosen_edge = GraphHelper.chooseEdge(nucleotide, dinucl_count)
                edgeList.append([nucleotide, choosen_edge])
        return edgeList

    def dinuclShuffle(self):
        """Provide a shuffled sequence."""
        last_edge_List = self.build_last_edge_graph()
        current_edgelists = copy.deepcopy(self.edgelists_dict)
        #remove last edges from each vertex list, shuffle, then add back
        #the removed edges at end of vertex lists.
        for [x, y] in last_edge_List:
            current_edgelists[x].remove(y)
        for nucleotide in self.sequence_composition_list:
            GraphHelper.shuffleEdgeList(current_edgelists[nucleotide])
        for [x, y] in last_edge_List:
            current_edgelists[x].append(y)

        #construct the eulerian path
        L = [self.firstCh]
        prevCh = self.firstCh
        for _ in range(self.sequence_length - 2):
            ch = current_edgelists[prevCh][0]
            L.append(ch)
            del current_edgelists[prevCh][0]
            prevCh = ch
        L.append(self.lastCh)
        return "".join(L)

    def getNshuffles(self, times):
        """Get a given amount of shuffles."""
        for _ in range(times):
            yield self.dinuclShuffle()


class ShuffleException(Exception):

    """Exception for Shuffle operations."""

    pass


class DataStructuresConversion:

    """Helper methods to convert data structures."""

    @staticmethod
    def intricate2plain(intricate):
        """Convert an "intricate" nucleotides data structure to a "plain" one.

        {'A': {'A': 0, 'C': 0, 'U': 1, 'G': 1},
        'C': {'A': 0, 'C': 0, 'U': 1, 'G': 0},
        'U': {'A': 1, 'C': 0, 'U': 0, 'G': 1},
        'G': {'A': 1, 'C': 1, 'U': 0, 'G': 0}}

        """
        result = {}
        for first, simple in intricate.items():
            for second, value in simple.items():
                dinucl = first + second
                result[dinucl] = value
        return result

    @staticmethod
    def plain2intricate(plain):
        """Convert a "plain" nucleotides data structure to an "intricate" one.

        plain = {'AA': 2, 'AT': 2, 'AG': 2, 'AC': 5,
                  'TA': 2, 'TT': 2, 'TG': 2, 'TC': 5,
                  'GA': 2, 'GT': 2, 'GG': 2, 'GC': 5,
                  'CA': 2, 'CT': 2, 'CG': 2, 'CC': 5}

        intricate = {'A': {'A': 2, 'C': 5, 'T': 2, 'G': 2},
                     'C': {'A': 2, 'C': 5, 'T': 2, 'G': 2},
                     'G': {'A': 2, 'C': 5, 'T': 2, 'G': 2},
                     'T': {'A': 2, 'C': 5, 'T': 2, 'G': 2}}

        """
        result = {}
        for key in set([x[0:1] for x in plain.keys()]):
            result[key] = {}
        for dinucl, value in plain.items():
            first, second = dinucl[:1], dinucl[1:]
            result[first][second] = value
        return result

    @staticmethod
    def intricate2edgelist(double):
        """Convert a "double" nucleotides data structure to an edge-list.

        {'A': {'A': 0, 'C': 0, 'U': 1, 'G': 1},
          'C': {'A': 0, 'C': 0, 'U': 1, 'G': 0},
          'U': {'A': 1, 'C': 0, 'U': 0, 'G': 1},
          'G': {'A': 1, 'C': 1, 'U': 0, 'G': 0}}

        edgelists_dict = {'A': ['U', 'U'], 'C': [],
                          'U': ['U', 'A', 'U', 'A'],
                          'G': []}

        """
        result = {}
        for first, simple in double.items():
            result[first] = []
            for second, value in simple.items():
                for _ in range(value):
                    result[first].append(second)
        return result


class GraphHelper:

    """Methods used when building the eulerian walk."""

    @staticmethod
    def chooseEdge(nucleotide, dinuclCnt):
        """Choose an edge, for the given vertex.

        Args:
            nucleotide: The vertex to use
            dinuclCnt: A count of each dinucleotide, as a dictionary

        Returns:
            The randomly choosed vertex/nucleotide

        """
        z = random.random()
        denom = sum(dinuclCnt[nucleotide].values())
        numerator = dinuclCnt[nucleotide]['A']
        if z < float(numerator) / float(denom):
            dinuclCnt[nucleotide]['A'] -= 1
            return 'A'
        numerator += dinuclCnt[nucleotide]['C']
        if z < float(numerator) / float(denom):
            dinuclCnt[nucleotide]['C'] -= 1
            return 'C'
        numerator += dinuclCnt[nucleotide]['G']
        if z < float(numerator) / float(denom):
            dinuclCnt[nucleotide]['G'] -= 1
            return 'G'
        dinuclCnt[nucleotide]['U'] -= 1
        return 'U'

    @staticmethod
    def is_connected_to_last(edgeList, vertices, lastCh):
        """Determine whether all vertices are connected to a given vertex.

        Args:
            edgeList: a collection of edges
            vertices: a collection of vertices
            lastCh: the last vertex to use

        Returns:
            True if all vertices are connected to the given vertex,
            false otherwise

        """
        D = {}
        for vertex in vertices:
            D[vertex] = 0
        for edge in edgeList:
            a = edge[0]
            b = edge[1]
            if b == lastCh:
                D[a] = 1
        for _ in range(2):
            for edge in edgeList:
                a = edge[0]
                b = edge[1]
                if D[b] == 1:
                    D[a] = 1
        for vertex in vertices:
            if vertex != lastCh and D[vertex] == 0:
                return False
        return True

    @staticmethod
    def shuffleEdgeList(L):
        """Shuffle the given edge list.

        Returns: the given collection, randomly shuffled

        """
        length = len(L)
        barrier = length
        for _ in range(length - 1):
            z = int(random.random() * barrier)
            tmp = L[z]
            L[z] = L[barrier - 1]
            L[barrier - 1] = tmp
            barrier -= 1
        return L


def shuffleNtimes(times, sequence):
    """Shuffle n times the given sequence.

    times - the number of times to shuffle the sequence
    sequence - the sequence to shuffle

    """
    myShuffle = Shuffle()
    myShuffle.init_with_sequence(sequence)
    for _ in range(times):
        yield myShuffle.dinuclShuffle()


def main(times, sequence_file):
    """Main method."""
    with open(sequence_file, 'r') as f:
        s = f.read().rstrip()
    myShuffle = Shuffle(s)
    for _ in range(times):
        print myShuffle.dinuclShuffle()


if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(int(sys.argv[1]), sys.argv[2])
    else:
        print "Usage: %s nb RNAsequence.raw" % sys.argv[0]
        sys.exit(1)
#        main(10,'dummy.raw')
#        print "Usage: %s nb RNAsequence.raw"
#        sys.exit(1)
