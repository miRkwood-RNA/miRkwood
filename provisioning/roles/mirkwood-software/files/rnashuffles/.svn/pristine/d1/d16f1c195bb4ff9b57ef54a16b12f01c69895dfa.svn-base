"""Test Shuffle.py."""

import unittest
from shuffles.Shuffle import GraphHelper
from shuffles.Shuffle import Shuffle
from shuffles.Shuffle import DataStructuresConversion as DSC
from shuffles.Shuffle import ShuffleException


class TestGraphFuncs(unittest.TestCase):

    def test_shuffleEdgeList(self):
        """Testing shuffleEdgeList."""
        L = ["A", "B", "C", "D"]
        self.assertEqual(sorted(L),
                         sorted(GraphHelper.shuffleEdgeList(L)))
        L2 = ['A', 'B']
        self.assertEqual(sorted(L2),
                         sorted(GraphHelper.shuffleEdgeList(L2)))
        L3 = ['A']
        self.assertEqual(L3, GraphHelper.shuffleEdgeList(L3))


class TestShuffle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.intricate = {'A': {'A': 0, 'C': 0, 'U': 3, 'G': 0},
                         'C': {'A': 0, 'C': 0, 'U': 0, 'G': 0},
                         'U': {'A': 2, 'C': 0, 'U': 0, 'G': 0},
                         'G': {'A': 0, 'C': 0, 'U': 0, 'G': 0}}
        cls.plain = {'AA': 0, 'AC': 0, 'AG': 0, 'AU': 3,
                     'CC': 0, 'CA': 0, 'CG': 0, 'CU': 0,
                     'GU': 0, 'GG': 0, 'GA': 0, 'GC': 0,
                     'UU': 0, 'UC': 0, 'UG': 0, 'UA': 2}
        cls.sequence = "AUAUAU"
        super(TestShuffle, cls).setUpClass()

    def test_check_valid_input(self):
        """Testing check_valid_input."""
        L1 = {'AA': 0, 'AC': 2,
              'CA': 1, 'CC': 0}
        L2 = {'AA': 0, 'AC': 1,
              'CA': 1, 'CC': 0}
        L3 = {'AA': 0, 'AC': 2,
              'CA': 0, 'CC': 0}
        L4 = {'AA': 5, 'AC': 1,
              'CA': 1, 'CC': 0}
        self.assertTrue(Shuffle.check_valid_input(L1, None, None))
        self.assertTrue(Shuffle.check_valid_input(L2, None, None))
        self.assertFalse(Shuffle.check_valid_input(L3, None, None))
        self.assertTrue(Shuffle.check_valid_input(L4, None, None))

    def test_init_with_sequence(self):
        """Test init_with_sequence."""
        s = Shuffle()
        s.init_with_sequence(self.sequence)
        self.assertEquals(s._dinucl_count,
                          self.intricate)
        self.assertEqual(s.sequence_length, 6)

    def test_init_with_sequence_with_N(self):
        """Test init_with_sequence with something else than ATGC."""
        s = Shuffle()
        with self.assertRaises(ShuffleException):
            s.init_with_sequence('atgcatgcnatgcatgc')

    def test_shuffle_sequence(self):
        """Test shuffling based on a sequence."""
        s = Shuffle()
        s.init_with_sequence(self.sequence)
        self.assertListEqual(list(s.getNshuffles(5)),
                             [self.sequence] * 5)

    def test_init_with_data(self):
        """Test init_with_data."""
        s = Shuffle()
        s.init_with_data(self.plain, "A", "U")
        self.assertEquals(s._dinucl_count,
                          self.intricate)
        self.assertEqual(s.sequence_length, 6)
        self.assertListEqual(s.sequence_composition_list,
                             ['A', 'U'])

    def test_shuffle_data(self):
        """Test shuffling based on data."""
        s = Shuffle()
        s.init_with_data(self.plain, "A", "U")
        self.assertListEqual(list(s.getNshuffles(5)),
                             [self.sequence] * 5)


class TestDataStructuresConversion(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.plain = {'AA': 2, 'AT': 2, 'AG': 2, 'AC': 5,
                     'TA': 2, 'TT': 2, 'TG': 2, 'TC': 5,
                     'CA': 2, 'CT': 2, 'CG': 2, 'CC': 5,
                     'GA': 2, 'GT': 2, 'GG': 2, 'GC': 5}

        cls.intricate = {'A': {'A': 2, 'C': 5, 'G': 2, 'T': 2},
                         'C': {'A': 2, 'C': 5, 'G': 2, 'T': 2},
                         'G': {'A': 2, 'C': 5, 'G': 2, 'T': 2},
                         'T': {'A': 2, 'C': 5, 'G': 2, 'T': 2}}
        super(TestDataStructuresConversion, cls).setUpClass()

    def test_plain2intricate(self):
        """Test plain2intricate."""
        self.assertEqual(DSC.plain2intricate(self.plain),
                         self.intricate)

    def test_intricate2plain(self):
        """Test intricate2plain."""
        self.assertEqual(DSC.intricate2plain(self.intricate),
                         self.plain)

    def test_consistency(self):
        """Test intricate2plain and plain2intricate consistency."""
        self.assertEqual(DSC.plain2intricate(
                         DSC.intricate2plain(self.intricate)),
                         self.intricate)
        self.assertEqual(DSC.intricate2plain(
                         DSC.plain2intricate(self.plain)),
                         self.plain)
