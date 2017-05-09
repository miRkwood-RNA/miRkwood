"""Test RNAfold.py."""

import unittest
import re
import os
from StringIO import StringIO
from shuffles.RNAfold import (RNAfold, RNAfoldParser, RNAFoldParserBis,
                              get_mfe, get_mfe_and_structure,
                              get_several_mfe_and_structure)


def skipIfRNAfoldUnavailable(some_func):
    """Skip the test if RNAfold is not available."""
    def decorator(func):
        """The custom  unittest.skipIf decorator."""
        condition = True
        try:
            rnafold = RNAfold()
            condition = False
        except:
            pass
        reason = "Executable for RNAfold is not available."
        return unittest.skipIf(condition, reason)(func)
    return decorator(some_func)


class TestRNAFoldBase(unittest.TestCase):

    """Base class for testing RNAfold related methods."""

    @classmethod
    def setUpClass(cls):
        cls.query_sequence = "GAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTT"
        cls.reference = "...(((..........))).................... (-1.7)"
        cls.reference2 = """
length = 39
GAAGGGCUAAUUCACUCCCAACGAAGACAAGAUAUCCUU
...(((..........)))....................
 minimum free energy =  -1.70 kcal/mol
"""
        cls.expected = (-1.7, '...(((..........)))....................')


class TestRNAFoldParsing(TestRNAFoldBase):

    """Testing RNAfold parsing related methods."""

    def test_structure_re(self):
        """Test structure_re definition."""
        structure_match = re.search(RNAfoldParser().structure_re,
                                    self.reference)
        self.assertIsNotNone(structure_match,
                             msg=('Regular expression for structure line '
                                  'in RNAFold parser does not match'))

    def test_structure(self):
        """Test _structure()."""
        structure_match = re.search(RNAfoldParser().structure_re,
                                    self.reference)
        result = RNAfoldParser()._structure(structure_match)
        self.assertEqual(result,
                         self.expected,
                         msg=('Structure parsing in RNAFold parser '
                              'did not return the correct values'))

    def test_parse(self):
        """Test parse()."""
        string_to_parse = self.query_sequence + "\n" + self.reference + "\n"
        result = RNAfoldParser().parse(StringIO(string_to_parse))
        self.assertEqual(result, self.expected)


class TestRNAFoldParsingBis(TestRNAFoldBase):

    """Testing RNAfold parsing related methods."""

    def test_structure_re(self):
        """Test structure_re definition."""
        structure_match = re.search(RNAFoldParserBis().structure_re,
                                    self.reference2.split('\n')[3])
        self.assertIsNotNone(structure_match,
                             msg=("Regular expression for structure line "
                                  "in RNAFold parser does not match"))

    def test_structure(self):
        """Test _structure()."""
        structure_match = re.search(RNAFoldParserBis().structure_re,
                                    self.reference2.split('\n')[3])
        result = RNAFoldParserBis()._structure(structure_match)
        self.assertEqual(result,
                         self.expected[1],
                         msg=("Structure parsing in RNAFold parser "
                              "did not return the correct values"))

    def test_energy_re(self):
        """Test energy_re definition."""
        energy_match = re.search(RNAFoldParserBis().energy_re,
                                 self.reference2.split('\n')[4])

        self.assertIsNotNone(energy_match,
                             msg=("Regular expression for energy line "
                                  "in RNAFold parser does not match"))

    def test_energy(self):
        """Test _energy() - Energy parsing in RNAFoldBis parser."""
        energy_match = re.search(RNAFoldParserBis().energy_re,
                                 self.reference2.split('\n')[4])
        result = RNAFoldParserBis()._energy(energy_match)
        self.assertEqual(result,
                         self.expected[0])

    def test_parse(self):
        """Test parse()."""
        string_to_parse = self.reference2
        result = RNAFoldParserBis().parse(StringIO(string_to_parse))
        self.assertEqual(result, self.expected)


class TestRNAFoldWrapper(TestRNAFoldBase):

    """Testing RNAfold wrapper related methods."""

    @skipIfRNAfoldUnavailable
    def test_executable(self):
        """Test RNAFold() _executable() method."""
        rnafold = RNAfold()
        self.assertEqual(rnafold._executable(), "RNAfold")

    @skipIfRNAfoldUnavailable
    def test_run(self):
        """Test RNAFold() run() method."""
        rnafold = RNAfold()
        result = rnafold.run(self.query_sequence, verbose=True)
        self.assertEqual(result, self.expected)

    @skipIfRNAfoldUnavailable
    def test_run_interactive(self):
        """Test RNAFold() run_interactive() method."""
        rnafold = RNAfold()
        sequences = [self.query_sequence] * 3
        expected = [self.expected] * 3
        result = list(rnafold.run_interactive(sequences))
        self.assertListEqual(result, expected)


class TestRNAFold(TestRNAFoldBase):

    """Testing RNAfold public methods."""

    @skipIfRNAfoldUnavailable
    def test_get_mfe_and_structure(self):
        """Test get_mfe_and_structure()."""
        result = get_mfe_and_structure(self.query_sequence)
        self.assertEqual(result, self.expected)

    @skipIfRNAfoldUnavailable
    def test_get_mfe(self):
        """Test get_mfe()."""
        result = get_mfe(self.query_sequence)
        self.assertEqual(result, self.expected[0])

    @skipIfRNAfoldUnavailable
    def test_get_several_mfe_and_structure(self):
        """Test get_several_mfe_and_structure()."""
        result = get_several_mfe_and_structure([self.query_sequence] * 10)
        self.assertListEqual(list(result), [self.expected] * 10)
