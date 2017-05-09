"""Test Shuffle.py."""

import unittest
from shuffles.ParameterGenerator import read_dataset_from_file_old, ParameterGenerator
from shuffles.ParameterGenerator import read_dataset_from_file, _parse_dataset_line
from StringIO import StringIO
from collections import Counter


class TestParameterGenerator(unittest.TestCase):

    """Base class for testing ParameterGenerator related methods."""

    @classmethod
    def setUpClass(cls):
        """Set up the class with a few useful things."""
        cls.dinucleotides = ['AC', 'AG', 'AU',
                             'CA', 'CG', 'CU',
                             'GA', 'GC', 'GU',
                             'UA', 'UC', 'UG',
                             'AA', 'CC', 'GG', 'UU']
        cls.expected = [Counter({'AA': 8, 'AC': 3, 'AG': 3, 'CC': 12, 'CA': 3,
                                 'CG': 3, 'GU': 9, 'GG': 12, 'AU': 3, 'GA': 3,
                                 'CU': 7, 'UU': 9, 'GC': 3, 'UG': 9, 'UA': 3,
                                 'UC': 7}),
                        Counter({'AA': 12, 'GG': 12, 'UU': 9, 'GU': 9, 'UG': 9,
                                 'CC': 8, 'UC': 7, 'CU': 7, 'AC': 3, 'AG': 3,
                                 'CA': 3, 'CG': 3, 'AU': 3, 'GA': 3, 'GC': 3,
                                 'UA': 3})]


class TestParameterGeneratorParsing(TestParameterGenerator):

    """Testing ParameterGenerator methods related to parsing files."""

    @classmethod
    def setUpClass(cls):
        """Set up the class with a few useful things."""
        super(TestParameterGeneratorParsing, cls).setUpClass()

    def test_parse_dataset_line(self):
        """Test _parse_dataset_line()."""
        input_strings = [("[('AA', 8), ('AC', 3), ('AG', 3), ('AU', 3),"
                         "('CA', 3), ('CC', 12), ('CG', 3), ('CU', 7),"
                         "('GA', 3), ('GC', 3), ('GG', 12), ('GU', 9),"
                         "('UA', 3), ('UC', 7), ('UG', 9), ('UU', 9)]"),
                         ("[('AA', 12), ('AC', 3), ('AG', 3), ('AU', 3),"
                         "('CA', 3), ('CC', 8), ('CG', 3), ('CU', 7),"
                         "('GA', 3), ('GC', 3), ('GG', 12), ('GU', 9),"
                         "('UA', 3), ('UC', 7), ('UG', 9), ('UU', 9)]")]
        results = [_parse_dataset_line(x) for x in input_strings]
        self.assertListEqual(results, self.expected)

    def test_read_dataset_old_from_file(self):
        """Test read_dataset_from_file_old()."""
        input_string = ("{"
                        "'AA': 8, 'AC': 3, 'AG': 3, 'CC': 12, 'CA': 3, "
                        "'CG': 3, 'GU': 9, 'GG': 12,'AU': 3, 'GA': 3, "
                        "'CU': 7, 'UU': 9, 'GC': 3, 'UG': 9, 'UA': 3, 'UC': 7"
                        "}\n{"
                        "'AA': 12, 'AC': 3, 'AG': 3, 'CC': 8, 'CA': 3, "
                        "'CG': 3, 'GU': 9, 'GG': 12, 'AU': 3, 'GA': 3, "
                        "'CU': 7, 'UU': 9, 'GC': 3, 'UG': 9, 'UA': 3, 'UC': 7"
                        "}\n")
        result = list(read_dataset_from_file_old(StringIO(input_string)))
        self.assertListEqual(result, self.expected)

    def test_read_dataset_from_file(self):
        """Test read_dataset_from_file()."""
        input_string = ("[('AA', 8), ('AC', 3), ('AG', 3), ('AU', 3),"
                        "('CA', 3), ('CC', 12), ('CG', 3), ('CU', 7),"
                        "('GA', 3), ('GC', 3), ('GG', 12), ('GU', 9),"
                        "('UA', 3), ('UC', 7), ('UG', 9), ('UU', 9)]"
                        "\n"
                        "[('AA', 12), ('AC', 3), ('AG', 3), ('AU', 3),"
                        "('CA', 3), ('CC', 8), ('CG', 3), ('CU', 7),"
                        "('GA', 3), ('GC', 3), ('GG', 12), ('GU', 9),"
                        "('UA', 3), ('UC', 7), ('UG', 9), ('UU', 9)]")
        result = list(read_dataset_from_file(StringIO(input_string)))
        self.assertListEqual(result, self.expected)


class TestParameterGeneratorInternals(TestParameterGenerator):

    """Testing ParameterGenerator methods related to internal coding."""

    @classmethod
    def setUpClass(cls):
        """Set up the class with a few useful things."""
        super(TestParameterGeneratorInternals, cls).setUpClass()
        cls.partial = Counter({'AC': 2, 'AG': 2, 'AU': 2,
                               'CA': 2, 'CG': 2, 'CU': 2,
                               'GA': 2, 'GC': 2})
        cls.complete = Counter({'AC': 2, 'AG': 2, 'AU': 2,
                                'CA': 2, 'CG': 2, 'CU': 2,
                                'GA': 2, 'GC': 2,
                                'GU': 30, 'UG': 30,  'UC': 2, 'UA': 2})

    def test_already_affected(self):
        """Testing _already_affected()."""
        parameter_generator = ParameterGenerator(None, None, None)
        parameter_generator.values = self.expected[0]
        self.assertEqual(parameter_generator._already_affected(0), 0)
        self.assertEqual(parameter_generator._already_affected(1), 3)
        self.assertEqual(parameter_generator._already_affected(2), 6)
        self.assertEqual(parameter_generator._already_affected(16),
                         sum(self.expected[0].values()))

    def test_compute_end_of_set(self):
        """Testing _compute_end_of_set()."""
        parameter_generator = ParameterGenerator(None, None, None)
        parameter_generator.values = self.partial
        parameter_generator._compute_end_of_set(80)
        self.assertEqual(parameter_generator.values, self.complete)

class TestParameterGeneratorGeneration(TestParameterGenerator):

    """Testing ParameterGenerator methods related to parameter generation."""

    def test_make_with_uniform_frequency(self):
        """Test make_with_uniform_frequency()."""
        expected_result = Counter()
        value = 5
        for din in self.dinucleotides:
            expected_result[din] = value
        parameter_generator = ParameterGenerator(None, None, None)
        result = parameter_generator.make_with_uniform_frequency(value)
        self.assertEqual(result, expected_result)
        self.assertEqual(value * 16, sum(result.values()))

    def test_make_with_enriched_AU(self):
        """Test make_with_enriched_AU()."""
        expected_result = Counter()
        value = 5
        for din in self.dinucleotides:
            expected_result[din] = value
        for din in ['AA', 'AU', 'UU', 'UA']:
            expected_result[din] = 2 * value
        parameter_generator = ParameterGenerator(None, None, None)
        result = parameter_generator.make_with_enriched_AU(value)
        self.assertEqual(result, expected_result)
        self.assertEqual(value * 20, sum(result.values()))

    def test_make_with_enriched_CG(self):
        """Test make_with_enriched_CG()."""
        expected_result = Counter()
        value = 5
        for din in self.dinucleotides:
            expected_result[din] = value
        for din in ['CC', 'CG', 'GC', 'GG']:
            expected_result[din] = 2 * value
        parameter_generator = ParameterGenerator(None, None, None)
        result = parameter_generator.make_with_enriched_CG(value)
        self.assertEqual(result, expected_result)
        self.assertEqual(value * 20, sum(result.values()))

if __name__ == "__main__":
    unittest.main()
