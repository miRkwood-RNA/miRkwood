"""Test Task.py."""

import unittest
from shuffles.Task import chop


class TestTask(unittest.TestCase):

    """Testing Task methods."""

    def test_chop(self):
        """Test chop()."""
        input_chop = xrange(30)
        result = chop(input_chop, 10)
        expected = [range(10), range(10, 20), range(20, 30)]
        self.assertListEqual(list(result), expected)


if __name__ == "__main__":
    unittest.main()
