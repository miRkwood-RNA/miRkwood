"""
RNA secondary structure folding

RNAFold reads RNA sequences, calculates their minimum free energy (mfe)
and returns the mfe structure in bracket notation and its free energy.
"""

from subprocess import Popen, PIPE
import pexpect
import os.path
import os
import re
from StringIO import StringIO
import sys


class RNAfoldParser():

    """A parser for RNAfold output."""

    def __init__(self):
        """Constructor."""
        self.structure_re = re.compile(r"""
        ^                            # Begin of line
        (?P<structure>[\.\(\)]+?)    # Either point or parenthesis, captured.
        \s+?                         # Whitespace
        \(                           # Opening parenthesis
        \s*?                         # Possible whitespace
        (?P<energy>[-.\d]*?)         # A negative number
        \s*?                         # Possible whitespace
        \)                           # Closing parenthesis
        \s*$                         # Whitespace until end of line
        """, re.X)

    @staticmethod
    def _structure(structure_match):
        """Parse a structure line.and return the (structure, energy) couple."""
        energy = float(structure_match.group('energy'))
        structure = structure_match.group('structure')
        return (energy, structure)

    def parse(self, handle):
        """Parse the output file and return the (structure, energy) couple."""
        for line in handle:
            if line.strip() == "":
                continue
            else:
                structure_match = re.search(self.structure_re, line)
                if structure_match:
                    return self._structure(structure_match)
                else:
                    continue


class RNAFoldParserBis():

    """Parser of the RNAFold output, when used in interactive mode."""

    def __init__(self):
        """Constructor."""
        self.structure_re = re.compile(r"""
        ^                            # Begin of line
        (?P<structure>[\.\(\)]+?)    # Either point or parenthesis, captured.
        \s*?$                        # Whitespace until end of line
        """, re.X)
        self.energy_re = re.compile(r"""
        ^                            # Begin of line
        [\w\s=]*                     # Either whitespace, letters or equal sign
        (?P<energy>[-.\d]*?)         # A negative number
        \skcal\/mol                  # Whitespace then "kcal/mol"
        \s*?$                        # Whitespace until end of line
        """, re.X)

    @staticmethod
    def _structure(structure_match):
        """Parse a structure line and return the structure."""
        structure = structure_match.group('structure')
        return structure

    @staticmethod
    def _energy(energy_match):
        """Parse an energy line and return the energy."""
        energy = float(energy_match.group('energy'))
        return energy

    def parse(self, handle):
        """Parse the output file and return the (structure, energy) couple."""
        for line in handle:
            if line.strip() == "":
                continue
            else:
                structure_match = re.search(self.structure_re, line)
                if structure_match:
                    structure = self._structure(structure_match)
                else:
                    energy_match = re.search(self.energy_re, line)
                    if energy_match:
                        energy = self._energy(energy_match)
                        return (energy, structure)
                    else:
                        continue


class RNAFoldException(Exception):

    """Exception for RNAfold execution."""

    pass


class RNAFoldUnavailableException(RNAFoldException):

    """Exception for RNAfold unavailability."""

    pass

class RNAFoldExecutionException(RNAFoldException):

    """Exception for RNAfold execution failure."""

    pass


class RNAfold():

    """RNA secondary structure folding."""

    def __init__(self):
        """Constructor."""
        self.program_name = "RNAfold"
        self.program_dir = "/usr/bin/"
        executable = self.program_name
        if not self._is_tool([executable, '--version']):
            executable = os.path.join(self.program_dir, self.program_name)
            if not self._is_tool([executable, '--version']):
                raise RNAFoldUnavailableException()
        self.executable = executable

    @staticmethod
    def _is_tool(command):
        """Return whether the given command is an executable."""
        try:
            devnull = open(os.devnull)
            Popen(command, stdout=devnull, stderr=devnull).communicate()
        except OSError, error:
            if error.errno == os.errno.ENOENT:
                return False
        return True

    def _executable(self):
        """Return the program executable."""
        return self.executable

    def run(self, sequence, verbose=False):
        """Run RNAfold and return the (energy, structure) couple."""
        child_process = Popen([self._executable(), '--noPS'], close_fds=True,
                              stdin=PIPE,
                              stdout=PIPE,
                              stderr=PIPE)
        stdout_str, stderr_str = child_process.communicate(input=sequence)
        return_code = child_process.returncode
        if(verbose):
            print str(stdout_str)
        if return_code is not 0:
            raise RNAFoldExecutionException(stderr_str)
        child_process.stdin.close()
        child_process.stdout.close()
        child_process.stderr.close()
        parser = RNAfoldParser()
        return parser.parse(StringIO(stdout_str))

    def run_interactive(self, sequences):
        """Run RNAfold on several sequences.

        Yield the (energy, structure) couple for each sequence.

        """

        parser = RNAFoldParserBis()
        child_process = pexpect.spawn(self._executable(), ['--noPS'])
        child_process.delaybeforesend = 0
        child_process.setecho(False)
        for sequence in sequences:
            child_process.expect('...8\r\n')
            child_process.sendline(sequence)
            child_process.expect('\r\n\r\n')
            yield (parser.parse(StringIO(child_process.before)))
        child_process.terminate()


def get_several_mfe_and_structure(sequences):
    """Return the mfe and structures for the given sequences."""
    wrapper = RNAfold()
    return wrapper.run_interactive(sequences)


def get_mfe(sequence):
    """Return the mininum free energy (mfe) for the given sequence."""
    (energy, _) = get_mfe_and_structure(sequence)
    return energy


def get_mfe_and_structure(sequence, verbose=False):
    """Return the mininum free energy and structure for the given sequence."""
    wrapper = RNAfold()
    return wrapper.run(sequence=sequence, verbose=verbose)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        print get_mfe_and_structure(sys.argv[1])
    else:
        print "Usage: %s sequence" % sys.argv[0]
        sys.exit(1)
