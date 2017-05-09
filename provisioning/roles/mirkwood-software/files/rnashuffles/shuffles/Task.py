#! /usr/bin/env python
# -*- coding: latin-1 -*-
"""Generating mfe distributions for all possible dinucleotides 16-uplets
admitting an Eulerian walk.

(December 2012)
"""
__author__ = "Concept by H�l�ne Touzet, implementation by Jean-Fr�d�ric."


from ParameterGenerator import ParameterGenerator, _parse_dataset_line
from Shuffle import Shuffle
import RNAfold
import multiprocessing
import sys
import itertools


class TaskException(Exception):

    """Exception for Task execution."""

    pass


class Task():

    """Generating mfe distribution."""

    def __init__(self, nb_shuffles=1000, processes=1,
                 generator_step=4, sequence_length=100):
        """Constructor."""
        self.nb_shuffles = nb_shuffles
        self.edge_nucleotide = 'U'
        self.processes = processes
        self.sequence_length = sequence_length
        self.generator_step = generator_step

    def info(self):
        """Return some information about the Task parameters."""
        result = "Sequence %s" % (self.sequence_length)
        result += "\nParameters set generator step: %s" % self.generator_step
        result += "\nShuffles: %s" % self.nb_shuffles
        result += "\nUsing %s processes" % self.processes
        return result

    def execute_threaded_bis(self):
        """Execute the task alternatively."""
        alpha = int(self.sequence_length / 32)
        beta = int(self.sequence_length / 8)
        generator = ParameterGenerator(self.sequence_length, alpha, beta)
        generator.step = self.generator_step
        filename = 'task-results-%s-%s.txt' % (self.sequence_length,
                                               self.generator_step)
        parameters_sets = generator.generate_dinucleotides_sets_with_loops()
        self.process_parameters_set_iterable(parameters_sets, filename)

    def process_parameters_set_iterable_old(self, parameters_sets, filename):
        """Process an iterable of parameters set (the old way)."""
        with open(filename, 'w') as f:
            for dinucleotides_set in parameters_sets:
                shuffler = Shuffle()
                shuffler.init_with_data(dinucleotides_set,
                                        self.edge_nucleotide, self.edge_nucleotide)
                f.write(str(dict(dinucleotides_set)) + "\n")
                try:
                    workers = multiprocessing.Pool(processes=self.processes)
                    mfe_list = workers.map_async(RNAfold.get_mfe,
                                                 shuffler.getNshuffles(self.nb_shuffles))
                    f.write(str(sorted(mfe_list.get())) + "\n")
                    workers.terminate()
                except Exception, e:
                    raise TaskException(e)

    def process_parameters_set_iterable(self, parameters_sets, filename):
        """Process an iterable of parameters set.

        - Build an iterable with the parameters_sets and the shuffle data
        - Chop this iterable into smaller iterables chunks
        - Using a pool of N workers, map the parallelised worker_unit
          to each chunk - meaning we apply the code to each parameters set
        - Write the worker_unit result on disk under the given filename

        """

        chop_size = 500
        shuffle_data = (self.edge_nucleotide, self.nb_shuffles)
        iterable = itertools.izip(parameters_sets,
                                  itertools.repeat(shuffle_data))
        workers = multiprocessing.Pool(processes=self.processes)
        with open(filename, 'w') as f:

            def write_on_callback(res):
                f.write("\n".join([str(x) for x in res]) + "\n")

            iterable_chunks = chop(iterable, self.processes * chop_size)
            for chunk in iterable_chunks:
                res = workers.map_async(_worker_unit_bis_star, chunk,
                                        callback=write_on_callback)
                res.wait()
        workers.close()
        workers.join()


def worker_unit_bis(dinucleotides_set, (edge_nucleotide, nb_shuffles)):
    """The parallelised piece of code.

    Given a dinucleotide set:
    - generate a given number of sequences
    - calculate the mfe of each
    - return as a string the parameters set and the list of mfe

    """

    res = ""
    shuffler = Shuffle()
    shuffler.init_with_data(dinucleotides_set,
                            edge_nucleotide, edge_nucleotide)
    res += str(sorted(dinucleotides_set.items())) + '\n'
    try:
        shuffles = shuffler.getNshuffles(nb_shuffles)
        mfe_list = RNAfold.get_several_mfe_and_structure(shuffles)
        res += str(sorted([x[0] for x in mfe_list]))
        return res
    except Exception, e:
        raise TaskException(e)


def _worker_unit_bis_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return worker_unit_bis(*a_b)


def chop(iterable, size):
    """Chop an iterable in pieces of a given size.

    Return an iterator over the pieces of the iterable.

    """
    i = iter(iterable)
    piece = list(itertools.islice(i, size))
    while piece:
        yield piece
        piece = list(itertools.islice(i, size))


def read_result_from_file(file_object):
    """Read a result from a given file.

    Return a set line by line.

    """
    while True:
        line1 = file_object.readline()
        line2 = file_object.readline()
        if not line2:
            break  # EOF
        din_set = _parse_dataset_line(line1)
        shuffles = eval(line2)
        yield (din_set, shuffles)


def main():
    """Main function."""

    from argparse import ArgumentParser
    parser = ArgumentParser(description="Generates a Parameters set")
    parser.add_argument('--gen-step', dest="generator_step",
                        type=int, required=True,
                        help='The step of the generator')
    parser.add_argument('--length', dest="length", type=int, required=True,
                        help='The length of the desired sequence')
    parser.add_argument('--shuffles', dest="nb_shuffles",
                        type=int, required=False, default=1000,
                        help='The number of shuffles to perform')
    parser.add_argument('--processes', dest="processes",
                        type=int, required=False, default=1,
                        help='The number of processes to use')
    args = parser.parse_args()
    task = Task(nb_shuffles=args.nb_shuffles, processes=args.processes,
                generator_step=args.generator_step, sequence_length=args.length)
    print task.info()
    task.execute_threaded_bis()


if __name__ == '__main__':
    main()
