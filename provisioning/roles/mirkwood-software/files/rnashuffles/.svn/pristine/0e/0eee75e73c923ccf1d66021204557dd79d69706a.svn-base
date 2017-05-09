#! /usr/bin/env python
# -*- coding: latin-1 -*-
"""Algorithm to generate parameters sets - 16-uplets of dinucleotides.

(December 2012)
"""
__author__ = "Algorithm by Hélène Touzet, implementation by Jean-Frédéric."

from collections import Counter
from Shuffle import Shuffle
import random
import logging

logging.basicConfig(filename='execution.log', level=logging.INFO)


class ParameterGenerator:

    """Generate Eulerian dinucleotides sets (as 16-uplets).

    Generate all possible dinucleotides sets (as 16-uplets) which as a graph
    admit a Eulerian walk

    """

    def __init__(self, length, alpha, beta, step=1):
        """Initialise the generator.

        length - length of the desired sequence
        alpha - the lower bound for each dinucleotide count
        beta - the upper bound for each dinucleotide count

        """
        self.length = length
        self.n = self.length
        self.alpha = alpha
        self.beta = beta
        self.nucleotides = ['A', 'U', 'G', 'C']
        self.dinucleotides = ['AC', 'AG', 'AU',
                              'CA', 'CG', 'CU',
                              'GA', 'GC',
                              'UA', 'UC',
                              'GU', 'UG',
                              'AA', 'CC', 'GG', 'UU']
        self.uplets = 16
        self.values = Counter()
        self.step = step
        self.positive_results = 0
        self.negative_results = 0
        self.invalid_input = 0
        self.odd_value = 0
        self.ignored = Counter()

    def _already_affected(self, order):
        """Return the number of vertices already affected."""
        return sum([self.values[x] for x in self.dinucleotides[:order]])

    def _smart_range_computer(self, n, already_affected, depth):
        lower_bound = n - already_affected - depth * self.beta
        upper_bound = n - already_affected - depth * self.alpha
        result_range = range(max(lower_bound, self.alpha),
                             min(self.beta, upper_bound) + 1,
                             self.step)
        return result_range

    def _smart_range(self, uplets, level, n):
        """Return the iteration range given the iteration level and the n."""
        depth = uplets - level
        order = uplets - depth - 1
        already_affected = self._already_affected(order)
        lower_bound = n - already_affected - depth * self.beta
        upper_bound = n - already_affected - depth * self.alpha
#        print "n %s - level %s * self.alpha %s - self.beta %s - already_affected %s" %\
#               (n, level, self.alpha, self.beta, already_affected)
#        print "%s / %s < #%s < %s / %s" %\
#               (self.alpha, lower_bound,
#                self.dinucleotides[order],
#                upper_bound, self.beta)
        result_range = range(max(lower_bound, self.alpha),
                             min(self.beta, upper_bound) + 1,
                             self.step)
        return result_range

    def _dumb_range(self, uplets=None, level=None, n=None):
        dinucleotide_range = range(self.alpha, self.beta + 1)
        return dinucleotide_range

    def generate_dinucleotides_sets_with_smartness(self, n):
        """Generate all sets of possible dinucleotides parameters for 12 parameters.

        Generate a parameters set such as a dictionary
        {'AA': 5, 'AT': 2, etc.}

        """
        uplets = 12
        n = self.n
        print "== Length %s" % n
        for self.values[self.dinucleotides[0]] in self._smart_range(uplets, 1, n):
            for self.values[self.dinucleotides[1]] in self._smart_range(uplets, 2, n):
                for self.values[self.dinucleotides[2]] in self._smart_range(uplets, 3, n):
                    for self.values[self.dinucleotides[3]] in self._smart_range(uplets, 4, n):
                        for self.values[self.dinucleotides[4]] in self._smart_range(uplets, 5, n):
                            for self.values[self.dinucleotides[5]] in self._smart_range(uplets, 6, n):
                                for self.values[self.dinucleotides[6]] in self._smart_range(uplets, 7, n):
                                    for self.values[self.dinucleotides[7]] in self._smart_range(uplets, 8, n):
                                        try:
                                            self._compute_end_of_set(n)
                                            yield self.values
                                        except ParameterGenerationException, _:
                                            pass

    def _compute_end_of_set(self, n, fail_if_out_of_range=False):
        """Based on a 8-uplets, compute the rest."""
#        print "_compute_end_of_set()"
        uplets = 12

        def _exception_message(dinucl, dinucl_value, dinucl_range):
            return"%s = %s - out of range [%s %s]" % (dinucl,
                                                      dinucl_value,
                                                      dinucl_range[0],
                                                      dinucl_range[-1])

        # Computing UA
        value_UA = (self.values['AC'] + self.values['AG'] + self.values['AU'] -
                    self.values['CA'] - self.values['GA'])

        if fail_if_out_of_range:
            range_UA = self._smart_range(uplets, 9, n)
            if value_UA not in range_UA and fail_if_out_of_range:
                self.ignored['UA'] += 1
                raise ParameterGenerationException(_exception_message('UA', value_UA, range_UA))
        self.values['UA'] = value_UA

        # Computing UC
        value_UC = (self.values['CA'] + self.values['CG'] + self.values['CU'] -
                    self.values['AC'] - self.values['GC'])

        if fail_if_out_of_range:
            range_UC = self._smart_range(uplets, 10, n)
            if value_UC not in range_UC:
                self.ignored['UC'] += 1
                raise ParameterGenerationException(_exception_message('UC', value_UC, range_UC))

        self.values['UC'] = value_UC

        # Computing GU
        term = (self.values['AC'] + self.values['AG'] +
                self.values['CA'] + self.values['CG'] +
                self.values['GA'] + self.values['GC'] +
                2 * self.values['AU'] + 2 * self.values['CU'])
        numerator = n - term

        if numerator % 2 != 0:
            self.odd_value += 1
            raise ParameterGenerationException("Odd value")
        value_GU = numerator / 2

        if fail_if_out_of_range:
            range_GU = self._smart_range(uplets, 11, n)
            if value_GU not in range_GU:
                self.ignored['GU'] += 1
                raise ParameterGenerationException(_exception_message('GU', value_GU, range_GU))

        self.values['GU'] = value_GU

        # Computing UG
        value_UG = (self.values['GU'] + self.values['AU'] + self.values['CU'] -
                    self.values['UA'] - self.values['UC'])

        if fail_if_out_of_range:
            range_UG = self._smart_range(uplets, 12, n)
            if value_UG not in range_UG:
                self.ignored['UG'] += 1
                raise ParameterGenerationException(_exception_message('UG', value_UG, range_UG))
        self.values['UG'] = value_UG

        if all([x >= 0 for x in self.values.values()]):
            self.positive_results += 1
            return self.values
        else:
            raise ParameterGenerationException("Invalid parameters set")

    def generate_dinucleotides_sets_without_loops(self):
        """Generate all dinucleotides sets, without loops.

        Generator function which yields a parameter set as a dictionary:
        {'AA': 5, 'AT': 2, etc.}

        """
        for values in self.generate_dinucleotides_sets_with_smartness(self.length):
            values['AA'] = values['CC'] = values['GG'] = values['UU'] = 0
            yield values

    def generate_dinucleotides_sets_with_loops(self):
        """Generate all dinucleotides sets, with loops.

        Generator function which yields a parameter set as a dictionary:
        {'AA': 5, 'AT': 2, etc.}

        """
        method = ParameterGenerator.generate_dinucleotides_sets_with_smartness
        return self.generate_loops_wrapper(method)

    def generate_dinucleotides_sets_with_dumbness(self):
        values = {}

        def din(rank):
            return sorted(self.dinucleotides)[rank]

        for values[din(0)] in self._dumb_range():
            for values[din(1)] in self._dumb_range():
                for values[din(2)] in self._dumb_range():
                    for values[din(3)] in self._dumb_range():
                        for values[din(4)] in self._dumb_range():
                            for values[din(5)] in self._dumb_range():
                                for values[din(6)] in self._dumb_range():
                                    for values[din(7)] in self._dumb_range():
                                        for values[din(8)] in self._dumb_range():
                                            for values[din(9)] in self._dumb_range():
                                                for values[din(10)] in self._dumb_range():
                                                    for values[din(11)] in self._dumb_range():
                                                        for values[din(12)] in self._dumb_range():
                                                            for values[din(13)] in self._dumb_range():
                                                                for values[din(14)] in self._dumb_range():
                                                                    for values[din(15)] in self._dumb_range():
                                                                        yield values
                print values

    def generate_dinucleotides_sets_with_smarter_dumbness(self):
        print "generate_dinucleotides_sets_with_smarter_dumbness"

        n = self.n

        def din(rank):
            return self.dinucleotides[rank]

        for self.values[din(0)] in self._dumb_range(1):
            for self.values[din(1)] in self._dumb_range(2):
                for self.values[din(2)] in self._dumb_range(3):
                    for self.values[din(3)] in self._dumb_range(4):
                        for self.values[din(4)] in self._dumb_range(5):
                            for self.values[din(5)] in self._dumb_range(6):
                                for self.values[din(6)] in self._dumb_range(7):
                                    for self.values[din(7)] in self._dumb_range(8):
                                        try:
                                            self._compute_end_of_set(n)
                                            yield self.values
                                        except ParameterGenerationException, _:
                                            pass

    def generate_loops_wrapper(self, parameters_set_method):
        print "generate_loops_wrapper(%s)" % parameters_set_method
        non_looping_range = range(self.length - 4 * self.beta,
                                  self.length - 4 * self.alpha + 1)
        non_looping_range = filter(lambda x: x > 0, non_looping_range)
        for non_looping in non_looping_range:
            self.n = non_looping
            for res in parameters_set_method(self):
                shuffler = Shuffle()
                if not shuffler.check_valid_input(res, 'U', 'U'):
                    pass
                else:
                    range_AA = self._dumb_range(self.uplets, 13, self.length)
                    for self.values['AA'] in range_AA:
                        range_CC = self._dumb_range(self.uplets, 14, self.length)
                        for self.values['CC'] in range_CC:
                            range_GG = self._dumb_range(self.uplets, 15, self.length)
                            for self.values['GG'] in range_GG:
                                value_UU = (self.length -
                                            self._already_affected(15))
                                if value_UU not in self._dumb_range():
                                    pass
                                else:
                                    self.values['UU'] = value_UU
#                                    print "Accepted"
#                                    print self.values
                                    yield self.values.copy()
        pass

    def make_with_randomness(self):
        uplets = 12
        n = self.n
        def random_dinucleotide(level):
            dinucleotide_range = self._dumb_range(uplets, level, n)
#            return random.randint(self.alpha, self.beta + 1)
            return random.choice(dinucleotide_range)

        for level, dinucleotide in enumerate(self.dinucleotides[:8]):
            self.values[dinucleotide] = random_dinucleotide(level + 1)

        try:
            self._compute_end_of_set(self.n)
            return self.values
        except Exception, exception:
#            logging.warning(exception)
            raise exception

    def make_full_with_randomness(self):
        values = {}
        uplets = 16
        n = self.length
        level = 0
        def random_dinucleotide(level):
            loops_affected = sum(values.values())
            depth = uplets - level
#            dinucleotide_range = self._smart_range_computer(self.length, loops_affected, depth)
#            dinucleotide_range = self._smart_range(uplets, level, n)
            dinucleotide_range = self._dumb_range()
            return random.choice(dinucleotide_range)

        for level, dinucleotide in enumerate(self.dinucleotides[-4:]):
            values[dinucleotide] = random_dinucleotide(level + 1)
        already_affected = sum([values[x] for x in self.dinucleotides[-4:]])
        self.n = self.length - already_affected
        values.update(self.make_with_randomness().copy())
        return values

    def generate_amount_with_randomness(self, amount=1000000):
        print "generate_amount_with_randomness"
        counter = 0
        try:
            while counter < amount:
                try:
                    item = self.make_full_with_randomness()
                    counter += 1
                    yield item
                except ParameterGenerationException, e:
#                    print "Bad set, ignoring"
                    continue
        except Exception, e:
            print "Problem with master for: " + str(e)
            raise e

    def make_with_uniform_frequency(self, count):
        """All dinucleotides have the same count x.
        #AA=#AC=#AG=#AT=#CA=#CC= ....= x
        """
        for din in self.dinucleotides:
            self.values[din] = count
        return self.values

    def make_with_enriched_AU(self, count):
        """Sequences enriched in AT
        #AA = #AU = #UU = #UA = 2x
        #XX = x
        """
        AU_dinucleotides = ['AA', 'AU', 'UU', 'UA']
        other_dinucleotides = ['AC', 'AG',
                               'CA', 'CG', 'CU',
                               'GA', 'GC',
                               'UC',
                               'GU', 'UG',
                               'CC', 'GG']
        for din in AU_dinucleotides:
            self.values[din] = 2 * count
        for din in other_dinucleotides:
            self.values[din] = count
        return self.values

    def make_with_enriched_CG(self, count):
        """Sequences enriched in CG
        #CC = #CG = #GG = #GC = 2x
        #XX = x
        """
        CG_dinucleotides = ['CC', 'CG', 'GG', 'GC']
        other_dinucleotides = ['AC', 'AG', 'AU',
                               'CA', 'CU',
                               'GA',
                               'UA', 'UC',
                               'GU', 'UG',
                               'AA', 'UU']
        for din in CG_dinucleotides:
            self.values[din] = 2 * count
        for din in other_dinucleotides:
            self.values[din] = count
        return self.values

    def generate_lots(self):
        """A special generation process."""
        def generate1(self):
            for count in range(4, 21):
                yield self.make_with_uniform_frequency(count)

        def generate2(self):
            for count in range(3, 17):
                yield self.make_with_enriched_AU(count)

        def generate3(self):
            for count in range(3, 17):
                yield self.make_with_enriched_CG(count)

#        self.generate_and_dump('point-a-point-1.txt', generate1(self))
#        self.generate_and_dump('point-a-point-2.txt', generate2(self))
        self.generate_and_dump('point-a-point-3.txt', generate3(self))

    def generate_and_dump(self, fileName, parameter_generator):
        """Dump a given set of parameters on disk."""
        print "n = %s" % self.n
        counter = 0
        ignored = 0
        ignoredbis = 0
        try:
            with open(fileName, 'w') as f:
                print "File opened"
                for item in parameter_generator:
                    shuffler = Shuffle()
                    if not shuffler.check_valid_input(item, 'U', 'U'):
                        self.invalid_input += 1
#                    elif sum(item.values()) != self.n:
#                        ignoredbis += 1
                    else:
#                        f.write(str(dict(item)) + '\n')
                        f.write(str(sorted(item.items())) + '\n')
                        counter += 1
            print "Values %s" % counter
            print "Ignored %s" % ignored
            print "Ignored bis %s" % ignoredbis
#            self._post_mortem()
        except Exception, e:
            print "Could not write value"
            print e
        self._post_mortem()

    def _post_mortem(self):
        print "Invalid %s" % self.invalid_input
        print "Positive: %s " % self.positive_results
        print "Negative: %s " % self.negative_results
        print "Odd 'GU': %s " % self.odd_value


class ParameterGenerationException(Exception):
    pass


def read_dataset_from_file_old(file_object):
    """Read a parameters set from a given file.

    Return a set line by line.

    """
    for line in file_object:
        if line == "\n":
            continue
        res_i = Counter()
        splitted = line[1:-2].split(',')
        for item in splitted:
            (key, value) = item.split(':')
            res_i[key.strip()[1:-1]] = int(value.strip())
        yield res_i


def read_dataset_from_file(file_object):
    """Read a parameters set from a given file.

    Return a set line by line.

    """
    for line in file_object:
        if line == "\n":
            continue
        yield _parse_dataset_line(line)


def _parse_dataset_line(line):
    """Parse a parameters dataset line.

    Return a Counter.

    """
    res_i = Counter()
    for (x, y) in eval(line):
        res_i[x] = y
    return res_i


def main():
    """Main method."""
    actions = {'smart': ParameterGenerator.generate_dinucleotides_sets_with_smartness,
               'dumb': ParameterGenerator.generate_dinucleotides_sets_with_dumbness,
               'notsodumb': ParameterGenerator.generate_dinucleotides_sets_with_smarter_dumbness,
               'random': ParameterGenerator.generate_amount_with_randomness,
#               'pointbypoint': ParameterGenerator.generate_lots,
               }

    from argparse import ArgumentParser
    parser = ArgumentParser(description="Generates a Parameters set")
    parser.add_argument('--step', dest="step", type=int, required=False,
                        help='The step of the generator')
    parser.add_argument('--length', dest="length", type=int, required=True,
                        help='The length of the desired sequence')
    parser.add_argument('--method', dest="method", type=str, required=True,
                        choices=actions.keys(),
                        help='The parameter generation method')
    args = parser.parse_args()

    alpha = int(args.length / 32)
    beta = int(args.length / 8) + 1
    generator = ParameterGenerator(args.length, alpha, beta, step=args.step)

    #Running the relevant parameter generation method
    try:
        parameters_set = actions[args.method](generator)
        fileNameRoot = "parameters"
        fileName = "%s-%s-%s.txt" % (fileNameRoot, args.method, args.length)
        generator.generate_and_dump(fileName, parameters_set)
    except KeyboardInterrupt:
        print "==== KeyboardInterrupt ===="
        print generator._post_mortem()

if __name__ == '__main__':
    main()
