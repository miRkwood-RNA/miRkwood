#!/usr/bin/python
# -*- coding: latin-1 -*-

"""Implementation of a Randfold."""

import Shuffle
import RNAfold


def read_simple_fasta(file_handler):
    """A simple FASTA reader, expecting only one sequence."""
    file_handler.next()
    return "".join(file_handler).replace('\n', '').replace('\r', '')


def randfold(sequence_text, iterations):
    """Randfold-like method.

    Return the structure and the p-value for the given sequence.

    """

    (sequence_mfe, sequence_structure) = RNAfold.get_mfe_and_structure(sequence_text)

    shuffled_generator = Shuffle.shuffleNtimes(iterations, sequence_text)
    rnafold_generator = RNAfold.get_several_mfe_and_structure(shuffled_generator)

    below_threshold = len([x[0]
                           for x in rnafold_generator
                           if x[0] < sequence_mfe]
                          )
    pvalue = below_threshold / float(iterations)
    return (sequence_structure, sequence_mfe, pvalue)


def fast_randfold(sequence_text, iterations, p_value=0.1):
    """Randfold-like method."""

    (sequence_mfe, sequence_structure) = RNAfold.get_mfe_and_structure(sequence_text)

    shuffled_generator = Shuffle.shuffleNtimes(iterations, sequence_text)
    rnafold_generator = RNAfold.get_several_mfe_and_structure(shuffled_generator)
    limit = p_value * iterations
    mfe_under_sequence_mfe = []
    while(len(mfe_under_sequence_mfe) <= limit):
        try:
            (mfe, _) = next(rnafold_generator)
            if (mfe < sequence_mfe):
                mfe_under_sequence_mfe.append(mfe)
        except StopIteration:
            break

    p_value = 1
    if len(mfe_under_sequence_mfe) <= limit:
        p_value = len(mfe_under_sequence_mfe) / float(iterations)
    return (sequence_structure, sequence_mfe, p_value)


def main():
    """Main method."""
    from argparse import ArgumentParser
    description = ('Randfold computes the probability that, for a given RNA'
                   'sequence, the Minimum Free Energy (MFE) of the secondary'
                   'structure is different from a distribution of MFE computed'
                   'with random sequences')
    parser = ArgumentParser(description=description)
    source_group = parser.add_argument_group(title='Source')
    group = source_group.add_mutually_exclusive_group(required=True)
    group.add_argument('sequence_file', type=file, nargs='?',
                       help='The file with the sequence to test')
    group.add_argument('--seq', dest='sequence_text', type=str,
                       help='The sequence to test')
    group.add_argument('--fasta', dest='fasta_file', type=file,
                       help='The FASTA file sequence to test')
    parser.add_argument('--iterations', dest="iterations",
                        type=int, required=False, default=100,
                        help='Number of iteration for Randfold (default:100)')
    fast_group = parser.add_argument_group(title='Fast mode')
    fast_group.add_argument('--fast', dest="fast_mode", action='store_true',
                            help='Whether to use the fast Randfold method')
    fast_group.add_argument('--pvalue', dest="pvalue",
                            type=float, required=False, default=0.1,
                            help='The target p-value, in fast mode (default:0.1)')
    args = parser.parse_args()
    if args.sequence_file:
        sequence_text = args.sequence_file.read().replace('\n', '').replace('\r', '')
    elif args.fasta_file:
        sequence_text = read_simple_fasta(args.fasta_file)
    else:
        sequence_text = args.sequence_text
    if args.fast_mode:
        result = fast_randfold(sequence_text, args.iterations, args.pvalue)
        print '\t'.join([str(x) for x in result])
    else:
        result = randfold(sequence_text, args.iterations)
        print '\t'.join([str(x) for x in result])


if __name__ == '__main__':
    main()
