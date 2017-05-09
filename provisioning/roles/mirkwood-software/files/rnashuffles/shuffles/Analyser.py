#!/usr/local/bin/python2.7
# encoding: utf-8
"""Sequence analyser
"""

import sys
import os
import logging
from collections import Counter
from Bio import SeqRecord
from Bio import SeqIO
from Shuffle import DataStructuresConversion
from Shuffle import Shuffle
import matplotlib.pyplot as plt

dinucleotides = ['AC', 'AG', 'AU',
                 'CA', 'CG', 'CU',
                 'GA', 'GC',
                 'UA', 'UC',
                 'GU', 'UG',
                 'AA', 'CC', 'GG', 'UU']


def sequence_analyser(sequence):
    """Return the dinucleotides counts of a given sequence."""
    shuffle = Shuffle()
    shuffle.init_with_sequence(sequence)
    dinucl_count = DataStructuresConversion.intricate2plain(shuffle._dinucl_count)
    return dinucl_count


def fasta_file_analyser(fasta_file_handler):
    """Return the dinucleotides counts of a given sequence."""

    normalised_counter = {}
    for dinucl in dinucleotides:
        normalised_counter[dinucl] = []
    for sequence_record in SeqIO.parse(fasta_file_handler, "fasta"):
        sequence = str(sequence_record.seq)
        try:
            dinucl_count = sequence_analyser(sequence)
        except KeyError, error:
            logging.warning("Ignore sequence, contains unknown character %s",
                            error)
            continue
        sequence_length = len(sequence)
        for key, value in dinucl_count.items():
#            print "%s / %s" % (value, sequence_length - 1)
            normalised = (100 * value) / (sequence_length - 1)
            normalised_counter[key].append(normalised)

#    values = normalised_counter['AA']
#    print sorted(values)
#    print len(values)
#    print sum(values)
#    for key, values in normalised_counter.items():
#        print values
    plt.hist(normalised_counter.values(),
             label=normalised_counter.keys(),
             bins=(0, 100/32, 100/8, 30),
             histtype='bar', normed=True, alpha=0.5)

    plt.title("Dinucleotides counts in MiRbase hairpins")
    plt.xlabel("Dinucleotide count")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()
#    for key, values in normalised_counter.items():
#        print values
#        print min(values)
#        print max(values)


def main(args=None):
    """Main method."""
    for fasta_file_handler in args.fasta_files:
        fasta_file_analyser(fasta_file_handler)


if __name__ == "__main__":
    from argparse import ArgumentParser, FileType
    parser = ArgumentParser(description="Analyses a FASTA file for dinucleotides counts")
    parser.add_argument(dest="fasta_files", metavar="files", nargs='+',
                        help="FASTA files to analyse",
                        type=FileType('r'))
    args = parser.parse_args()
    main(args)
