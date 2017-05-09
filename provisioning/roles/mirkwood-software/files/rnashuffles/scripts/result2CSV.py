#!/usr/local/bin/python
# encoding: utf-8

"""Convert a result set into CSV file."""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from argparse import ArgumentParser
from shuffles.Task import read_result_from_file
from os.path import splitext
import csv
import numpy


def result_to_csv(result, csv_file_name, pvalue=10):
    """Write a result as CSV."""
    with open(csv_file_name, 'w') as csv_file_object:
        writer = csv.writer(csv_file_object, delimiter=';')
        dinucleotides = ['AC', 'AG', 'AU',
                         'CA', 'CG', 'CU',
                         'GA', 'GC',
                         'UA', 'UC',
                         'GU', 'UG',
                         'AA', 'CC', 'GG', 'UU']
        header = sorted(dinucleotides)
        header.append('MFE_frontier')
        writer.writerow(header)
        for (parameters, distribution) in result:
            values = [x[1] for x in sorted(parameters.items())]
            mfe_frontier = float(numpy.percentile(distribution, pvalue))
            values.append(mfe_frontier)
            writer.writerow(values)


def main():
    """Main method."""
    parser = ArgumentParser(description="Convert a result set into CSV")
    parser.add_argument("-p", '--pvalue', dest='pvalue', type=int,
                        help='The pvalue to use (in percentage)')
    parser.add_argument("result_file", type=file,
                        help='The results file')
    args = parser.parse_args()
    csv_file_name = splitext(args.result_file.name)[0] + '.csv'

    result = read_result_from_file(args.result_file)
    result_to_csv(result, csv_file_name, pvalue=args.pvalue)

if __name__ == "__main__":
    main()
