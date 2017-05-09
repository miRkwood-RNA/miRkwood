#!/usr/local/bin/python2.7
# encoding: utf-8

"""Evaluate prediction results quality."""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from argparse import ArgumentParser
from shuffles.Task import read_result_from_file
import csv
import itertools
from collections import Counter
import os.path
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import numpy

def read_prediction_results(csv_file_object):
    """Parser of the R output CSV file."""
    csv_file = csv.DictReader(csv_file_object)
    for line in csv_file:
        yield line


def compare_results(csv_file_object, results_file_object):
    """Compare the results line by line."""
    results_gen = read_result_from_file(results_file_object)
    prediction_gen = read_prediction_results(csv_file_object)
    master_gen = itertools.izip(results_gen, prediction_gen)
    for (counter, distribution), prediction_dict in master_gen:
        # Checking if we are considering the same set
        for key, value in counter.items():
            assert float(prediction_dict[key]) == float(value)
        mfe_frontier_predicted = float(prediction_dict['MFE_frontier_computed'])
        inf = len([x for x in distribution if x < mfe_frontier_predicted])
        yield inf / float(len(distribution)) * 100

def make_graph(comp):
    """Make a fancy graph."""
    fig1 = plt.figure()
    fig1.suptitle("Results comparison", fontsize=14, verticalalignment='center')
    bins = list(numpy.arange(0, 20, 0.1))
#    bins.extend([3, 4, 5, 10, int(max(comp)) + 1])
#    plt.hist(comp, histtype='stepfilled', alpha=0.5,
#             bins=bins)
#    plt.xticks([0.2, 0.6, 1, 1.4, 1.8, 2, 3, 4, 5, 10, int(max(comp)) + 1],
#               rotation=30, size='small')
#    plt.xlabel('Pseudo p-value found')
#    plt.ylabel('Counts')
#    plt.axvline(1)
#    plt.savefig('Comparison.pdf')
#
#    fig2 = plt.figure()
#    fig2.suptitle("Results comparison", fontsize=14, verticalalignment='center')
#    plt.hist(comp, histtype='stepfilled', alpha=0.5,
#             bins=[0, 1, int(max(comp)) + 1])
#    plt.axvline(1)
#    plt.savefig('Comparison2.pdf')
    

def main():
    """Main method."""
    parser = ArgumentParser(description="Compare results")
    parser.add_argument('--predict', dest="prediction_csv_file", type=file,
                        required=True,
                        help='The CSV file of the predictions')
    parser.add_argument('--shuf', dest="shuf_file", type=file,
                        required=True,
                        help='The .shuf file')
    args = parser.parse_args()
    comp_gen = compare_results(args.prediction_csv_file, args.shuf_file)
#    comp = list(itertools.islice(comp_gen, 0, 1000))
#    comp = list(comp_gen)
#    make_graph(comp)
    classes = Counter()
    filtered_gen = (abs(1 - item) for item in comp_gen)
    for item in filtered_gen:
        if item <= 0.2:
            classes[0.2] += 1
        elif item <= 0.4:
            classes[0.4] += 1
        elif item <= 0.6:
            classes[0.6] += 1
        else:
            classes['rest'] += 1
    print classes


if __name__ == "__main__":
    main()
