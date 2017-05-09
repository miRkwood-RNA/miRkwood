#!/usr/local/bin/python2.7
# encoding: utf-8

"""Compute fancy graphs on point-by-point energy distributions."""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from shuffles.Task import Task, read_result_from_file
from shuffles.ParameterGenerator import read_dataset_from_file


def make_shuffles(set_file):
    """Compute the shuffles based on the given parameters file."""
    parameters_set = read_dataset_from_file(set_file)
    task = Task(nb_shuffles=5000, processes=4)
    outfile = set_file.name[:-4] + ".shuf"
    task.process_parameters_set_iterable(parameters_set, outfile)


def make_stats(results_file, count_value='AA', front_text=""):
    """Compute stats and graphs based on the given result file."""
    outfile = results_file.name[:-5] + ".pdf"
    pp = PdfPages(outfile)
    fig1 = plt.figure()
    fig1.suptitle(front_text, fontsize=14, verticalalignment='center')
    pp.savefig(fig1)

    medians = dict()
    first_deciles = dict()
    for (din_set, shuffles) in read_result_from_file(results_file):
        count = din_set[count_value]
        length = sum(din_set.values())
        fig = plt.figure()
        plt.hist(shuffles, histtype='bar', alpha=0.5,
                 bins=numpy.arange(int(min(shuffles)), int(max(shuffles)) + 1, 1))
        plt.title("Energy distribution (#%s = %s)" % (count_value, count))
        plt.xlabel("Energy")
        plt.ylabel("Counts")
        pp.savefig(fig)
        medians[length] = numpy.median(shuffles)
        first_deciles[length] = numpy.percentile(shuffles, 10)

    mean_fig = plt.figure()
    plt.plot(medians.keys(), medians.values(), 'go',
             label="Median")
    plt.plot(first_deciles.keys(), first_deciles.values(), 'ro',
             label="First decile")
    plt.title("Median and first decile evolution")
    plt.xlabel("Sequence length")
    plt.ylabel("Energy")
    plt.legend()
    pp.savefig(mean_fig)
    pp.close()


def main():
    """Main method."""
    parser = argparse.ArgumentParser(description="Computes from a parameters set file")
    parser.add_argument('--results', dest="results", type=argparse.FileType('r'),
                        help='A results file')
    parser.add_argument('--set', dest="set", type=argparse.FileType('r'),
                        help='A parameters set file')
    args = parser.parse_args()

    if args.set:
        make_shuffles(args.set)
    elif args.results:
        text1 = "Uniform frequency : all dinucleotides have the same count x."
        text2 = "Sequences enriched in AT"
        text3 = "Sequences enriched in CG"
        make_stats(args.results, count_value='AU', front_text=text3)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
