#!python
# encoding: utf-8

"""Compute shuffles on a given file."""

import argparse


def make_shuffles(set_file, processes=1):
    """Compute the shuffles based on the given parameters file."""
    parameters_set = read_dataset_from_file(set_file)
    task = Task(nb_shuffles=1000, processes=processes)
    outfile = set_file.name[:-4] + ".shuf"
    task.process_parameters_set_iterable(parameters_set, outfile)


def main():
    """Main method."""
    parser = argparse.ArgumentParser(description="Run a full RandomForest thingy")
    parser.add_argument('--threads', dest="threads", type=int,
                        help='Number of threads to use')
    parser.add_argument('--set', dest="set", type=argparse.FileType('r'),
                        help='A parameters set file')
    args = parser.parse_args()

    if args.set:
        make_shuffles(args.set, processes=args.threads)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
