'''
Checking Shuffles project in various manners.
'''
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from shuffles.ParameterGenerator import read_dataset_from_file
from collections import Counter


def second_check(input_file):
    """Implement the second check method.

    Return the sums of each dinucleotides counts for all parameters sets
    (for the given length/step), as a Counter.

    """

    counter = Counter()
    dataset = read_dataset_from_file(input_file)
    [counter.update(x) for x in dataset]
    input_file.close()
    return counter


def main(args):
    if args.second:
        print second_check(args.second)
    if args.third:
        print sorted(dict(eval(args.third)).items())

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Checks the values of a parameters set")
    parser.add_argument("--second", type=file,
                        help="Second check method")
    parser.add_argument("--third",
                        help="Third check method")
    args = parser.parse_args()
    main(args)
