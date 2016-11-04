#!/usr/bin/python

import re


def list_variables(filename):
    """Read the output of CLASS and return a list of column numbers and their
    associated variables."""
    with open(filename) as f:
        prev_line = ""
        for line in f:
            if ('#' in line) is False:
                break
            prev_line = line

    h_list = re.split('\d+:', prev_line)[1:]  # Remove '#'
    h_list = map(str.strip, h_list)

    return h_list


def __print_beautifully(h_list):

    h_list = zip(range(len(h_list)), h_list)

    print "('Column', 'Variable')\n"

    for entry in h_list:
        print entry

    print "\nRecall that (.) = 8piG/3"

    return 1


def check_two_header_lists(filename1, filename2):

    if list_variables(filename1) == list_variables(filename2):
        return True
    else:
        return False


def main():

    parser = argparse.ArgumentParser(
        description="Read the output of CLASS and return a list of column\
        numbers and their associated variables.")

    parser.add_argument("filepath", help="Path to file.")

    args = parser.parse_args()

    __print_beautifully(list_variables(args.filepath))

    return 1


if __name__ == "__main__":
    import argparse
    main()
