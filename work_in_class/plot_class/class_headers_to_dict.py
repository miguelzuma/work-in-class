#!/usr/bin/python

import variable_title_dics as vtd
import list_class_headers as lch
import re
import sys
from collections import defaultdict


def class_headers_list_to_dict(h_list):
    """Take a list of column headers and return a dictionary {variable: column}
    """

    title_col_list = zip(h_list, range(len(h_list)))
    title_col_list.sort()  # Make sure var[N] are in order.

    var_col_dic = defaultdict(list)

    for title, column in title_col_list:
        try:
            var = vtd.title_var_dic[title]
            var_col_dic[var] = column
        except KeyError:
            title_split = filter(None, re.split('\W', title))

            if 'ncdm' in title or len(title_split) == 3:
                var = vtd.title_var_dic[title[:-3]]  # [:-3] removes [N]
                var_col_dic[var].append(column)

            elif len(title_split) == 4:
                var = vtd.title_var_dic[title_split[0] + title_split[2]]
                x = title_split[1]

                try:
                    var_col_dic[var][x].append(column)
                except:
                    var_col_dic[var][x] = [column]

            else:
                print title
                sys.exit("Column title " + title + " not implemented. Please\
                    see variable_title_dics.py and consider adding it")

    return var_col_dic


def class_headers_to_dict(filename):
    """Read the CLASS output and return a dictionary {variable: column}"""

    h_list = lch.list_variables(filename)
    var_col_dic = class_headers_list_to_dict(h_list)

    return var_col_dic


def class_file_title_var_dic(filename):
    """Read the CLASS output and return a dictionary {title: var} with just the
    elements in the present file"""

    h_list = lch.list_variables(filename)
    title_var_dic = {}

    for title in h_list:
        var = vtd.title_var_dic[title]
        title_var_dic[title] = var

    return title_var_dic  # Just the titles present in the read file


def main():

    parser = argparse.ArgumentParser(
        description="Read the output of CLASS and return a dictionary relating\
        a variable name with its column number.")

    parser.add_argument("filepath", help="Path to file.")

    args = parser.parse_args()

    print class_headers_to_dict(args.filepath)

    return 1


if __name__ == "__main__":
    import argparse
    import list_class_headers as lch
    main()
