#!/usr/bin/env/python
"""Renames files according to a naming map file."""

import argparse
import os
import sys
import warnings

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

FILE_NEWLINE = "\n"
FIELD_DELIM = "\t"
WHITESPACE = " "


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="")

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-m", "--map_file", type=str, required=True,
                        help='Mapping .txt file containing a header line, input filepath in first column, and '
                             'output filepath in second column.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(map_file):
    """Runs the file renaming workflow.

    :param str map_file: name mapping file
    """

    with open(map_file, "r") as map_fh:

        for i, line in enumerate(map_fh):

            # Skip the header
            if i == 0:
                continue

            # Strip newline, any whitespace, then split on tab
            line_strip = line.strip(FILE_NEWLINE + WHITESPACE)
            line_split = line_strip.split(FIELD_DELIM)

            if len(line_split) != 2:
                raise NotImplementedError("Exactly two fields not detected for line: %s." % line_strip)

            input_filename, output_filename = line_split

            if not os.path.exists(input_filename):
                warnings.warn("File %s does not exist. Skipping" % input_filename)
                continue

            if os.path.exists(output_filename):
                raise NotImplementedError("File %s already exists." % output_filename)

            os.rename(input_filename, output_filename)


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])
    workflow(map_file=parsed_args["map_file"])


if __name__ == "__main__":
    main()
