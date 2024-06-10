#!/usr/bin/env python3
"""Filters CCLE matrix based on cell line depMap profile IDs and genes of interest."""

# See https://depmap.org/portal/download/all/

import argparse
import logging
import os
import sys

FILE_DELIM = "\t"
FILE_NEWLINE = "\n"
LOG_FORMATTER = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
DEPMAP_ID = "depMap_ID"

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"


def add_extension(filename, ext):
    """Adds an extension.

    :param str filename: file path
    :param str ext: extension to add
    """

    ext_res = ".".join((filename, ext,))
    return ext_res


def replace_extension(filename, ext, ignore_exts=(".gz", ".bz", ".bz2",)):
    """Replaces extension of a filename.

    :param str filename: file path
    :param str ext: extension to add
    :param tuple ignore_exts: extensions to strip before replacing the extension
    :return str: filepath with new extension
    """

    split = os.path.splitext(filename)

    if split[1] in set(ignore_exts):
        split = os.path.splitext(split[0])

    ext_res = add_extension(split[0], ext)
    return ext_res


LOGFILE = replace_extension(os.path.basename(__file__), "log")
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(LOG_FORMATTER)
logger.addHandler(console_handler)


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-m", "--CCLE_matrix", type=str, required=True,
                        help='CCLE transcript expression matrix with depMap IDs in first column and gene names in first row.')

    parser.add_argument("-g", "--genes", type=str, required=True, help='Genes, one per line.')

    parser.add_argument("-d", "--depMap_IDs", type=str, required=True,
                        help='depMap IDs for cell lines of interest (e.g. ACH-001188)')

    parser.add_argument("-o", "--output_dir", type=str, default=".",
                        help='Optional output directory. Default current working directory.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(ccle_matrix, genes, cell_line_ids, output_dir="."):
    """Runs the CCLE matrix filtering workflow.

    :param str ccle_matrix: CCLE matrix
    :param str genes: genes to select
    :param str cell_line_ids: cell lines to select (depMap IDs)
    :param str output_dir: optional output directory
    """

    output_file = os.path.join(output_dir, replace_extension(os.path.basename(ccle_matrix), "filt.txt"))

    with open(genes, "r") as gene_fh, open(cell_line_ids, "r") as cell_line_fh, \
            open(ccle_matrix, "r") as ccle_fh, open(output_file, "w") as output_fh:

        gene_set = set(gene_fh.read().split(FILE_NEWLINE))
        cell_line_set = set(cell_line_fh.read().split(FILE_NEWLINE))

        # Remove any empty genes do to empty lines
        gene_set -= {""}
        cell_line_set -= {""}

        counter = 0
        for i, line in enumerate(ccle_fh):

            line_split = line.split(",")

            # Determine which columns to extract
            if i == 0:
                col_indices = set([i for i, name in enumerate(line_split) if name.split(" ")[0] in gene_set])
                line_selected = [field for i, field in enumerate(line_split) if i in col_indices]
                output_fh.write(FILE_DELIM.join([DEPMAP_ID] + line_selected) + FILE_NEWLINE)
                continue

            # Extract the fields for genes of interest
            if line_split[0] in cell_line_set:
                line_selected = [field for i, field in enumerate(line_split) if i in col_indices]
                output_fh.write(FILE_DELIM.join([line_split[0]] + line_selected) + FILE_NEWLINE)
                counter += 1

        logger.info("Filtered %i cell-line rows." % counter)


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    outdir = parsed_args["output_dir"]
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    log_handler = logging.FileHandler(os.path.join(outdir, LOGFILE))
    log_handler.setFormatter(LOG_FORMATTER)
    logger.addHandler(log_handler)

    logger.info("Started %s" % sys.argv[0])

    workflow(ccle_matrix=parsed_args["CCLE_matrix"], genes=parsed_args["genes"],
             cell_line_ids=parsed_args["depMap_IDs"], output_dir=parsed_args["output_dir"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
