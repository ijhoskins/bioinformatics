#!/usr/bin/env python3
"""Trims reads from a BAM based on flanking sequences."""

import argparse
import logging
import os
import pysam
import regex
import sys

FASTQ_QNAME_CHAR = "@"
FILE_DELIM = "\t"
FILE_NEWLINE = "\n"
LOG_FORMATTER = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

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

MM_ALLOWANCE = 1
SAM_QUAL_FIELD = 10


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-b", "--bam", type=str, required=True, help='BAM file to trim.')

    parser.add_argument("-f", "--flank_sequences", type=str, required=True,
                        help='Comma-separated sequences flanking the sequence of interest. Will retain flanking sequences.')

    parser.add_argument("-e", "--mm_allowance", type=int, default=MM_ALLOWANCE,
                        help='Mismatch allowance for matching the flanking sequences. Default %i.' % MM_ALLOWANCE)

    parser.add_argument("-o", "--output_dir", type=str, default=".",
                        help='Optional output directory. Default current working directory.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(bam, flank_sequences, mm_allowance=MM_ALLOWANCE, output_dir="."):
    """Runs the BAM trimming workflow.

    :param str bam: input BAM
    :param str flank_sequences: comma-separated flanking sequences
    :param int mm_allowance: mismatch allowance for matching the flanking sequences, default 3
    :param str output_dir: optional output directory
    """

    flank_sequences_split = flank_sequences.split(",")

    # Use regex to enable fuzzy matching
    flank_left_re = regex.compile("(%s){s<=%i}" % (flank_sequences_split[0].upper(), mm_allowance))
    flank_right_re = regex.compile("(%s){s<=%i}" % (flank_sequences_split[1].upper(), mm_allowance))

    output_fn = os.path.join(output_dir, replace_extension(os.path.basename(bam), "trim.fq"))

    with pysam.AlignmentFile(bam, mode="rb", check_sq=False) as input_af, open(output_fn, "w") as output_fh:

        filtered_seqs = 0
        for i, align_seg in enumerate(input_af.fetch(until_eof=True)):

            # Need to search each sequence for the flanking nucleotides
            flank_left = flank_left_re.search(align_seg.query_sequence)
            flank_right = flank_right_re.search(align_seg.query_sequence)

            if not (flank_left and flank_right):
                filtered_seqs += 1
                continue

            # If we match both the left and right flank we can extracte the sequence
            flank_left_idx = flank_left.span()[0]
            flank_right_idx = flank_right.span()[1]

            trim_seq = align_seg.query_sequence[flank_left_idx:flank_right_idx]
            trim_quals_ascii = align_seg.to_string().split(FILE_DELIM)[SAM_QUAL_FIELD][flank_left_idx:flank_right_idx]
            fastq_entry = FILE_NEWLINE.join((FASTQ_QNAME_CHAR + str(i), trim_seq, "+", trim_quals_ascii))

            output_fh.write(fastq_entry + FILE_NEWLINE)

        logger.warning(
            "Filtered out %i reads that did not have matches to both flanking sequences." % filtered_seqs)


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

    workflow(bam=parsed_args["bam"], flank_sequences=parsed_args["flank_sequences"],
             mm_allowance=parsed_args["mm_allowance"], output_dir=parsed_args["output_dir"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
