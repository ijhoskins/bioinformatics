# !/usr/bin/env/python
"""Runs filtering of error-free reads."""

import argparse
import logging
import os
import pysam
import sys

FILT_SUFFIX = "filt.bam"
EDIT_DIST = "NM"
MD_TAG = "MD"

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = None
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"


__logger = logging.getLogger(__name__)
__logger.setLevel(logging.DEBUG)
__fhandler = logging.FileHandler("stderr.log")
__fhandler.setLevel(logging.DEBUG)
__formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
__fhandler.setFormatter(__formatter)
__logger.addHandler(__fhandler)


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-a", "--alignments", required=True, type=str, help='Input BAM file.')

    parser.add_argument("-o", "--outdir", type=str, help='Output directory.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


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

    ext_res = "{}.{}".format(split[0], ext)
    return ext_res


def filter_alignments(am, outdir):
    """Filters alignnments that are error-free.

    :param str am: input SAM/BAM file
    :param str outdir: output directory name
    """

    am_filt = replace_extension(am, FILT_SUFFIX)
    output_bam_name = os.path.join(outdir, am_filt)

    with open(output_bam_name, mode="wb") as output_bam, \
            pysam.AlignmentFile(am, "rb") as input_af, \
            pysam.AlignmentFile(output_bam, "wb", header=input_af.header) as output_af:

        for read_aln in input_af.fetch(until_eof=True):

            # NM denotes the edit operations (SNPs, InDels) relative to the reference
            # This value is contained within the optional alignment "tags"
            # if there are no edits, the read is error-free
            if int(read_aln.get_tag("NM")) == 0:

                output_af.write(read_aln)


def workflow(alignments, outdir="."):
    """Filters the reads without error in the alignments.

    :param str alignments: input BAM file
    :param str outdir: Optional output dir for the results
    """

    filter_alignments(am=alignments, outdir=outdir)


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(alignments=parsed_args["alignments"], outdir=parsed_args["outdir"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])

