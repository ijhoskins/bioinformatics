#!/usr/bin/env/python
"""Filters a GFF/GTF by transcript IDs."""

import argparse
import logging
import os
import pybedtools
import sys

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

__logger = logging.getLogger(__name__)

# Consider putting the following in a logging config file
__logger.setLevel(logging.DEBUG)
__fhandler = logging.FileHandler("stderr.log")
__fhandler.setLevel(logging.DEBUG)
__formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
__fhandler.setFormatter(__formatter)
__logger.addHandler(__fhandler)

DEFAULT_EXT = "filt.gtf"
DEFAULT_OUTDIR = "."
FILE_NEWLINE = "\n"
FILE_DELIM = "\t"
GFF_ATTR_TRANSCRIPT_ID = "transcript_id"


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    parser.add_argument("-g", "--gff", type=str, required=True, help="Gencode GTF to filter.")

    parser.add_argument("-i", "--ids", type=str, required=True, help="Text file of transcript IDs, one per line.")

    parser.add_argument("-e", "--ext", type=str, default=DEFAULT_EXT, help="Extension for output GFF (e.g. set_A.gff).")

    parser.add_argument("-o", "--outdir", type=str, default=DEFAULT_OUTDIR,
                        help='Output directory. Default current working directory.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


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


def extract_gff_records(gff, ids, outdir=".", ext=DEFAULT_EXT):
    """Extracts GFF records for specific transcripts.

    :param str gff: Ensembl GFF/GTF filename
    :param str ids: Ensembl transcript IDs, without minor version number, one per line
    :param str ext: optional output extension for the filtered GFF
    :param str outdir: optional output directory.
    :return str: filtered GFF filepath
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    outfile = os.path.join(outdir, replace_extension(os.path.basename(gff), ext))

    with open(ids, "r") as ids_fh:
        trx_ids = {e.strip(FILE_NEWLINE) for e in ids_fh}

    gff_bedtool = pybedtools.BedTool(gff)

    with open(outfile, "w") as out_gff:
        for interval in gff_bedtool:

            if GFF_ATTR_TRANSCRIPT_ID not in interval.attrs:
                continue

            trx_id = interval.attrs[GFF_ATTR_TRANSCRIPT_ID].split(".")[0]
            
            if trx_id in trx_ids:
                out_gff.write(FILE_DELIM.join(interval.fields) + FILE_NEWLINE)

    return outfile


def workflow(gff, ids, ext=DEFAULT_EXT, outdir="."):
    """Filter a GFF by transcript IDs.

    :param str gff: Ensembl GFF/GTF filename
    :param str ids: Ensembl transcript IDs, without minor version number, one per line
    :param str ext: optional output extension for the GFF (e.g. set_A.gff)
    :param str outdir: optional output dir for the results
    :return str: filtered GFF filepath
    """

    # Uses default file extension from the GFF, written to the outdir by side-effect
    filt_gff = extract_gff_records(gff=gff, ids=ids, ext=ext, outdir=outdir)
    return filt_gff


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(gff=parsed_args["gff"], ids=parsed_args["ids"], ext=parsed_args["ext"], outdir=parsed_args["outdir"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
