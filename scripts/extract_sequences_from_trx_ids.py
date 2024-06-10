#!/usr/bin/env python3
"""Extracts sequences from a FASTA file based on Ensembl transcript IDs."""

import argparse
import os
import pysam
import sys

HEADER_DELIM = "|"
TRX_IDX = 0
TRX_EXT = "."
FOREGROUND_EXT = ".extracted.fa"
DEFAULT_OUTDIR = "."
FILE_NEWLINE = "\n"


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    parser.add_argument("-i", "--trx_ids_file", type=str, required=True,
                        help='Ensembl transcript IDs to select, with no minor version extension. One gene per line.')

    parser.add_argument("-f", "--fasta", type=str, required=True,
                        help='Curated APPRIS transcriptome FASTA with pipe-delimited FASTQ header with region information.')

    parser.add_argument("-r", "--make_rna", action="store_true", required=False,
                        help='Flag to convert DNA sequences to RNA sequences.')

    parser.add_argument("-m", "--minimal_name", action="store_true", required=False,
                        help='Flag to output FASTA headers/names using a minimal Ensembl transcript ID.')

    parser.add_argument("-o", "--output_dir", type=str, required=False, default=DEFAULT_OUTDIR,
                        help='Optional output directory. Default current working directory.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    outdir = parsed_args["output_dir"]
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    trx_file = parsed_args["trx_ids_file"]

    foreground_outfile = os.path.join(outdir, os.path.basename(os.path.splitext(trx_file)[0]) + FOREGROUND_EXT)

    with open(trx_file, "r") as trx_fh:
        transcripts = {e.strip() for e in trx_fh}

    fasta_file = parsed_args["fasta"]

    with pysam.FastxFile(fasta_file, "r") as in_fh, open(foreground_outfile, "w") as foreground_fh:

        for rec in in_fh:

            rec_split = rec.name.split(HEADER_DELIM)
            rec_trx = rec_split[TRX_IDX].split(TRX_EXT)[0]

            if parsed_args["make_rna"]:
                rna_seq = "".join(e if e != "T" else "U" for e in rec.sequence)
                rec.sequence = rna_seq

            if rec_trx in transcripts:

                if parsed_args["minimal_name"]:
                    rec.name = rec_trx

                foreground_fh.write(str(rec) + FILE_NEWLINE)


if __name__ == "__main__":
    main()
