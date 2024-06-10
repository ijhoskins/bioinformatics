#!/usr/bin/env python3
"""Extracts sequences from a FASTA file based on gene names."""

import argparse
import os
import pysam
import sys

HEADER_DELIM = "|"
GENE_IDX = 5
FOREGROUND_EXT = ".extracted.fa"
DEFAULT_OUTDIR = "."
FILE_NEWLINE = "\n"


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    parser.add_argument("-i", "--gene_ids_file", type=str, required=True,
                        help='Gene IDs to select for enumeration. One gene per line.')

    parser.add_argument("-f", "--fasta", type=str, required=True,
                        help='Curated APPRIS transcriptome FASTA with pipe-delimited FASTQ header with region information.')

    parser.add_argument("-r", "--make_rna", action="store_true", required=False,
                        help='Flag to convert DNA sequences to RNA sequences.')

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

    gene_file = parsed_args["gene_ids_file"]

    foreground_outfile = os.path.join(outdir, os.path.basename(os.path.splitext(gene_file)[0]) + FOREGROUND_EXT)

    with open(gene_file, "r") as gene_fh:
        genes = {e.strip() for e in gene_fh}

    fasta_file = parsed_args["fasta"]

    with pysam.FastxFile(fasta_file, "r") as in_fh, open(foreground_outfile, "w") as foreground_fh:

        for rec in in_fh:

            rec_split = rec.name.split(HEADER_DELIM)
            rec_gene = rec_split[GENE_IDX]

            if parsed_args["make_rna"]:
                rna_seq = "".join(e if e != "T" else "U" for e in rec.sequence)
                rec.sequence = rna_seq

            if rec_gene in genes:
                foreground_fh.write(str(rec) + FILE_NEWLINE)


if __name__ == "__main__":
    main()
