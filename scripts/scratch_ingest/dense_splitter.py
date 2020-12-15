"""Dense splitter

DESCRIPTION
Takes a dense matrix file and subset it by cells (columnwise) into X smaller dense files.
Subsetting is purely columnwise and does not take into account the sparsity of gene expression values.

By default, splits on tabs; use "--delimiter" to specify comma for csv files.

INPUT ASSUMPTIONS
Assumes input filename format of <filename>.<filetype>.<if compressed, .gz>
Checks for .gz suffix to handle compressed input files.

Outputs sub_<file number>.<filetype> files wherever the specified input file lives.

MANUAL POSTPROCESSING ASSUMPTIONS
rename subfiles to reflect input filename
compress subfiles before uploading to portal

USAGE
python dense_splitter.py <file> <number of files>

EXAMPLE - subset RNA_expression.csv.gz into four files, splitting on commas
    outputs sub_1.csv, sub_2.csv, sub_3.csv, sub_4.csv, in the troubleshooting/Z54407_SCP1064_cfrangie directory
python dense_splitter.py --delimiter "," troubleshooting/Z54407_SCP1064_cfrangie/RNA_expression.csv.gz 4
"""

import argparse
import gzip
import sys
from glob import glob


def create_parser():
    """Parse command line values
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--delimiter", help="file delimiter", default="\t")
    parser.add_argument("input", help="input TSV file")
    parser.add_argument("num_files", help="number of resulting dense files")
    return parser


def check_if_old_output():
    """Exit if old output files found
    """
    output_filename_format = "sub_[0-9]*"

    old_output = False

    if glob(output_filename_format):
        print(f"dense split files already exist, please delete files and try again")
        print(f'found {glob("sub_[0-9]*")}')
        old_output = True
    if old_output:
        exit(1)


if __name__ == "__main__":
    args = create_parser().parse_args()
    arguments = vars(args)

    args.num_files = int(args.num_files)
    if args.num_files <= 0:
        print("num_files must be positive.")
        sys.exit(1)

    check_if_old_output()

    # matrix sub-files
    filename_parts = args.input.split(".")
    if filename_parts[-1] == "gz":
        file_is_gzipped = True
    else:
        file_is_gzipped = False
    if file_is_gzipped:
        file_suffix = filename_parts[-2]
        file_prefix = ".".join(filename_parts[:-2])
    else:
        file_suffix = filename_parts[-1]
        file_prefix = ".".join(filename_parts[:-1])

    if file_is_gzipped:
        print("opening gzipped file")
        f = gzip.open(args.input, "rt")
    else:
        f = open(args.input)
    line = f.readline()
    # If line starts with "GENE", split on space is correct but
    # if R-style file (ie. starts with empty field) "fix" empty field
    if line[0].upper() != "G":
        line = "GENE" + line
    num_fields = len(line.split(args.delimiter)) - 1
    bin = num_fields // args.num_files
    last_bin = bin + (num_fields % args.num_files)
    print(f"number of fields: {num_fields}, bin size: {bin}, last_bin: {last_bin}")
    f.close()

    if file_is_gzipped:
        f = gzip.open(args.input, "rt")
    else:
        f = open(args.input)

    start = bin + 1
    end = start + bin
    for line in f:
        values = line.split(args.delimiter)
        for sub_index in range(0, args.num_files):
            sub_filename = f"sub_{sub_index + 1}.{file_suffix}"
            start = (bin * sub_index) + 1
            end = start + bin
            if sub_index == (args.num_files - 1):
                is_last_bin = True
                end = start + last_bin
            else:
                is_last_bin = False
            with open(sub_filename, "a+") as sub:
                extracted_values = values[0:1] + values[start:end]
                if is_last_bin:
                    subsetted_values = args.delimiter.join(extracted_values)
                else:
                    subsetted_values = args.delimiter.join(extracted_values) + "\n"
                sub.write(subsetted_values)
    f.close()
