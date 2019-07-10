'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Alistair Legione, 11 Jun 2019
License     : MIT
Maintainer  : legionea@unimelb.edu.au
Portability : POSIX

This program is a basic python coversion of Mick Watson's Ideel.
It reads one or more input FASTA files and for each file it will use
prodigal for rapid annotation, then run diamond blast, then compare the
query length to hit length.

It was built with the help of 'Bionitio'
'''

import argparse
import os
import sys
import subprocess
import logging
import pkg_resources
import pandas
import altair
import seaborn
import selenium
import datetime
#from Bio import SeqIO


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_FASTA_FILE_ERROR = 3
EXIT_OUTDIR_EXISTS_ERROR = 4
#DEFAULT_MIN_LEN = 0
DEFAULT_VERBOSE = False
#HEADER = 'FILENAME\tNUMSEQ\tTOTAL\tMIN\tAVG\tMAX'
PROGRAM_NAME = "pydeel"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args(prefix):
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Pydeel: a tool to investigate bacterial or viral genome assembly based on protein lengths. Provide a fasta file and protein database as input and pydeel will provide gene completeness ratios'
    parser = argparse.ArgumentParser(description=description)
    # parser.add_argument('--minlen',
    #                     metavar='N',
    #                     type=int,
    #                     default=DEFAULT_MIN_LEN,
    #                     help='Minimum length sequence to include in stats (default {})'.format(DEFAULT_MIN_LEN))
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE, will be saved in outdir')
    parser.add_argument('-i', '--input',
                        required=True,
                        metavar='Path/to/input.fasta',
                        type=str,
                        help='File containing sequence, either in fasta format to be annotated, or pre-annotated gff or gbk')
    parser.add_argument('-c', '--code',
                        required = False,
                        default = 11,
                        type = int,
                        choices = range(1, 25),
                        metavar = 11,
                        help = 'Translation table for input sequence (default: 11)')
    parser.add_argument('-d', '--database',
                        metavar='Uniprot.dmnd',
                        type=str,
                        help='Protein database in diamond format')
    parser.add_argument('-o', '--outdir',
                        required = True,
                        type = str,
                        metavar = 'Path/to/output',
                        help = 'Name of output directory (required)')
    parser.add_argument('-t', '--title',
                        required = False,
                        default = 'pydeel',
                        type = str,
                        help='Prefix/title for files (default: "' + prefix + '-pydeel")')

    args = parser.parse_args()

    if args.code < 1 or args.code > 25:
        print('--code must be between 1 and 25 (inclusive)')
        exit
    else:
        args.code = str(args.code)

    # Change some arguments to full paths.

    args.outdir = os.path.abspath(args.outdir)

    if os.path.exists(args.database):
        args.database = os.path.abspath(args.database)
    else:
        print(args.database, 'does not exist!')
        exit


    if args.input:
        args.input = os.path.abspath(args.input)

    return args





def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))



def run_prodigal(input_file, gcode, protein_file):
    '''
    Runs the annotation tool prodigal on the input fasta file, saves the output
    to an .faa file
    '''
    print('Running prodigal on', input_file)
    print('Saving prodigal results to', protein_file)
    subprocess.run(["prodigal",
                    "-i", input_file,
                    "-g", gcode,
                    "-a", protein_file])


def make_diamond_db(database):
    subprocess.run(["diamond", "makedb",
                    "--in", database,
                    "-d", database])

def run_diamond(database, protein_file, lengths_file):
    format = str("6 slen qlen")
    print(format)
    print(format.strip('"\''))
    print("running:", "diamond", "blastp",
                    "--threads", "8",
                    "--max-target-seqs", "1",
                    "--db", database,
                    "--query", protein_file,
                    "--outfmt", "6 slen qlen",
                    "--out", lengths_file)
    command = "diamond blastp --threads 8 --max-target-seqs 1 --outfmt 6 qseqid qstart qend qlen sseqid sstart send slen nident pident evalue stitle qtitle --db " + database + " --query " + protein_file + " --out " + lengths_file
    subprocess.run(command, shell = True)
    #subprocess.run(["diamond", "blastp --outfmt 6 slen qlen", "--threads", "8", "--max-target-seqs", "1", "--db", database, "--query", protein_file, "--out", lengths_file], shell = True)


def convert_dataframe(lengths_data, convert_output):
    # import the diamond blast output file as a pandas dataframe, add headings
    df = pandas.read_csv(lengths_data,
                        sep = "\t",
                        names = ["qseqid", "qstart", "qend", "qlen", "sseqid", \
                        "sstart", "send", "slen", "nident", "pident", "evalue", \
                        "stitle", "qtitle"])
    # the 'slen' value is always one less than the query because database doesn't include stop codons
    df.head()
    df['slen'] = df['slen'] + 1

    # add a value for the query length divided by the sequence length)
    df['codingRatio'] = df['qlen'] / df['slen']

    # add a 'gene number', assumes proteins were searched in order
    df['geneNumber'] = df.index + 1

    #Finished with the file, time to write it back to the original location
    df.to_csv(convert_output, sep = "\t", header = True)

def plot_ratio(lengths_data, fullpath):
    lengths_data = pandas.read_csv(lengths_data,
                        sep = "\t")

    hist = altair.Chart(lengths_data)\
        .mark_bar(clip = True)\
        .encode(x = altair.X('codingRatio:Q',
                             bin = altair.Bin(step = 0.1),
                             scale = altair.Scale(domain=(0, 2)),
                             axis = altair.Axis(title='Query/Reference Ratio')
                            ),
                y = 'count()',
               tooltip = 'count()')\
        .configure_mark(
            fill = 'red',
            stroke = 'black')

    histzoom = altair.Chart(lengths_data)\
    .mark_bar(clip = True)\
    .encode(x = altair.X('codingRatio:Q',
                         bin = altair.Bin(step = 0.1),
                         scale = altair.Scale(domain=(0, 2)),
                         axis = altair.Axis(title='Query/Reference Ratio')
                        ),
            y = altair.Y('count()',
                        scale = altair.Scale(domain = (0,100))
                        ),
           tooltip = 'count()')\
    .configure_mark(
        fill = ' #c658dd ',
        stroke = 'black')

    genomeRatio = altair.Chart(lengths_data)\
    .mark_line(clip = True)\
    .encode(x = altair.X('geneNumber:Q',
                         scale = altair.Scale(domain = (0, len('geneNumber')))
                        ),
            y = altair.Y('codingRatio:Q',
                        scale = altair.Scale(type = 'log')),
            tooltip = 'codingRatio:Q')\
    .interactive()

    # save outputs
    print("saving image files")
    hist.save(fullpath + '-ratioplot-full' + '.html')
    histzoom.save(fullpath + '-ratioplot-zoom' + '.html')
    genomeRatio.save(fullpath + '-ratioplot-genome' + '.html')

    return None

def plot_ratio_seaborn(lengths_data, fullpath):
    lengths_data = pandas.read_csv(lengths_data,
                        sep = "\t")
    #Draw a plot of ratios
    seaborn.set_style(style = "ticks")
    hist = seaborn.distplot(lengths_data['codingRatio'],
                            hist = True,
                            hist_kws = {"color":"red"}
                            )

    fig = hist.get_figure()
    fig.savefig(fullpath + '-seaborn.png')

# =============================================================================
#     hist = altair.Chart(lengths_data)\
#         .mark_bar(clip = True)\
#         .encode(x = altair.X('codingRatio:Q',
#                              bin = altair.Bin(step = 0.1),
#                              scale = altair.Scale(domain=(0, 2)),
#                              axis = altair.Axis(title='Query/Reference Ratio')
#                             ),
#                 y = 'count()',
#                tooltip = 'count()')\
#         .configure_mark(
#             fill = 'red',
#             stroke = 'black')
#
# =============================================================================
    #save plot
    return None

def main():
    '''
    Orchestrate the execution of the program
    '''
    time=datetime.datetime.now()
    prefix = time.strftime("%Y%m%d-%H%M%S")
    options = parse_args(prefix)

    # # Create target Directory if don't exist
    # if not os.path.exists(options.outdir):
    #     os.mkdir(options.outdir)
    #     print("Directory " , options.outdir ,  " Created ")
    # else:
    #     exit_with_error(print("Directory" , options.outdir ,  "already exists"), EXIT_OUTDIR_EXISTS_ERROR)

    init_logging(options.log)


    #    if options.database is not None

    fullpath = options.outdir + "/" + options.title


    protein_file = fullpath + '.faa'

    print("input file:", options.input)
    print("output directory:", options.outdir)
    print("prefix:", options.title)
    print("genetic code:", options.code)

    if not os.path.exists(protein_file):
        run_prodigal(options.input, options.code, protein_file)
    else:
        print(protein_file, 'detected, skipping prodigal')

    lengths_file = fullpath + '.tsv'
    if not os.path.exists(lengths_file):
        run_diamond(options.database, protein_file, lengths_file)
    else:
        print(lengths_file, 'detected, skipping diamond blast')

    #error occurs if tsv has already been converted
    pandas_file = fullpath + '-pandas.tsv'
    if not os.path.exists(pandas_file):
        print("converting dataframe from diamond to pandas")
        convert_dataframe(lengths_file, pandas_file)
    else:
        print(pandas_file, 'detected, skipping data conversion')


    print("plotting coding ratios")
    plot_ratio(pandas_file, fullpath)
    plot_ratio_seaborn(pandas_file, fullpath)






# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()


# class FastaStats(object):
#     '''Compute various statistics for a FASTA file:
#
#     num_seqs: the number of sequences in the file satisfying the minimum
#        length requirement (minlen_threshold).
#     num_bases: the total length of all the counted sequences.
#     min_len: the minimum length of the counted sequences.
#     max_len: the maximum length of the counted sequences.
#     average: the average length of the counted sequences rounded down
#        to an integer.
#     '''
#     #pylint: disable=too-many-arguments
#     def __init__(self,
#                  num_seqs=None,
#                  num_bases=None,
#                  min_len=None,
#                  max_len=None,
#                  average=None):
#         "Build an empty FastaStats object"
#         self.num_seqs = num_seqs
#         self.num_bases = num_bases
#         self.min_len = min_len
#         self.max_len = max_len
#         self.average = average
#
#     def __eq__(self, other):
#         "Two FastaStats objects are equal if their attributes are equal"
#         if type(other) is type(self):
#             return self.__dict__ == other.__dict__
#         return False
#
#     def __repr__(self):
#         "Generate a printable representation of a FastaStats object"
#         return "FastaStats(num_seqs={}, num_bases={}, min_len={}, max_len={}, " \
#             "average={})".format(
#                 self.num_seqs, self.num_bases, self.min_len, self.max_len,
#                 self.average)
#
#     def from_file(self, fasta_file, minlen_threshold=DEFAULT_MIN_LEN):
#         '''Compute a FastaStats object from an input FASTA file.
#
#         Arguments:
#            fasta_file: an open file object for the FASTA file
#            minlen_threshold: the minimum length sequence to consider in
#               computing the statistics. Sequences in the input FASTA file
#               which have a length less than this value are ignored and not
#               considered in the resulting statistics.
#         Result:
#            A FastaStats object
#         '''
#         num_seqs = num_bases = 0
#         min_len = max_len = None
#         for seq in SeqIO.parse(fasta_file, "fasta"):
#             this_len = len(seq)
#             if this_len >= minlen_threshold:
#                 if num_seqs == 0:
#                     min_len = max_len = this_len
#                 else:
#                     min_len = min(this_len, min_len)
#                     max_len = max(this_len, max_len)
#                 num_seqs += 1
#                 num_bases += this_len
#         if num_seqs > 0:
#             self.average = int(floor(float(num_bases) / num_seqs))
#         else:
#             self.average = None
#         self.num_seqs = num_seqs
#         self.num_bases = num_bases
#         self.min_len = min_len
#         self.max_len = max_len
#         return self
#
#     def pretty(self, filename):
#         '''Generate a pretty printable representation of a FastaStats object
#         suitable for output of the program. The output is a tab-delimited
#         string containing the filename of the input FASTA file followed by
#         the attributes of the object. If 0 sequences were read from the FASTA
#         file then num_seqs and num_bases are output as 0, and min_len, average
#         and max_len are output as a dash "-".
#
#         Arguments:
#            filename: the name of the input FASTA file
#         Result:
#            A string suitable for pretty printed output
#         '''
#         if self.num_seqs > 0:
#             num_seqs = str(self.num_seqs)
#             num_bases = str(self.num_bases)
#             min_len = str(self.min_len)
#             average = str(self.average)
#             max_len = str(self.max_len)
#         else:
#             num_seqs = num_bases = "0"
#             min_len = average = max_len = "-"
#         return "\t".join([filename, num_seqs, num_bases, min_len, average,
#                           max_len])
#
#
#
#
#
#
#
# def process_files(options):
#     '''Compute and print FastaStats for each input FASTA file specified on the
#     command line. If no FASTA files are specified on the command line then
#     read from the standard input (stdin).
#
#     Arguments:
#        options: the command line options of the program
#     Result:
#        None
#     '''
#     if options.fasta_files:
#         for fasta_filename in options.fasta_files:
#             logging.info("Processing FASTA file from %s", fasta_filename)
#             try:
#                 fasta_file = open(fasta_filename)
#             except IOError as exception:
#                 exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
#             else:
#                 with fasta_file:
#                     stats = FastaStats().from_file(fasta_file, options.minlen)
#                     print(stats.pretty(fasta_filename))
#     else:
#         logging.info("Processing FASTA file from stdin")
#         stats = FastaStats().from_file(sys.stdin, options.minlen)
#         print(stats.pretty("stdin"))
