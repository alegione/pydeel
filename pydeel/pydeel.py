'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Alistair Legione, 11 Jun 2019
License     : MIT
Maintainer  : legionea@unimelb.edu.au
Portability : POSIX

This program is a basic python conversion of Mick Watson's Ideel.
It reads one or more input FASTA files and for each file it will use
prodigal for rapid annotation, then run diamond blast, then compare the
query length to hit length.

It was built with the help of 'Bionitio'
'''

import altair
import argparse
from Bio import SeqIO
import numpy
import statsmodels.stats.api as sms
import datetime
import logging
import os
import pandas
import pkg_resources
import seaborn
#import selenium # doesn't need to be imported??
import subprocess
import sys

#Local modules
#import Convert2fasta

EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_FASTA_FILE_ERROR = 3
EXIT_OUTDIR_EXISTS_ERROR = 4
DEFAULT_VERBOSE = False
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
    parser.add_argument('--version',
                        action = 'version',
                        version = '%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('-i', '--input',
                        required = True,
                        metavar = 'Path/to/input.fasta',
                        type = str,
                        help = 'File (or directory) containing sequence, either in fasta format to be annotated or pre-annotated faa')
    parser.add_argument('-c', '--code',
                        required = False,
                        default = 11,
                        type = int,
                        choices = range(1, 25),
                        metavar = 11,
                        help = 'Translation table for input sequence (default: 11)')
    parser.add_argument('-d', '--database',
                        metavar = 'Uniprot.dmnd',
                        type = str,
                        default = None,
                        help = 'Protein database in diamond format')
    parser.add_argument('-o', '--outdir',
                        required = True,
                        type = str,
                        metavar = 'Path/to/output',
                        help = 'Name of output directory (required)')
    parser.add_argument('-n', '--name',
                        required = False,
                        default = 'pydeel',
                        type = str,
                        help = 'Prefix/name for files (default: "YYYYMMDD-hhmmss-pydeel")')
    parser.add_argument('-p', '--proteins',
                        metavar = 'Path/to/RefProtein.faa',
                        required = False,
                        type = str,
                        default = None,
                        help = 'Input protein reference to compare annotations against'
                        )
    parser.add_argument('-f','--force',
                        required = False,
                        action='store_true',
                        help = 'Overwrite any directories/files with the same names present at the target'
                        )
    parser.add_argument('-r', '--resume',
                        required = False,
                        action = 'store_true',
                        help = 'Continue run where last completed')
    parser.add_argument('--single',
                        required = False,
                        type = float,
                        default = 0.90,
                        help = 'Report pass/fail of single gene ratios'
                        )
    parser.add_argument('--group',
                        required = False,
                        type = float,
                        default = 0.20,
                        help = 'Report pass/fail of genome based on proportion of genes that fail individual test'
                        )

# TODO: add option to NOT output png files directory to avoid need for selenium and chromedriver
    args = parser.parse_args()

    if args.code < 1 or args.code > 25:
        Error = '--code must be between 1 and 25 (inclusive)'
        exit_with_error(Error, EXIT_COMMAND_LINE_ERROR)
    else:
        args.code = str(args.code)

    # Change some arguments to full paths.

    args.outdir = os.path.abspath(args.outdir)

    if args.database is not None and args.proteins is not None:
        Error = "Please only provide one of --database or --proteins, not both"
        exit_with_error(Error, EXIT_COMMAND_LINE_ERROR)
    elif args.database is not None and os.path.exists(args.database):
        args.database = os.path.abspath(args.database)
        print("database:", args.database)
    elif args.proteins is not None and os.path.exists(args.proteins):
        args.proteins = os.path.abspath(args.proteins)
        print("proteins:", args.proteins)
    else:
        if args.database is not None:
            Error = args.database + " does not exist!"
            exit_with_error(Error, EXIT_COMMAND_LINE_ERROR)
        elif args.proteins is not None:
            Error = args.proteins + " does not exist!"
            exit_with_error(Error, EXIT_COMMAND_LINE_ERROR)
        else:
            Error = "I'm not sure how I got here"
            exit_with_error(Error, EXIT_COMMAND_LINE_ERROR)

    if args.input:
        args.input = os.path.abspath(args.input)

    return args

def init_logging(log_filename):
    '''Initialise the logging facility, and write log statement
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
        logging.getLogger().addHandler(logging.StreamHandler())

        logging.info('Program started')
        logging.info('Command line: %s', ' '.join(sys.argv))


def is_fasta(filename):
    '''
    taken from a helpful stackoverflow comment (https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta)
    '''
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file



def run_prodigal(input_file, gcode, protein_file, output_file):
    '''
    Runs the annotation tool prodigal on the input fasta file, saves the output
    to an .faa file
    '''
    logging.info('Running prodigal on %s', input_file)
    subprocess.run(["prodigal",
                    "-i", input_file,
                    "-g", gcode,
                    "-a", protein_file,
                    "-o", output_file])
    ## TODO: add a check here to see if the file was actually made
    logging.info('Saved prodigal results to %s', protein_file)


    return None

def make_diamond_db(proteins, output_file):
    logging.info("Building diamond database from %s protein database", proteins)

    subprocess.run(["diamond", "makedb",
                    "--in", proteins,
                    "--db", output_file,
                    "--quiet"])
    ## TODO: add a check here to see if the file was actually made
    logging.info("Saved database to %s", output_file)
    return None

def run_diamond(database, protein_file, lengths_file):
    logging.info("Running diamond on %s against %s to generate lengths ratio", protein_file, database)
    dmnd_command = "diamond blastp --quiet --threads 8 --max-target-seqs 1 --outfmt 6 qseqid qstart qend qlen sseqid sstart send slen nident pident evalue stitle qtitle --db " \
                    + database \
                    + " --query " + protein_file + \
                    " --out " + lengths_file
    logging.info("diamond command: %s", str(dmnd_command))
    subprocess.run(dmnd_command, shell = True)
    return None

def convert_dataframe(lengths_data, convert_output, single_threshold):
    # import the diamond blast output file as a pandas dataframe, add headings
    logging.info('Loading diamond results to pandas dataframe and outputting %s', convert_output)

    df = pandas.read_csv(lengths_data,
                        sep = "\t",
                        names = ["qseqid", "qstart", "qend", "qlen", "sseqid", \
                        "sstart", "send", "slen", "nident", "pident", "evalue", \
                        "stitle", "qtitle"])

    # the 'slen' value is always one less than the query because database doesn't include stop codons
    df['slen'] = df['slen'] + 1

    # add a value for the query length divided by the sequence length)
    df['codingRatio'] = df['qlen'] / df['slen']

    # add a pass/fail to each ORF depending on the provided threshold
    df = df.assign(ORF_pass = numpy.where(df['codingRatio'] >= single_threshold, 'yes', 'no'))

    # add a 'gene number', assumes proteins were searched in order
    df['geneNumber'] = df.index + 1

    #Finished with the file, time to write it back to the original location
    df.to_csv(convert_output, sep = "\t", header = True)

    return None

def plot_ratio(lengths_data, fullpath):
    logging.info("Generating plots")

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
    .mark_line()\
    .encode(x = altair.X('geneNumber:Q',
                         scale = altair.Scale(domain = (0, len('geneNumber')))
                        ),
            y = altair.Y('codingRatio:Q',
                        scale = altair.Scale(type = 'log')),
            tooltip = 'codingRatio:Q')\
    .interactive()

    # violin = altair.Chart(lengths_data) \
    #     .transform_bin(['bin_max', 'bin_min'], field='codingRatio', bin=altair.Bin(maxbins=20)) \
    #     .transform_calculate(binned=(altair.datum.bin_max + altair.datum.bin_min) / 2) \
    #     .transform_aggregate(value_count='count()', groupby=['qseqid', 'binned']) \
    #     .transform_impute(impute='value_count', groupby=['qseqid'], key='binned', value=0) \
    #     .mark_area( interpolate='monotone',
    #                 orient='horizontal') \
    #     .encode(x=altair.X('value_count:Q',
    #                         title=None,
    #                         stack='center',
    #                         axis=altair.Axis(labels=False, values=[0],grid=False, ticks=True),
    #                         ),
    #             y=altair.Y('binned:Q', bin='binned', title='Coding Ratio'),
    #             color=altair.Color('qseqid:N', legend=None),
    #             column=altair.Column('qseqid:N',
    #                 header=altair.Header(titleOrient='bottom',
    #                                     labelOrient='bottom',
    #                                     labelPadding=0,
    #                                     ),
    #                                 )) \
    #     .interactive()

    # save outputs
    logging.info("Saving image files to html")
    hist.save(fullpath + '-ratioplot-full' + '.html')
    histzoom.save(fullpath + '-ratioplot-zoom' + '.html')
    genomeRatio.save(fullpath + '-ratioplot-genome' + '.html')
    # violin.save(fullpath + '-codingRatioViolin' + '.html')

#    print("Saving image files to png")
#    hist.save(fullpath + '-ratioplot-full' + '.png', webdriver='firefox')
#    histzoom.save(fullpath + '-ratioplot-zoom' + '.png', webdriver='firefox')
#    genomeRatio.save(fullpath + '-ratioplot-genome' + '.png', webdriver='firefox')

    return None

def plot_ratio_seaborn(lengths_data, fullpath):
    lengths_data = pandas.read_csv(lengths_data,
                        sep = "\t")
    #Draw a plot of ratios
    seaborn.set_style(style = "ticks")
    ratio = lengths_data[lengths_data['codingRatio'] <= 2]
    hist = seaborn.distplot(ratio['codingRatio'],
                            hist = True,
                            hist_kws = {"color":"red"},
                            kde = False,
                            rug = True,
                            bins = 10
                            )

    fig = hist.get_figure()
    fig.savefig(fullpath + '-seaborn.png')
    return None


#TODO: Basic stats that can also be used for unit testing


def genome_stats(lengths_data, single_threshold, group_threshold):
    #Run some stats?
    title = os.path.splitext(os.path.basename(lengths_data))[0]
    lengths_data = pandas.read_csv(lengths_data,
                        sep = "\t")
    # full length divided by 'count_orfs'

    count_orfs = len(lengths_data)
    average_ratio = numpy.average(lengths_data['codingRatio'])
    perfect_orfs_count = (lengths_data.codingRatio.values == 1).sum()
    perfect_orfs_percent = perfect_orfs_count / count_orfs
    below_singleThres_count = len(lengths_data[lengths_data.ORF_pass.str.contains("no")])
    below_singleThres_percent = below_singleThres_count / count_orfs
    q25, median, q75 = numpy.percentile(lengths_data.codingRatio, [25, 50, 75])
    min = numpy.min(lengths_data.codingRatio)
    max = numpy.max(lengths_data.codingRatio)
    std = numpy.std(lengths_data.codingRatio)
    lower, upper = sms.DescrStatsW(lengths_data.codingRatio).tconfint_mean()
    if below_singleThres_percent > group_threshold:
        pass_group_threshold = False
    else:
        pass_group_threshold = True
    results = [title, count_orfs, perfect_orfs_count, perfect_orfs_percent * 100, \
                below_singleThres_count, below_singleThres_percent * 100, pass_group_threshold, \
                average_ratio, std, lower, upper, min, q25, median, q75, max]
    print(results)
    return results




def main():
    '''
    Orchestrate the execution of the program
    '''
    # some initial variables set
    time = datetime.datetime.now()
    prefix = time.strftime("%Y%m%d-%H%M%S")
    options = parse_args(prefix)
    logfile = options.outdir + "/" + options.name + ".log"

    # Replace target directory if it exists and 'force' set
    if os.path.exists(options.outdir) and options.force == True:
        import shutil
        shutil.rmtree(options.outdir, ignore_errors=True)
        print("Directory", options.outdir,  "will be replaced")

    # Create target directory
    if not os.path.exists(options.outdir):
         os.mkdir(options.outdir)
         init_logging(logfile)

         logging.info("Directory %s created", options.outdir)
    elif os.path.exists(options.outdir) and options.resume == False:
         exit_with_error(logging.info("Directory %s already exists", options.outdir), EXIT_OUTDIR_EXISTS_ERROR)
    else:
        logging.info("Directory %s exists. Resuming workflow", options.outdir)


    if options.proteins is not None:
        ref_protein = options.proteins
        ref_database = str(ref_protein.split('.', 1)[0]) + '.dmnd'
        if not os.path.exists(ref_database) or options.force == True:
            make_diamond_db(ref_protein, ref_database)
    else:
        ref_database = options.database


    if os.path.isdir(options.input) == True:
        input_data = os.listdir(options.input)
        multi = True
    else:
        input_data = [options.input]
        multi = False

    annotation_stats = options.outdir + "/" + options.name + "_annotation_stats.tsv"

    for sequence in input_data:
        if multi == True:
            sequence = os.path.abspath(options.input + "/" + sequence)

        else:
            sequence = os.path.abspath(options.input)
        # check if input is a fasta file
        if is_fasta(sequence) == False:
            # If an error is thrown here it won't work with a folder with multiple types in it?
            # exit_with_error('The input must be in fasta format', EXIT_FASTA_FILE_ERROR)
            continue
        #TODO: check if file is FASTA or gbk and convert if necessary

        full_outpath = options.outdir + "/" + options.name + "_" + os.path.splitext(os.path.basename(sequence))[0]

        protein_file = full_outpath + '.faa'

        if not os.path.exists(protein_file) or options.force == True:
            logging.info('Running prodigal on %s using genetic code #%s', sequence, options.code)
            prodigal_output = full_outpath + '.prodigal.txt'
            run_prodigal(sequence, options.code, protein_file, prodigal_output)
        else:
            logging.info("%s detected, skipping prodigal", protein_file)

        lengths_file = full_outpath + '.tsv'
        if not os.path.exists(lengths_file) or options.force == True:
            logging.info('Running Diamond BLAST on generated open reading frames')
            run_diamond(ref_database, protein_file, lengths_file)
        else:
            logging.info("%s detected, skipping diamond blast", lengths_file)

        #error occurs if tsv has already been converted
        pandas_file = full_outpath + '-pandas.tsv'
        if not os.path.exists(pandas_file) or options.force == True:
            logging.info("converting dataframe from diamond to pandas")
            convert_dataframe(lengths_file, pandas_file, options.single)
        else:
            logging.info("%s detected, skipping data conversion", pandas_file)

            # TODO need to set up to overwrite old file if running with --force, if added to the if statement it overwrites every loop
        if not os.path.exists(annotation_stats):
            logging.info("Creating statistics file")
            genome_stat_df = pandas.DataFrame([genome_stats(pandas_file, options.single, options.group)], columns = ['title', 'count_orfs', 'perfect_orfs', 'perfect_orfs_%', 'below_' + str(options.single), 'below_' + str(options.single) + '_%', 'Pass_' + str(options.group) + '_threshold', \
                        'average_ratio', 'Standard Deviation', '95% CI lower', '95% CI upper', 'Minimum', 'Q25', 'median', 'Q75', 'Max'])
            genome_stat_df.to_csv(annotation_stats, index = False, sep = '\t', columns = ["title", "count_orfs", "perfect_orfs", "perfect_orfs_%", "below_" + str(options.single), "below_" + str(options.single) + "_%", "Pass_" + str(options.group) + "_threshold", \
                        "average_ratio", "Standard Deviation", "95% CI lower", "95% CI upper", "Minimum", "Q25", "median", "Q75", "Max"])
        else:
            logging.info("Appending annotation stats to existing file")
            genome_stat_df = pandas.DataFrame([genome_stats(pandas_file, options.single, options.group)], columns = ['title', 'count_orfs', 'perfect_orfs', 'perfect_orfs_%', 'below_' + str(options.single), 'below_' + str(options.single) + '_%', 'Pass_' + str(options.group) + '_threshold', \
                        'average_ratio', 'Standard Deviation', '95% CI lower', '95% CI upper', 'Minimum', 'Q25', 'median', 'Q75', 'Max'])
            old_file = pandas.read_csv(annotation_stats, sep = '\t')
            new_file = old_file.append(genome_stat_df, ignore_index = True, sort = False)
            new_file.to_csv(annotation_stats, sep = '\t', index = False)


        logging.info("plotting coding ratios")
        if not os.path.exists(full_outpath + '-ratioplot-genome' + '.png') or options.force == True:
            plot_ratio(pandas_file, full_outpath)
            plot_ratio_seaborn(pandas_file, full_outpath)




# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
