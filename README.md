[![travis](https://travis-ci.org/alegione/pydeel.svg?branch=master)](https://travis-ci.org/alegione/pydeel)

# Overview

This program is a basic python conversion of Mick Watson's [Ideel](https://github.com/mw55309/ideel). It reads an input FASTA file and uses prodigal for rapid annotation, then runs diamond blast on the output, then compares the query length to hit length.

It was built with the help of '[Bionitio](https://github.com/bionitio-team/bionitio)'


# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/alegione/pydeel/master/LICENSE).

# Installing

Clone this repository:
```
git clone https://github.com/alegione/pydeel
```

Pydeel can be installed using `pip` in a variety of ways:

1. Inside a virtual environment:
```
python3 -m venv pydeel_dev
source pydeel_dev/bin/activate
pip install -U /path/to/pydeel
```
2. Into the global package database for all users:
```
pip install -U /path/to/pydeel
```
3. Into the user package database (for the current user only):
```
pip install -U --user /path/to/pydeel
```


# General behaviour

Pydeel is a basic python conversion of Mick Watson's Ideel. It reads one or more input FASTA files and for each file it will use prodigal for rapid annotation, then run diamond blast, then compare the query length to hit length.

## Help message

Pydeel can display usage information on the command line via the `-h` or `--help` argument:

```
pydeel.py -h
usage: pydeel.py [-h] [--version] [--log LOG_FILE] -i Path/to/input.fasta
                 [-c 11] [-d Uniprot.dmnd] -o Path/to/output [-t TITLE]

Pydeel: a tool to investigate bacterial or viral genome assembly based on
protein lengths. Provide a fasta file and protein database as input and pydeel
will provide gene completeness ratios

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --log LOG_FILE        record program progress in LOG_FILE, will be saved in outdir
  -i Path/to/input.fasta, --input Path/to/input.fasta
                        File containing sequence, either in fasta format to be
                        annotated, or pre-annotated gff or gbk
  -c 11, --code 11      Translation table for input sequence (default: 11)
  -d Uniprot.dmnd, --database Uniprot.dmnd
                        Protein database in diamond format
  -o Path/to/output, --outdir Path/to/output
                        Name of output directory (required)
  -t TITLE, --title TITLE
                        Prefix/title for files (default: "YYYYMMDD-hhmmss-pydeel")
```

## Logging

If the ``--log FILE`` command line argument is specified, pydeel will output a log file containing information about program progress. The log file includes the command line used to execute the program, and a note indicating which files have been processes so far. Events in the log file are annotated with their date and time of occurrence.

# Exit status values

Pydeel returns the following exit status values:

* 0: The program completed successfully.
* 1: File I/O error. This can occur if at least one of the input FASTA files cannot be opened for reading. This can occur because the file does not exist at the specified path, or pydeel does not have permission to read from the file.
* 2: A command line error occurred. This can happen if the user specifies an incorrect command line argument. In this circumstance pydeel will also print a usage message to the standard error device (stderr).
* 3: Input FASTA file is invalid. This can occur if pydeel can read an input file but the file format is invalid.
* 4: Output directory exists already and cannot be overwritten

# Testing

## Unit tests

```
cd pydeel/python/pydeel
python -m unittest -v pydeel_test
```

## Test suite

A set of sample test input files is provided in the `test_data` folder.

# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[pydeel issue tracker](https://github.com/alegione/pydeel/issues)
