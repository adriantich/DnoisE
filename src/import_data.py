#!/usr/bin/env python3

"""
.. codeauthor:: Adrià Antich <adriantich@gmail.com>

This programme is called by the DnoisE.

import_data.py is designed to import fasta, fastq and csv files to be processed by DnoisE
Data is imported to pandas.DataFrame

"""

import pandas as pd


def import_data(de):
    print('start')
    if de.part == 1:
        print('reading input file')
        if de.fasta:
            input_file = pd.read_csv(de.MOTUfile, sep=';', header=None)
            seqs = list(input_file.loc[list(range(1, input_file.shape[0], 2)), 0])
            ids = list(input_file.loc[list(range(0, input_file.shape[0], 2)), 0])
            size = list(input_file.loc[list(range(0, input_file.shape[0], 2)), 1])
            de.data_initial = pd.DataFrame({'id': ids, de.count: size, de.seq: seqs})
            de.data_initial = de.data_initial.replace(to_replace='>', value='', regex=True)
            de.data_initial = de.data_initial.replace(to_replace=(de.count + '='), value='', regex=True)
            de.data_initial[de.count] = pd.to_numeric(de.data_initial[de.count])
            del input_file, seqs, ids, size

        else:
            de.data_initial = pd.read_csv(de.MOTUfile, sep=de.sep)
    else:
        de.read_variables2()

    print('input file read')
