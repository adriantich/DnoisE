#!/usr/bin/env python3

"""
.. codeauthor:: Adrià Antich <adriantich@gmail.com>

This programme is called by the DnoisE.

import_data.py is designed to import fasta, fastq and csv files to be processed by DnoisE
Data is imported to pandas.DataFrame

"""

import io
import numpy as np
import pandas as pd
import sys


def import_data(de):
    print('start')

    print('reading input file')

    if de.input_type == 'fastq':
        try:
            input_file = pd.read_csv(de.MOTUfile, sep=' ', header=None)
        except:
            print("ERROR! incorrect FASTQ format file. FASTQ files must have space between header tags")
            sys.exit()
        try:
            if str(input_file.loc[0, 0]).find('@') != 0:
                print('ERROR! FASTQ files sequence headers must begin with "@"')
                if str(input_file.loc[0, 0]).find('@') < 0:
                    print('no @ found in first line')
                sys.exit()
            if str(input_file.loc[0]).find(de.count) < 0:
                print('ERROR! --count_name not found in FASTQ')
                sys.exit()
        except:
            print("ERROR! incorrect FASTQ header format")
            sys.exit()
        try:
            seqs = list(input_file.loc[list(range(1, input_file.shape[0], 4)), 0])
            ids = list(input_file.loc[list(range(0, input_file.shape[0], 4)), 0])
            if np.isnan(list(input_file.loc[0].str.contains(str(de.count + '=')))).sum() > 0:
                print("ERROR! incorrect FASTQ header format. No space allowed at the end of header")
                sys.exit()
            else:
                size = list(input_file.values[list(range(0, input_file.shape[0], 4)),
                                              list(input_file.loc[0].str.contains(str(de.count + '=')))])
            de.data_initial = pd.DataFrame({'id': ids, de.count: size, de.seq: seqs})
            de.data_initial = de.data_initial.replace(to_replace='@', value='', regex=True)
            de.data_initial = de.data_initial.replace(to_replace=(de.count + '='), value='', regex=True)
            de.data_initial[de.count] = pd.to_numeric(de.data_initial[de.count])
            del input_file, seqs, ids, size
        except:
            print("ERROR! incorrect FASTQ format file. ")
            sys.exit()
    elif de.input_type == 'fasta':
        try:
            file = open(de.MOTUfile, mode='r')
            all_file = file.read()
            file.close()
            all_file = all_file.replace('\r', '')
            all_file = all_file.replace('\n', '')
            all_file = all_file.replace('>', '\n>')
            io_file = io.StringIO(all_file)
            de.data_initial = pd.read_csv(io_file, sep=';', header=None)
            del all_file, io_file
            de.data_initial.columns = ['id', de.count, de.seq]
            de.data_initial = de.data_initial.replace(to_replace='>', value='', regex=True)
            de.data_initial = de.data_initial.replace(to_replace=(de.count + '='), value='', regex=True)
            de.data_initial[de.count] = pd.to_numeric(de.data_initial[de.count])
        except:
            print("ERROR! incorrect FASTA format file")
            sys.exit()

    elif de.input_type == 'csv':
        de.data_initial = pd.read_csv(de.MOTUfile, sep=de.sep)

    if de.merge_from_info:
        de.merge_data = pd.read_csv(de.infofile, sep=',')

    print('input file read')
