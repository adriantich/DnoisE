#!/usr/bin/env python3

"""
.. codeauthor:: Adrià Antich <adriantich@gmail.com>

This programme is called by the DnoisE.

denoise_functions.py creates a class object which contains many of the functions used by DnoisE and objects used by the
programm.

"""

import getopt
import itertools
import Levenshtein as lv
import multiprocessing as mp
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm


class DnoisEFunctions:
    data_initial = pd.DataFrame()
    run_list = []
    run_list_2 = []
    alpha = 5
    Ad1 = 0.8379801130824722
    Ad2 = 0.357379606161045
    Ad3 = 0.8379801130824722
    max_ratio = (1 / 2) ** (alpha * 1 * min(Ad1, Ad2, Ad3) + 1)
    min_mother = 1/max_ratio
    MOTUfile = ""
    good_seq = []
    output_info = []
    denoised_d_output = []
    denoised_ratio_output = []
    denoised_ratio_d_output = []
    good_seq_2 = []
    output_info_2 = []
    denoised_d_output_2 = []
    denoised_ratio_output_2 = []
    denoised_ratio_d_output_2 = []
    entropy = True
    cores = 1
    MOTUoutfile = ""
    seq = 'sequence'
    count = 'size'
    q1 = 0
    q2 = 0
    q3 = 0
    q4 = 0
    q5 = 0
    q6 = 0
    q7 = 0
    q8 = 0
    q9 = 0
    q10 = 0
    initial_pos = 3
    justcount = True
    abund_col_names = []
    first_col_names = []
    output_file_type = ''
    input_type = ''
    sep = '\t'
    end = 1
    start = 1
    output_type = 'ratio_d'
    part = 1
    good_mothers = []
    new_output_part_2 = False
    new_fasta_output_part_2 = False
    modal_length_value = []
    compute_entropy = False
    infofile = []
    merge_from_info = False
    unique_length = False

    def __init__(self):
        print("starting to denoise")

    def read_parameters(self, argument_list):
        short_options = "hP:f:F:j:c:s:z:n:a:q:p:e:yx:m:u"
        long_options = ["help", "fasta_input=", "csv_input=", "fastq_input=",
                        "fasta_output=", "csv_output=", "part=", "joining_criteria=",
                        "cores=", "start_sample_cols=",
                        "end_sample_cols=", "count_name=", "alpha=", "sequence=", "separation=", "entropy=",
                        "entropy_correction", "first_nt_codon_position=", "modal_length=", "joining_file=",
                        "unique_length"]
        try:
            arguments, values = getopt.getopt(argument_list, short_options, long_options)
        except getopt.error as err:
            # Output error, and return with an error code
            print(str(err))
            sys.exit(2)
        for current_argument, current_value in arguments:
            # help
            if current_argument in ("-h", "--help"):
                print("\033[1mDisplaying help\033[0m\n"
                      "\t\t-h --help display help\n"
                      "\t\033[1mInput file options:\033[0m\n"
                      "\t\t--csv_input [path] input file path in csv format\n"
                      "\t\t--fasta_input [path] input file path in fasta format\n"
                      "\t\t--fastq_input [path] input file path in fastq format\n"
                      "\t\t--joining_file [path] file path of an info output from DnoisE. This option "
                      "allows to use the information of previous runs of DnoisE to return different joining criteria"
                      "outputs without running all the programm again\n"
                      "\t\t-n --count_name [size/reads/count...] count name column 'size' by default\n"
                      "\t\t-p --sep [1/2/3] separation in case of csv input file\n"
                      "\t\t\t\t1='\t' (tab)\n"
                      "\t\t\t\t2=','\n"
                      "\t\t\t\t3=';'\n"
                      "\t\t-q --sequence [sequence/seq...] sequence column name, 'sequence' by default\n"
                      "\t\t-s --start_sample_cols [number] first sample column (1 == 1st col) if not given, "
                      "just one column with total read count expected (see README.md)\n"
                      "\t\t-z --end_sample_cols [number] last sample column (n == nst col) if not given, "
                      "just one column with total read count expected (see README.md)\n"
                      "\t\033[1mOutput file options:\033[0m\n"
                      "\t\t--csv_output [path] common path for csv format\n"
                      "\t\t--fasta_output [path] common path for fasta format\n"
                      "\t\t-j --joining_criteria [1/2/3]\n"
                      "\t\t\t\t1-> will join by the lesser [abundance ratio / beta(d)]\n"
                      "\t\t\t\t2-> will join by the lesser abundance ratio (r criterion)\n"
                      "\t\t\t\t3-> will join by the lesser d value (d criterion)\n"
                      "\t\t\t\t4-> will provide all joining criteria in three "
                      "different outputs (r_d criterion) (default)\n"
                      "\t\033[1mOther options:\033[0m\n"
                      "\t\t-a --alpha [number] alpha value, 5 by default\n"
                      "\t\t-c --cores [number] number of cores, 1 by default\n"
                      "\t\t-e --entropy [number,number,number] entropy values (or any user-settable "
                      "measure of variability) of the different codon positions [0.47,0.23,1.02] by default\n"
                      "\t\t-m --modal_length [number] when running DnoisE with entropy correction, "
                      "sequence length expected can be set, if not, modal_length is used and sequences"
                      " with modal_length + or - 3*n are accepted\n"
                      "\t\t-u --unique_length only modal length is accepted as sequence length when running with "
                      "entropy correction\n"
                      "\t\t-x --first_nt_codon_position [number] as DnoisE has been developed for COI "
                      "sequences amplified with Leray-XT primers, default value is 3\n"
                      "\t\t-y --entropy_correction compute a distance correction "
                      "based on entropy is performed (see ENTROPY CORRECTION below). If set to F, "
                      "no correction for entropy is performed (corresponding to the standard Unoise formulation)\n")
                sys.exit()
            # input args
            elif current_argument == "--csv_input":
                print("Denoising %s file" % current_value)
                self.MOTUfile = current_value
                self.input_type = 'csv'
                arg_ci = True
            elif current_argument == "--fasta_input":
                print("Denoising %s file" % current_value)
                self.MOTUfile = current_value
                self.input_type = 'fasta'
                arg_fi = True
            elif current_argument == "--fastq_input":
                print("Denoising %s file" % current_value)
                self.MOTUfile = current_value
                self.input_type = 'fastq'
                arg_fqi = True
            elif current_argument == "--joining_file":
                print("Merging sequences from %s file info" % current_value)
                self.infofile = current_value
                self.merge_from_info = True
                arg_ii = True
            elif current_argument in ("-n", "--count_name"):
                self.count = current_value
                arg_n = True
                print("Count col name: %s" % current_value)
            elif current_argument in ("-p", "--sep"):
                print(current_value)
                if current_value == "1":
                    self.sep = '\t'
                if current_value == "2":
                    self.sep = ','
                if current_value == "3":
                    self.sep = ';'
                arg_p = True
                print("Sep: %s" % self.sep)
            elif current_argument in ("-q", "--sequence"):
                self.seq = current_value
                arg_q = True
                print("Seq: %s" % current_value)
            elif current_argument in ("-s", "--start_sample_cols"):
                self.start = int(current_value)
                arg_s = True
                print("Abundant cols starts in: %s" % current_value)
            elif current_argument in ("-z", "--end_sample_cols"):
                self.end = int(current_value)
                arg_z = True
                print("Abundant cols ends in: %s" % current_value)
            # output args
            elif current_argument == "--csv_output":
                print("Output files will be named %s*" % current_value)
                self.MOTUoutfile = current_value
                self.output_file_type = 'csv'
                arg_co = True
            elif current_argument == "--fasta_output":
                print("Output files will be named %s*" % current_value)
                self.MOTUoutfile = current_value
                self.output_file_type = 'fasta'
                arg_fo = True
            elif current_argument in ("-j", "--joining_criteria"):
                if current_value == "1":
                    self.output_type = 'ratio_d'
                elif current_value == "2":
                    self.output_type = 'ratio'
                elif current_value == "3":
                    self.output_type = 'd'
                elif current_value == "4":
                    self.output_type = 'all'
                arg_j = True
                print("output file: %s" % self.output_type)
            # other args
            elif current_argument in ("-a", "--alpha"):
                self.alpha = int(current_value)
                arg_a = True
                print("Alpha: %s" % current_value)
            elif current_argument in ("-c", "--cores"):
                print("Running with %s cores" % current_value)
                self.cores = int(current_value)
                arg_c = True
            elif current_argument in ("-e", "--entropy"):
                e1, e2, e3 = current_value.split(",")
                print(str("E1: " + e1 + " \nE2: " + e2 + " \nE3: " + e3))
                e1 = float(e1)
                e2 = float(e2)
                e3 = float(e3)
                self.entropy = True
                arg_e = True
            elif current_argument in ("-m", "--modal_length"):
                self.modal_length_value = [int(current_value)]
                print("modal_length set as %s" % current_value)
            elif current_argument in ("-u", "--unique_length"):
                self.unique_length = True
                arg_u = True
                print("Is entropy taken into account: %s" % self.entropy)
            elif current_argument in ("-x", "--first_nt_codon_position"):
                self.initial_pos = int(current_value)
                arg_x = True
                print("first nt is a position %s" % current_value)
            elif current_argument in ("-y", "--entropy_correction"):
                self.entropy = True
                arg_y = True
                print("Is entropy taken into account: %s" % self.entropy)
            elif current_argument in ("-P", "--Part"):
                print("Part %s running" % current_value)
                self.part = int(current_value)
                arg_pp = True

        if 'arg_ci' not in locals() and 'arg_fi' not in locals() and 'arg_fqi' not in locals():  # no input file
            print("Err: input file needed")
            sys.exit()
        if 'arg_ii' not in locals():
            self.merge_from_info = False
        if "arg_n" not in locals():  # no count name
            print("count_name not given, 'size' by default")
            self.count = "size"
        if self.input_type == 'csv':  # input file as csv
            if 'arg_p' not in locals():  # no sep specified
                print("Separation not given, '\t' by default")
                self.sep = '\t'
            if "arg_q" not in locals():  # no sequence tag name
                print("sequence tag name not given, 'sequence' by default")
                self.seq = 'sequence'
            if "arg_s" not in locals():  # no start of sample cols specified
                if "arg_z" not in locals():  # no end of sample cols specified
                    print("start and end not given, no samples in file")
                    self.abund_col_names = [self.count]
                    self.justcount = True
                else:
                    print("Err: end given but no start")
                    sys.exit()
            elif "arg_z" not in locals():  # no end of sample cols specified
                print("Err: start given but no end")
                sys.exit()
            else:  # sample cols specified
                self.justcount = False
        if 'arg_co' not in locals() and 'arg_fo' not in locals():  # no output file
            print("Err: output path needed")
            sys.exit()
        if 'arg_j' not in locals():  # no joining criteria specified
            self.output_type = 'ratio_d'
        if "arg_a" not in locals():  # alpha not specified
            if not self.merge_from_info:
                print("alpha not given, 5 by default")
            self.alpha = 5
        if 'arg_c' not in locals():
            print("cores not given, 1 core by default")
            self.cores = 1
        if "arg_y" not in locals() and "arg_e" not in locals():  # no entropy correction
            self.entropy = False
            if not self.merge_from_info:
                print("Ad correction not applied")
            self.Ad1, self.Ad2, self.Ad3 = (1, 1, 1)
        else:  # entropy correction
            if "arg_e" not in locals():
                self.compute_entropy = True
                if not self.merge_from_info:
                    print("Entropy will be computed from data")
            # defining Ad correction factor taking into account Entropy
            else:
                self.Ad1 = e1 / (e1 + e2 + e3) * 3
                self.Ad2 = e2 / (e1 + e2 + e3) * 3
                self.Ad3 = e3 / (e1 + e2 + e3) * 3
            if 'arg_x' not in locals():
                self.initial_pos = 3
            if 'arg_u' not in locals():
                self.unique_length = False

    def quartiles_run(self, data_initial):
        q = data_initial.shape[0]

        self.q1 = int(q * 1 / 10)
        self.q2 = int(q * 2 / 10)
        self.q3 = int(q * 3 / 10)
        self.q4 = int(q * 4 / 10)
        self.q5 = int(q * 5 / 10)
        self.q6 = int(q * 6 / 10)
        self.q7 = int(q * 7 / 10)
        self.q8 = int(q * 8 / 10)
        self.q9 = int(q * 9 / 10)
        self.q10 = q - 1

    def denoising(self, pos):
        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]
        position = len(self.run_list)

        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is daughter break to the next pM
            if self.run_list[a].get('daughter'):
                continue
            pM = self.run_list[a].get('id')
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > self.max_ratio:
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d = lv.distance(pDseq, pMseq)
            # if d:
            # if dd exist & d is higher than dd:
            # break go to next pM
            if 'dd' in locals():
                if dd == 1:
                    break
                if d >= dd:
                    continue
                # if dd doesn't exist | d is smaller than dd:
                # dd = d
                else:
                    dd = d
            else:
                dd = d
                # if Edgar's equation:
            if b_ratio <= (1 / 2) ** (self.alpha * d + 1):
                # TRUE:
                # add to Ml:
                xavier_criteria = b_ratio / ((1 / 2) ** (self.alpha * d + 1))
                df1 = [
                    {'pM': pM, 'pMpos': pMpos, 'pD': pD, 'ratio': b_ratio, 'd': d, 'xavier_criteria': xavier_criteria}]
                Ml = Ml.append(df1)
                # identification of the pM
                # ratio
                # d
                # Xavier's criteria ---> ratio/((1/2)^(alpha*d + 1))
                # (the less value the best pM)
                # FALSE ---> break this pM
            else:
                continue

            # editing variables and files
        if Ml.empty:  # means that is not a daughter
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_ratio_d': None, 'xavier_criteria': None}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': False}
            return [True], [info], [pD], [pD], [pD], [run_list]
        else:  # it is a daughter
            # print pD name to each pM depending on different criteria
            # _mothers_d
            pM_d = Ml.loc[(Ml['d'] == min(Ml.loc[:, 'd'])), 'pM'][0]
            pM_ratio = Ml.loc[(Ml['ratio'] == min(Ml.loc[:, 'ratio'])), 'pM'][0]
            if type(pM_ratio) is not str:
                pM_ratio = pM_ratio.values[-1]
            pM_ratio_d = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'pM'][0]
            if type(pM_ratio_d) is not str:
                pM_ratio_d = pM_ratio_d.values[-1]
            info = {'daughter': pD, 'mother_d': pM_d, 'd': min(Ml.loc[:, 'd']),
                    'mother_ratio': pM_ratio, 'ratio': min(Ml.loc[:, 'ratio']),
                    'mother_ratio_d': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria'])}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': True}
            return [False], [info], [pM_d], [pM_ratio], [pM_ratio_d], [run_list]

    def denoising_parallel(self, pos):
        if pos == self.q1:
            print('10% aprox')
        elif pos == self.q2:
            print('20% aprox')
        elif pos == self.q3:
            print('30% aprox')
        elif pos == self.q4:
            print('40% aprox')
        elif pos == self.q5:
            print('50% aprox')
        elif pos == self.q6:
            print('60% aprox')
        elif pos == self.q7:
            print('70% aprox')
        elif pos == self.q8:
            print('80% aprox')
        elif pos == self.q9:
            print('90% aprox')
        elif pos == self.q10:
            print('100% aprox')

        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]
        position = len(self.run_list)
        # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.run_list[a].get('daughter'):
                continue
            pM = self.run_list[a].get('id')
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > self.max_ratio:
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d = lv.distance(pDseq, pMseq)
            # if d:
            # if dd exist & d is higher than dd:
            # break go to next pM
            if 'dd' in locals():
                if dd == 1:
                    break
                if d >= dd:
                    continue
                # if dd doesn't exist | d is smaller than dd:
                # dd = d
                else:
                    dd = d
            else:
                dd = d
                # if Edgar's equation:
            if b_ratio <= (1 / 2) ** (self.alpha * d + 1):
                # TRUE:
                # add to Ml:
                xavier_criteria = b_ratio / ((1 / 2) ** (self.alpha * d + 1))
                df1 = [
                    {'pM': pM, 'pMpos': pMpos, 'pD': pD, 'ratio': b_ratio, 'd': d, 'xavier_criteria': xavier_criteria}]
                Ml = Ml.append(df1)
                # identification of the pM
                # ratio
                # d
                # Xavier's criteria ---> ratio/((1/2)^(alpha*d + 1))
                # (the less value the best pM)
                # FALSE ---> break this pM
            else:
                continue
            # editing variables and files
        if Ml.empty:  # means that is not a daughter
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_ratio_d': None, 'xavier_criteria': None}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': False}
            return True, info, pD, pD, pD, run_list
            # exist
        else:  # it is a daughter
            # print pD name to each pM depending on different criteria
            # pDinfo = np.transpose(pd.DataFrame(self.data_initial.loc[pos, abund_col_names]))
            pM_d = Ml.loc[(Ml['d'] == min(Ml.loc[:, 'd'])), 'pM'][0]
            pM_ratio = Ml.loc[(Ml['ratio'] == min(Ml.loc[:, 'ratio'])), 'pM'][0]
            if type(pM_ratio) is not str:
                pM_ratio = pM_ratio.values[-1]
            pM_ratio_d = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'pM'][0]
            if type(pM_ratio_d) is not str:
                pM_ratio_d = pM_ratio_d.values[-1]
            info = {'daughter': pD, 'mother_d': pM_d, 'd': min(Ml.loc[:, 'd']),
                    'mother_ratio': pM_ratio, 'ratio': min(Ml.loc[:, 'ratio']),
                    'mother_ratio_d': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria'])}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': True}
            return False, info, pM_d, pM_ratio, pM_ratio_d, run_list

    def denoising_ratio(self, pos):
        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]
        position = len(self.run_list)
            # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.run_list[a].get('daughter'):
                continue
            pM = self.run_list[a].get('id')
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > self.max_ratio:
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d = lv.distance(pDseq, pMseq)
            # if d:
            # if dd exist & d is higher than dd:
            # break go to next pM
            if 'dd' in locals():
                if dd == 1:
                    break
                if d >= dd:
                    continue
                # if dd doesn't exist | d is smaller than dd:
                # dd = d
                else:
                    dd = d
            else:
                dd = d
                # if Edgar's equation:
            if b_ratio <= (1 / 2) ** (self.alpha * d + 1):
                # TRUE:
                # add to Ml:
                xavier_criteria = b_ratio / ((1 / 2) ** (self.alpha * d + 1))
                df1 = [
                    {'pM': pM, 'pMpos': pMpos, 'pD': pD, 'ratio': b_ratio, 'd': d, 'xavier_criteria': xavier_criteria}]
                Ml = Ml.append(df1)
                # identification of the pM
                # ratio
                # d
                # Xavier's criteria ---> ratio/((1/2)^(alpha*d + 1))
                # (the less value the best pM)
                # FALSE ---> break this pM
                break
            else:
                continue

            # editing variables and files
        if Ml.empty:  # means that is not a daughter
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_ratio_d': None, 'xavier_criteria': None}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': False}
            return [True], [info], [pD], [run_list]
        else:  # it is a daughter
            # print pD name to each pM depending on different criteria
            # _mothers_d
            pM_d = Ml.loc[(Ml['d'] == min(Ml.loc[:, 'd'])), 'pM'][0]
            pM_ratio = Ml.loc[(Ml['ratio'] == min(Ml.loc[:, 'ratio'])), 'pM'][0]
            if type(pM_ratio) is not str:
                pM_ratio = pM_ratio.values[-1]
            pM_ratio_d = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'pM'][0]
            if type(pM_ratio_d) is not str:
                pM_ratio_d = pM_ratio_d.values[-1]
            info = {'daughter': pD, 'mother_d': pM_d, 'd': min(Ml.loc[:, 'd']),
                    'mother_ratio': pM_ratio, 'ratio': min(Ml.loc[:, 'ratio']),
                    'mother_ratio_d': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria'])}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': True}
            return [False], [info], [pM_ratio], [run_list]

    def denoising_parallel_ratio(self, pos):

        if pos == self.q1:
            print('10% aprox')
        elif pos == self.q2:
            print('20% aprox')
        elif pos == self.q3:
            print('30% aprox')
        elif pos == self.q4:
            print('40% aprox')
        elif pos == self.q5:
            print('50% aprox')
        elif pos == self.q6:
            print('60% aprox')
        elif pos == self.q7:
            print('70% aprox')
        elif pos == self.q8:
            print('80% aprox')
        elif pos == self.q9:
            print('90% aprox')
        elif pos == self.q10:
            print('100% aprox')

        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]
        position = len(self.run_list)
        # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.run_list[a].get('daughter'):
                continue
            pM = self.run_list[a].get('id')
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > self.max_ratio:
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d = lv.distance(pDseq, pMseq)
            # if d:
            # if dd exist & d is higher than dd:
            # break go to next pM
            if 'dd' in locals():
                if dd == 1:
                    break
                if d >= dd:
                    continue
                # if dd doesn't exist | d is smaller than dd:
                # dd = d
                else:
                    dd = d
            else:
                dd = d
                # if Edgar's equation:
            if b_ratio <= (1 / 2) ** (self.alpha * d + 1):
                # TRUE:
                # add to Ml:
                xavier_criteria = b_ratio / ((1 / 2) ** (self.alpha * d + 1))
                df1 = [
                    {'pM': pM, 'pMpos': pMpos, 'pD': pD, 'ratio': b_ratio, 'd': d, 'xavier_criteria': xavier_criteria}]
                Ml = Ml.append(df1)
                # identification of the pM
                # ratio
                # d
                # Xavier's criteria ---> ratio/((1/2)^(alpha*d + 1))
                # (the less value the best pM)
                # FALSE ---> break this pM
                break
            else:
                continue
            # editing variables and files
        if Ml.empty:  # means that is not a daughter
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_ratio_d': None, 'xavier_criteria': None}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': False}
            return True, info, pD, run_list
            # exist
        else:  # it is a daughter
            # print pD name to each pM depending on different criteria
            # pDinfo = np.transpose(pd.DataFrame(self.data_initial.loc[pos, abund_col_names]))
            pM_d = Ml.loc[(Ml['d'] == min(Ml.loc[:, 'd'])), 'pM'][0]
            pM_ratio = Ml.loc[(Ml['ratio'] == min(Ml.loc[:, 'ratio'])), 'pM'][0]
            if type(pM_ratio) is not str:
                pM_ratio = pM_ratio.values[-1]
            pM_ratio_d = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'pM'][0]
            if type(pM_ratio_d) is not str:
                pM_ratio_d = pM_ratio_d.values[-1]
            info = {'daughter': pD, 'mother_d': pM_d, 'd': min(Ml.loc[:, 'd']),
                    'mother_ratio': pM_ratio, 'ratio': min(Ml.loc[:, 'ratio']),
                    'mother_ratio_d': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria'])}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': True}
            return False, info, pM_ratio, run_list

    def denoising_Adcorrected(self, pos):
        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]
        position = len(self.run_list)
        # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria',
                                   'difpos1', 'difpos2', 'difpos3'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.run_list[a].get('daughter'):
                continue
            pM = self.run_list[a].get('id')
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > self.max_ratio:
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d, difpos1, difpos2, difpos3 = difference(seq1=pDseq, seq2=pMseq,
                                                           initial_pos=self.initial_pos,
                                                           Ad1=self.Ad1, Ad2=self.Ad2, Ad3=self.Ad3)
            # if d:
            # if dd exist & d is higher than dd:
            # break go to next pM
            if 'dd' in locals():
                if dd == min(self.Ad1, self.Ad2, self.Ad3):
                    break
                if d >= dd:
                    continue
                # if dd doesn't exist | d is smaller than dd:
                # dd = d
                else:
                    dd = d
            else:
                dd = d
                # if Edgar's equation:
            if b_ratio <= (1 / 2) ** (self.alpha * d + 1):
                # TRUE:
                # add to Ml:
                xavier_criteria = b_ratio / ((1 / 2) ** (self.alpha * d + 1))
                df1 = [
                    {'pM': pM, 'pMpos': pMpos, 'pD': pD, 'ratio': b_ratio, 'd': d, 'xavier_criteria': xavier_criteria,
                     'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}]
                Ml = Ml.append(df1)
                # identification of the pM
                # ratio
                # d
                # Xavier's criteria ---> ratio/((1/2)^(alpha*d + 1))
                # (the less value the best pM)
                # FALSE ---> break this pM
            else:
                continue

            # editing variables and files
        if Ml.empty:  # means that is not a daughter
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_ratio_d': None, 'xavier_criteria': None,
                    'difpos1': None, 'difpos2': None, 'difpos3': None}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': False}
            return [True], [info], [pD], [pD], [pD], [run_list]
        else:  # it is a daughter
            # print pD name to each pM depending on different criteria
            # _mothers_d
            pM_d = Ml.loc[(Ml['d'] == min(Ml.loc[:, 'd'])), 'pM'][0]
            pM_ratio = Ml.loc[(Ml['ratio'] == min(Ml.loc[:, 'ratio'])), 'pM'][0]
            if type(pM_ratio) is not str:
                pM_ratio = pM_ratio.values[-1]
            pM_ratio_d = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'pM'][0]
            difpos1 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos1'][0]
            difpos2 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos2'][0]
            difpos3 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos3'][0]
            if type(pM_ratio_d) is not str:
                pM_ratio_d = pM_ratio_d.values[-1]
            info = {'daughter': pD, 'mother_d': pM_d, 'd': min(Ml.loc[:, 'd']),
                    'mother_ratio': pM_ratio, 'ratio': min(Ml.loc[:, 'ratio']),
                    'mother_ratio_d': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria']),
                    'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': True}
            return [False], [info], [pM_d], [pM_ratio], [pM_ratio_d], [run_list]

    def denoising_Adcorrected_parallel(self, pos):
        if pos == self.q1:
            print('10% aprox')
        elif pos == self.q2:
            print('20% aprox')
        elif pos == self.q3:
            print('30% aprox')
        elif pos == self.q4:
            print('40% aprox')
        elif pos == self.q5:
            print('50% aprox')
        elif pos == self.q6:
            print('60% aprox')
        elif pos == self.q7:
            print('70% aprox')
        elif pos == self.q8:
            print('80% aprox')
        elif pos == self.q9:
            print('90% aprox')
        elif pos == self.q10:
            print('100% aprox')

        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]
        position = len(self.run_list)
        # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria',
                                   'difpos1', 'difpos2', 'difpos3'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.run_list[a].get('daughter'):
                continue
            pM = self.run_list[a].get('id')
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > self.max_ratio:
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d, difpos1, difpos2, difpos3 = difference(seq1=pDseq, seq2=pMseq,
                                                           initial_pos=self.initial_pos,
                                                           Ad1=self.Ad1, Ad2=self.Ad2, Ad3=self.Ad3)
            # if d:
            # if dd exist & d is higher than dd:
            # break go to next pM
            if 'dd' in locals():
                if dd == min(self.Ad1, self.Ad2, self.Ad3):
                    break
                if d >= dd:
                    continue
                # if dd doesn't exist | d is smaller than dd:
                # dd = d
                else:
                    dd = d
            else:
                dd = d
                # if Edgar's equation:
            if b_ratio <= (1 / 2) ** (self.alpha * d + 1):
                # TRUE:
                # add to Ml:
                xavier_criteria = b_ratio / ((1 / 2) ** (self.alpha * d + 1))
                df1 = [
                    {'pM': pM, 'pMpos': pMpos, 'pD': pD, 'ratio': b_ratio, 'd': d, 'xavier_criteria': xavier_criteria,
                     'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}]
                Ml = Ml.append(df1)
                # identification of the pM
                # ratio
                # d
                # Xavier's criteria ---> ratio/((1/2)^(alpha*d + 1))
                # (the less value the best pM)
                # FALSE ---> break this pM
            else:
                continue
            # editing variables and files
        if Ml.empty:  # means that is not a daughter
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_ratio_d': None, 'xavier_criteria': None,
                    'difpos1': None, 'difpos2': None, 'difpos3': None}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': False}
            return True, info, pD, pD, pD, run_list
            # exist
        else:  # it is a daughter
            # print pD name to each pM depending on different criteria
            # pDinfo = np.transpose(pd.DataFrame(self.data_initial.loc[pos, abund_col_names]))
            pM_d = Ml.loc[(Ml['d'] == min(Ml.loc[:, 'd'])), 'pM'][0]
            pM_ratio = Ml.loc[(Ml['ratio'] == min(Ml.loc[:, 'ratio'])), 'pM'][0]
            if type(pM_ratio) is not str:
                pM_ratio = pM_ratio.values[-1]
            pM_ratio_d = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'pM'][0]
            difpos1 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos1'][0]
            difpos2 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos2'][0]
            difpos3 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos3'][0]
            if type(pM_ratio_d) is not str:
                pM_ratio_d = pM_ratio_d.values[-1]
            info = {'daughter': pD, 'mother_d': pM_d, 'd': min(Ml.loc[:, 'd']),
                    'mother_ratio': pM_ratio, 'ratio': min(Ml.loc[:, 'ratio']),
                    'mother_ratio_d': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria']),
                    'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': True}
            return False, info, pM_d, pM_ratio, pM_ratio_d, run_list

    def denoising_Adcorrected_ratio(self, pos):
        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]
        position = len(self.run_list)
        # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria',
                                   'difpos1', 'difpos2', 'difpos3'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.run_list[a].get('daughter'):
                continue
            pM = self.run_list[a].get('id')
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > self.max_ratio:
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d, difpos1, difpos2, difpos3 = difference(seq1=pDseq, seq2=pMseq,
                                                           initial_pos=self.initial_pos,
                                                           Ad1=self.Ad1, Ad2=self.Ad2, Ad3=self.Ad3)
            # if d:
            # if dd exist & d is higher than dd:
            # break go to next pM
            if 'dd' in locals():
                if dd == min(self.Ad1, self.Ad2, self.Ad3):
                    break
                if d >= dd:
                    continue
                # if dd doesn't exist | d is smaller than dd:
                # dd = d
                else:
                    dd = d
            else:
                dd = d
                # if Edgar's equation:
            if b_ratio <= (1 / 2) ** (self.alpha * d + 1):
                # TRUE:
                # add to Ml:
                xavier_criteria = b_ratio / ((1 / 2) ** (self.alpha * d + 1))
                df1 = [
                    {'pM': pM, 'pMpos': pMpos, 'pD': pD, 'ratio': b_ratio, 'd': d, 'xavier_criteria': xavier_criteria,
                     'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}]
                Ml = Ml.append(df1)
                # identification of the pM
                # ratio
                # d
                # Xavier's criteria ---> ratio/((1/2)^(alpha*d + 1))
                # (the less value the best pM)
                # FALSE ---> break this pM
                break
            else:
                continue

            # editing variables and files
        if Ml.empty:  # means that is not a daughter
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_ratio_d': None, 'xavier_criteria': None,
                    'difpos1': None, 'difpos2': None, 'difpos3': None}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': False}
            return [True], [info], [pD], [run_list]
        else:  # it is a daughter
            # print pD name to each pM depending on different criteria
            # _mothers_d
            pM_d = Ml.loc[(Ml['d'] == min(Ml.loc[:, 'd'])), 'pM'][0]
            pM_ratio = Ml.loc[(Ml['ratio'] == min(Ml.loc[:, 'ratio'])), 'pM'][0]
            if type(pM_ratio) is not str:
                pM_ratio = pM_ratio.values[-1]
            pM_ratio_d = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'pM'][0]
            difpos1 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos1'][0]
            difpos2 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos2'][0]
            difpos3 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos3'][0]
            if type(pM_ratio_d) is not str:
                pM_ratio_d = pM_ratio_d.values[-1]
            info = {'daughter': pD, 'mother_d': pM_d, 'd': min(Ml.loc[:, 'd']),
                    'mother_ratio': pM_ratio, 'ratio': min(Ml.loc[:, 'ratio']),
                    'mother_ratio_d': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria']),
                    'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': True}
            return [False], [info], [pM_ratio], [run_list]

    def denoising_Adcorrected_parallel_ratio(self, pos):
        if pos == self.q1:
            print('10% aprox')
        elif pos == self.q2:
            print('20% aprox')
        elif pos == self.q3:
            print('30% aprox')
        elif pos == self.q4:
            print('40% aprox')
        elif pos == self.q5:
            print('50% aprox')
        elif pos == self.q6:
            print('60% aprox')
        elif pos == self.q7:
            print('70% aprox')
        elif pos == self.q8:
            print('80% aprox')
        elif pos == self.q9:
            print('90% aprox')
        elif pos == self.q10:
            print('100% aprox')

        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]
        position = len(self.run_list)
        # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria',
                                   'difpos1', 'difpos2', 'difpos3'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.run_list[a].get('daughter'):
                continue
            pM = self.run_list[a].get('id')
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > self.max_ratio:
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d, difpos1, difpos2, difpos3 = difference(seq1=pDseq, seq2=pMseq,
                                                           initial_pos=self.initial_pos,
                                                           Ad1=self.Ad1, Ad2=self.Ad2, Ad3=self.Ad3)
            # if d:
            # if dd exist & d is higher than dd:
            # break go to next pM
            if 'dd' in locals():
                if dd == min(self.Ad1, self.Ad2, self.Ad3):
                    break
                if d >= dd:
                    continue
                # if dd doesn't exist | d is smaller than dd:
                # dd = d
                else:
                    dd = d
            else:
                dd = d
                # if Edgar's equation:
            if b_ratio <= (1 / 2) ** (self.alpha * d + 1):
                # TRUE:
                # add to Ml:
                xavier_criteria = b_ratio / ((1 / 2) ** (self.alpha * d + 1))
                df1 = [
                    {'pM': pM, 'pMpos': pMpos, 'pD': pD, 'ratio': b_ratio, 'd': d, 'xavier_criteria': xavier_criteria,
                     'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}]
                Ml = Ml.append(df1)
                # identification of the pM
                # ratio
                # d
                # Xavier's criteria ---> ratio/((1/2)^(alpha*d + 1))
                # (the less value the best pM)
                # FALSE ---> break this pM
                break
            else:
                continue
            # editing variables and files
        if Ml.empty:  # means that is not a daughter
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_ratio_d': None, 'xavier_criteria': None,
                    'difpos1': None, 'difpos2': None, 'difpos3': None}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': False}
            return True, info, pD, run_list
            # exist
        else:  # it is a daughter
            # print pD name to each pM depending on different criteria
            # pDinfo = np.transpose(pd.DataFrame(self.data_initial.loc[pos, abund_col_names]))
            pM_d = Ml.loc[(Ml['d'] == min(Ml.loc[:, 'd'])), 'pM'][0]
            pM_ratio = Ml.loc[(Ml['ratio'] == min(Ml.loc[:, 'ratio'])), 'pM'][0]
            if type(pM_ratio) is not str:
                pM_ratio = pM_ratio.values[-1]
            pM_ratio_d = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'pM'][0]
            difpos1 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos1'][0]
            difpos2 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos2'][0]
            difpos3 = Ml.loc[(Ml['xavier_criteria'] == min(Ml.loc[:, 'xavier_criteria'])), 'difpos3'][0]
            if type(pM_ratio_d) is not str:
                pM_ratio_d = pM_ratio_d.values[-1]
            info = {'daughter': pD, 'mother_d': pM_d, 'd': min(Ml.loc[:, 'd']),
                    'mother_ratio': pM_ratio, 'ratio': min(Ml.loc[:, 'ratio']),
                    'mother_ratio_d': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria']),
                    'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}
            run_list = {'id': pD, self.count: pDabund, 'run': True, 'daughter': True}
            return False, info, pM_ratio, run_list

    def write_d_from_info(self, mother):
        row = [
            self.good_mothers[list(self.good_mothers.id == mother)][self.first_col_names].values.tolist()[0] +
            list(self.data_initial.loc[[[mother] +
                                        list(self.merge_data.daughter[self.merge_data.mother_d == mother])][0],
                                       self.abund_col_names].sum(0)) +
            self.good_mothers[list(self.good_mothers.id == mother)][self.seq].values.tolist()]

        return row

    def write_ratio_from_info(self, mother):
        row = [
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.first_col_names].values.tolist()[0] +
            list(self.data_initial.loc[[[mother] +
                                        list(self.merge_data.daughter[self.merge_data.mother_ratio == mother])][0],
                                       self.abund_col_names].sum(0)) +
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.seq].values.tolist()]

        return row

    def write_ratio_d_from_info(self, mother):
        row = [
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.first_col_names].values.tolist()[0] +
            list(self.data_initial.loc[[[mother] +
                                        list(self.merge_data.daughter[self.merge_data.mother_ratio_d == mother])][0],
                                       self.abund_col_names].sum(0)) +
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.seq].values.tolist()]

        return row

    def write_output_ratio(self, mother):
        row = [
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.first_col_names].values.tolist()[0] +
            list(self.data_initial.loc[list(pd.Series(self.denoised_ratio_output) == mother), self.abund_col_names].sum(0)) +
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.seq].values.tolist()]

        return row

    def write_output_d(self, mother):
        row = [
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.first_col_names].values.tolist()[0] +
            list(self.data_initial.loc[list(pd.Series(self.denoised_d_output) == mother), self.abund_col_names].sum(0)) +
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.seq].values.tolist()]

        return row

    def write_output_ratio_d(self, mother):
        row = [
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.first_col_names].values.tolist()[0] +
            list(self.data_initial.loc[list(pd.Series(self.denoised_ratio_d_output) == mother), self.abund_col_names].sum(0)) +
            self.good_mothers[list(self.good_mothers.loc[:, 'id'] == mother)][self.seq].values.tolist()]

        return row


def run_dnoise_testing(declass):
    if declass.cores > 1:
        if declass.output_type == 'ratio':
            while len(declass.run_list) < declass.data_initial.shape[0]:
                declass.min_mother = declass.run_list[-1].get(declass.count) * declass.max_ratio
                run_to = sum(declass.data_initial.loc[:, declass.count] > declass.min_mother)
                if run_to == len(declass.run_list):
                    run_to += 1
                declass.quartiles_run(declass.data_initial)
                print('running until %s reads' % declass.min_mother)
                print(len(declass.run_list) / declass.data_initial.shape[0] * 100, '%')
                pool = mp.Pool(declass.cores)
                [declass.good_seq_2,
                 declass.output_info_2,
                 declass.denoised_ratio_output_2,
                 declass.run_list_2] = zip(*pool.map(declass.denoising_parallel_ratio,
                                                [pos for pos in range(len(declass.run_list), run_to)]))
                pool.close()
                del pool
                declass.good_seq.extend(list(declass.good_seq_2))
                declass.output_info.extend(list(declass.output_info_2))
                declass.denoised_d_output.extend(list(declass.denoised_d_output_2))
                declass.denoised_ratio_output.extend(list(declass.denoised_ratio_output_2))
                declass.denoised_ratio_d_output.extend(list(declass.denoised_ratio_d_output_2))
                declass.run_list.extend(list(declass.run_list_2))
                print('run until %s reads' % declass.min_mother)
                print(len(declass.run_list) / declass.data_initial.shape[0] * 100, '%')
                del (declass.good_seq_2, declass.output_info_2, declass.denoised_ratio_output_2, declass.run_list_2)
        else:
            while len(declass.run_list) < declass.data_initial.shape[0]:
                declass.min_mother = declass.run_list[-1].get(declass.count) * declass.max_ratio
                run_to = sum(declass.data_initial.loc[:, declass.count] > declass.min_mother)
                if run_to == len(declass.run_list):
                    run_to += 1
                declass.quartiles_run(declass.data_initial)
                print('running until %s reads' % declass.min_mother)
                print(len(declass.run_list) / declass.data_initial.shape[0] * 100, '%')
                pool = mp.Pool(declass.cores)
                [declass.good_seq_2,
                 declass.output_info_2,
                 declass.denoised_d_output_2,
                 declass.denoised_ratio_output_2,
                 declass.denoised_ratio_d_output_2,
                 declass.run_list_2] = zip(*pool.map(declass.denoising_parallel,
                                                [pos for pos in range(len(declass.run_list), run_to)]))
                pool.close()
                del pool
                declass.good_seq.extend(list(declass.good_seq_2)[:])
                declass.output_info.extend(list(declass.output_info_2)[:])
                declass.denoised_d_output.extend(list(declass.denoised_d_output_2)[:])
                declass.denoised_ratio_output.extend(list(declass.denoised_ratio_output_2)[:])
                declass.denoised_ratio_d_output.extend(list(declass.denoised_ratio_d_output_2)[:])
                declass.run_list.extend(list(declass.run_list_2)[:])
                print('run until %s reads' % declass.min_mother)
                print(len(declass.run_list) / declass.data_initial.shape[0] * 100, '%')
                del (declass.good_seq_2, declass.output_info_2, declass.denoised_d_output_2, declass.denoised_ratio_output_2,
                     declass.denoised_ratio_d_output_2, declass.run_list_2)
        if 'declass.min_mother' in locals():
            del declass.min_mother
        del (declass.cores, declass.alpha)
    else:
        if declass.output_type == 'ratio':
            for pos in tqdm(range(1, declass.data_initial.shape[0])):
                [declass.good_seq[len(declass.good_seq):],
                 declass.output_info[len(declass.output_info):],
                 declass.denoised_ratio_output[len(declass.denoised_ratio_output):],
                 declass.run_list[len(declass.run_list):]] = declass.denoising_ratio(pos)
        else:
            for pos in tqdm(range(1, declass.data_initial.shape[0])):
                [declass.good_seq[len(declass.good_seq):],
                 declass.output_info[len(declass.output_info):],
                 declass.denoised_d_output[len(declass.denoised_d_output):],
                 declass.denoised_ratio_output[len(declass.denoised_ratio_output):],
                 declass.denoised_ratio_d_output[len(declass.denoised_ratio_d_output):],
                 declass.run_list[len(declass.run_list):]] = declass.denoising(pos)


def difference(self, seq1, seq2, initial_pos, Ad1, Ad2, Ad3):
    difcount = 0
    difpos1 = 0
    difpos2 = 0
    difpos3 = 0
    for i, o in enumerate(seq1):
        if o != seq2[i]:
            dif_position = (i + initial_pos) % 3
            if dif_position == 1:
                difcount += Ad1
                difpos1 += 1
            elif dif_position == 2:
                difcount += Ad2
                difpos2 += 1
            elif dif_position == 0:
                difcount += Ad3
                difpos3 += 1
    return difcount, difpos1, difpos2, difpos3


def mother_id(self, test, M):
    return test == M


def modal_length(array):
    pdarray = pd.Series(array, dtype="category")
    uniques = list(pdarray.unique())
    counts = pd.Series(list(map(array.count, uniques)))
    most = max(counts)
    return list(itertools.compress(uniques, [list(counts == most)][0]))


class DeSub:
    def __init__(self):
        print("starting to denoise a new seq length subset")


def copy_to_subset(declass, desub, seq_length, len_seq):
    print("starting to denoise a new seq length subset")

    desub.alpha = declass.alpha
    desub.data_initial = declass.data_initial.loc[(np.asarray(seq_length) == len_seq)]
    desub.data_initial.index = list(range(desub.data_initial.shape[0]))
    desub.abund_col_names = declass.abund_col_names[:]
    desub.seq = declass.seq
    desub.count = declass.count
    desub.output_type = declass.output_type
    desub.cores = declass.cores
    desub.run_list = declass.run_list
    desub.max_ratio = declass.max_ratio
    desub.initial_pos = declass.initial_pos
    desub.Ad1 = declass.Ad1
    desub.Ad2 = declass.Ad2
    desub.Ad3 = declass.Ad3



