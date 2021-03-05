
import pandas as pd
import Levenshtein as lv
import json
from ast import literal_eval
import os
import getopt
import sys


class denoise_functions:
    not_possibleMothers = 0
    possibleMothers = 0
    data_initial = pd.DataFrame()
    runned_list = pd.DataFrame()
    alpha = 5
    min_mother = 20
    Ad1 = 0.8379801130824722
    Ad2 = 0.357379606161045
    Ad3 = 0.8379801130824722
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
    count = 'count'
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
    fasta_output = True
    fasta = True
    sep = '\t'
    end = 1
    start = 1

    def __init__(self):
        print("starting to denoise")

    def read_parameters(self, argument_list):
        short_options = "hi:o:P:f:F:c:s:z:n:a:q:p:e:y:"
        long_options = ["help", "input=", "output=", "part=", "fasta_input=", "fasta_output=", "cores=",
                        "start_sample_cols=",
                        "end_sample_cols=", "count_name=", "alpha=", "sequence=", "separation=", "entropy=",
                        "entropy_influence="]
        try:
            arguments, values = getopt.getopt(argument_list, short_options, long_options)
        except getopt.error as err:
            # Output error, and return with an error code
            print(str(err))
            sys.exit(2)
        for current_argument, current_value in arguments:
            if current_argument in ("-h", "--help"):
                print("Displaying help\n"
                      " -h --help Display help\n"
                      " -i --input input file path\n"
                      " -o --output common output files path\n"
                      " -P --part Denoise can be runned per parts, part 1 run just possible mothers and part 2 run the others\n"
                      "                     In part 1 a directory as database is named as --output\n"
                      "                     In part 2 (parallellitzable) this directory is removed and outputs are printed\n"
                      "                         Part 2 requires --input, --output and --cores if necessary\n"
                      " -f --fasta_input logical, if T (default), fasta file as input, if F .csv as input\n"
                      " -F --fasta_output logical, if T (default), fasta file as input, if F .csv as input\n"
                      " -c --cores number of cores, 1 by default\n"
                      " -s --start_sample_cols first sample column (1 == 1st col) if not given, just total read count expected\n"
                      " -z --end_sample_cols first sample column (1 == 1st col) if not given, just total read count expected\n"
                      " -n --count_name count name column (count/reads/size..) 'count' by default\n"
                      " -a --alpha alpha value, 5 by default\n"
                      " -q --sequence sequence column name, 'sequence' by default\n"
                      " -p --sep separation 1='\t'\n"
                      "                     2=','\n"
                      "                     3=';'\n"
                      " -e --entropy entropy of different positions [0.4298,0.1833,0.9256] by default\n"
                      " -y --entropy_influence logical, if T, Ad correction parameter is used. When this happens, only "
                      "sequences with mode length are computed")
                exit()
            elif current_argument in ("-i", "--input"):
                print("Denoising %s file" % current_value)
                self.MOTUfile = current_value
                arg_i = True
            elif current_argument in ("-o", "--output"):
                print("Output files will be named %s*" % current_value)
                self.MOTUoutfile = current_value
                arg_o = True
            elif current_argument in ("-P", "--Part"):
                print("Part %s running" % current_value)
                self.part = int(current_value)
                self.parted = True
                arg_P = True
            elif current_argument in ("-c", "--cores"):
                print("Running with %s cores" % current_value)
                self.cores = int(current_value)
                arg_c = True
            elif current_argument in ("-f", "--fasta_input"):
                if current_value == "T":
                    self.fasta = True
                else:
                    self.fasta = False
                arg_f = True
                print("Fasta file: %s" % self.fasta)
            elif current_argument in ("-F", "--fasta_output"):
                if current_value == "T":
                    self.fasta_output = True
                else:
                    self.fasta_output = False
                arg_F = True
                print("Fasta output file: %s" % self.fasta_output)
            elif current_argument in ("-s", "--start_sample_cols"):
                self.start = int(current_value)
                arg_s = True
                print("Abundant cols starts in: %s" % current_value)
            elif current_argument in ("-z", "--end_sample_cols"):
                self.end = int(current_value)
                arg_z = True
                print("Abundant cols ends in: %s" % current_value)
            elif current_argument in ("-n", "--count_name"):
                self.count = current_value
                arg_n = True
                print("Count col name: %s" % current_value)
            elif current_argument in ("-a", "--alpha"):
                self.alpha = int(current_value)
                arg_a = True
                print("Alpha: %s" % current_value)
            elif current_argument in ("-q", "--sequence"):
                self.seq = current_value
                arg_q = True
                print("Seq: %s" % current_value)
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
            elif current_argument in ("-e", "--entropy"):
                E1, E2, E3 = current_value.split(",")
                print(str("E1: " + E1 + " \nE2: " + E2 + " \nE3: " + E3))
                E1 = float(E1)
                E2 = float(E2)
                E3 = float(E3)
                arg_e = True
            elif current_argument in ("-y", "--entropy_influence"):
                if current_value == "T":
                    self.entropy = True
                    arg_y = True
                    print("Is entropy taken into account: %s" % self.entropy)
                else:
                    print("Is entropy taken into account: False")

        if 'arg_i' not in locals():
            print("Err: input file needed")
            exit()
        if 'arg_o' not in locals():
            print("Err: output path needed")
            exit()
        if 'arg_P' not in locals():
            print("Running all at once")
            self.parted = False
            self.part = 1
        if self.part == 1:
            if 'arg_c' not in locals():
                print("cores not given, 1 core by default")
                self.cores = 1
            if 'arg_f' not in locals():
                # de moment per a la opcio amb fasta no hi ha la opció de posar samples
                print("by default, fasta file expected")
                self.fasta = True
                self.count = "count"
                self.justcount = True
                if 'arg_F' not in locals():
                    self.fasta_output = True
                    print("by default, fasta output")
            else:
                if "arg_n" not in locals():
                    print("count_name not given, 'count' by default")
                    self.count = "count"
                if "arg_s" not in locals():
                    if "arg_z" not in locals():
                        print("start and end not given, no samples in file")
                        abund_col_names = [self.count]
                        self.justcount = True
                    else:
                        print("Err: end given but no start")
                        exit()
                elif "arg_z" not in locals():
                    print("Err: start given but no end")
                    exit()
                else:
                    self.justcount = False
                    if self.fasta_output:
                        print('WARNING!!!! .csv output when samples information is given')
                    self.fasta_output = False
                    print("output as .csv")
                if "arg_q" not in locals():
                    print("sequence not given, 'sequence' by default")
                    self.seq = 'sequence'
                if "arg_p" not in locals():
                    print("Separation not given, '\t' by default")
                    self.sep = '\t'
            if "arg_a" not in locals():
                print("alpha not given, 5 by default")
                self.alpha = 5
            if "arg_y" not in locals():
                self.entropy = False
                print("Ad correction not applied")
                self.Ad1, self.Ad2, self.Ad3 = (1, 1, 1)
            else:
                if "arg_e" not in locals():
                    E1, E2, E3 = (0.4298, 0.1833, 0.9256)
                    print("Entropy set as 0.4298, 0.1833 and 0.9256 by default")
                # defining Ad correction factor taking into account Entropy
                self.Ad1 = E1 / (E1 + E2 + E3) * 3
                self.Ad2 = E2 / (E1 + E2 + E3) * 3
                self.Ad3 = E3 / (E1 + E2 + E3) * 3

    def quartiles_runned(self):
        q = self.not_possibleMothers - 1 - self.possibleMothers

        self.q1 = int(q * 1 / 10) + self.possibleMothers
        self.q2 = int(q * 2 / 10) + self.possibleMothers
        self.q3 = int(q * 3 / 10) + self.possibleMothers
        self.q4 = int(q * 4 / 10) + self.possibleMothers
        self.q5 = int(q * 5 / 10) + self.possibleMothers
        self.q6 = int(q * 6 / 10) + self.possibleMothers
        self.q7 = int(q * 7 / 10) + self.possibleMothers
        self.q8 = int(q * 8 / 10) + self.possibleMothers
        self.q9 = int(q * 9 / 10) + self.possibleMothers
        self.q10 = self.not_possibleMothers - 1

    def denoising(self, pos):
        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]
        if pos == 0:
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_xavier_criteria': None, 'xavier_criteria': None}
            self.runned_list.loc[pos, 'runned'] = True
            # the return: good_seq / executed / info / denoised_d_output / denoised_ratio_output / denoised_ratio_d_output
            return [True], [info], [pD], [pD], [pD]
        if pD in list(self.runned_list.loc[:, 'id']):
            position = self.runned_list[self.runned_list['id'] == pD].index.tolist()[0]
        else:
            position = self.runned_list.shape[0]
            # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.runned_list.loc[a].at['daughter']:
                continue
            pM = self.runned_list.iloc[a, 0]
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > (1 / self.min_mother):
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
                    'mother_xavier_criteria': None, 'xavier_criteria': None}
            self.runned_list.loc[pos, 'runned'] = True
            return [True], [info], [pD], [pD], [pD]
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
                    'mother_xavier_criteria': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria'])}
            self.runned_list.loc[pos, 'daughter'] = True
            self.runned_list.loc[pos, 'runned'] = True
            return [False], [info], [pM_d], [pM_ratio], [pM_ratio_d]

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
        position = self.runned_list.shape[0]
        # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.runned_list.iloc[a].at['daughter']:
                continue
            pM = self.runned_list.iloc[a, 0]
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > (1 / self.min_mother):
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
                    'mother_xavier_criteria': None, 'xavier_criteria': None}
            return True, info, pD, pD, pD
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
                    'mother_xavier_criteria': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria'])}
            return False, info, pM_d, pM_ratio, pM_ratio_d

    def difference(self, seq1, seq2):
        # està adaptat per al nostre fragment el qual comença en la posició 2
        difcount = 0
        difpos1 = 0
        difpos2 = 0
        difpos3 = 0
        for i, o in enumerate(seq1):
            if o != seq2[i]:
                dif_position = (i + self.initial_pos) % 3
                if dif_position == 1:
                    difcount += self.Ad1
                    difpos1 += 1
                elif dif_position == 2:
                    difcount += self.Ad2
                    difpos2 += 1
                elif dif_position == 0:
                    difcount += self.Ad3
                    difpos3 += 1
        return difcount, difpos1, difpos2, difpos3

    def denoising_Adcorrected(self, pos):
        pD = self.data_initial.loc[pos, 'id']
        pDseq = self.data_initial.loc[pos, self.seq]
        pDabund = self.data_initial.loc[pos, self.count]

        if pos == 0:
            info = {'daughter': pD, 'mother_d': None, 'd': None,
                    'mother_ratio': None, 'ratio': None,
                    'mother_xavier_criteria': None, 'xavier_criteria': None,
                    'difpos1': None, 'difpos2': None, 'difpos3': None}]
            self.runned_list.loc[pos, 'runned'] = True
            # the return: good_seq / executed / info / denoised_d_output / denoised_ratio_output / denoised_ratio_d_output
            return [True], [info], [pD], [pD], [pD]
        if pD in list(self.runned_list.loc[:, 'id']):
            position = self.runned_list[self.runned_list['id'] == pD].index.tolist()[0]
        else:
            position = self.runned_list.shape[0]
            # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria',
                     'difpos1', 'difpos2', 'difpos3'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.runned_list.loc[a].at['daughter']:
                continue
            pM = self.runned_list.iloc[a, 0]
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > (1 / self.min_mother):
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d, difpos1, difpos2, difpos3 = self.difference(seq1=pDseq, seq2=pMseq)
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
                     'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}
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
                    'mother_xavier_criteria': None, 'xavier_criteria': None,
                    'difpos1': None, 'difpos2': None, 'difpos3': None}
            self.runned_list.loc[pos, 'runned'] = True
            return [True], [info], [pD], [pD], [pD]
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
                    'mother_xavier_criteria': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria']),
                    'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}
            self.runned_list.loc[pos, 'daughter'] = True
            self.runned_list.loc[pos, 'runned'] = True
            return [False], [info], [pM_d], [pM_ratio], [pM_ratio_d]

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
        position = self.runned_list.shape[0]
        # create void list for pM info ---> Motherslist (Ml)
        Ml = pd.DataFrame(columns=['pM', 'pMpos', 'pD', 'ratio', 'd', 'xavier_criteria',
                     'difpos1', 'difpos2', 'difpos3'])
        # compare with each bigger seq possible Mother (pM).
        for a in range(position):
            # if the pM is running wait (this should be done with a while loop)
            # if the pM is daughter break to the next pM
            if self.runned_list.iloc[a].at['daughter']:
                continue
            pM = self.runned_list.iloc[a, 0]
            pMpos = a
            pMseq = self.data_initial.loc[pMpos, self.seq]
            pMabund = self.data_initial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > (1 / self.min_mother):
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d, difpos1, difpos2, difpos3 = self.difference(seq1=pDseq, seq2=pMseq)
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
                    'mother_xavier_criteria': None, 'xavier_criteria': None,
                    'difpos1': None, 'difpos2': None, 'difpos3': None}
            return True, info, pD, pD, pD
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
                    'mother_xavier_criteria': pM_ratio_d, 'xavier_criteria': min(Ml.loc[:, 'xavier_criteria']),
                    'difpos1': difpos1, 'difpos2': difpos2, 'difpos3': difpos3}
            return False, info, pM_d, pM_ratio, pM_ratio_d

    def mother_id(self, test, M):
        return test == M

    def write_variables(self):
        variables = {"entropy": self.entropy,
                     "cores": self.cores,
                     "alpha": self.alpha,
                     "min_mother": self.min_mother,
                     "possibleMothers": self.possibleMothers,
                     "not_possibleMothers": self.not_possibleMothers,
                     "MOTUfile": self.MOTUfile,
                     "MOTUoutfile": self.MOTUoutfile,
                     "justcount": self.justcount,
                     "abund_col_names": self.abund_col_names,
                     "first_col_names": self.first_col_names,
                     "seq": self.seq,
                     "count": self.count,
                     "fasta_output": self.fasta_output}

        if self.entropy:
            variables['Ad1'] = self.Ad1
            variables['Ad2'] = self.Ad2
            variables['Ad3'] = self.Ad3

        directory = str(self.MOTUoutfile + '_dir')

        print(str('data_base in ' + directory))

        if not os.path.isdir(directory):
            os.mkdir(directory)

        with open(str(directory + '/variables.json'), 'w') as doc_json:
            json.dump(str(variables), doc_json)

        with open(str(directory + '/data_initial.json'), 'w') as doc_json:
            self.data_initial.to_json(doc_json, orient='table')

        with open(str(directory + '/runned_list.json'), 'w') as doc_json:
            self.runned_list.to_json(doc_json, orient='table')

        with open(str(directory + '/good_seq.json'), 'w') as doc_json:
            json.dump(str(self.good_seq), doc_json)

        with open(str(directory + '/output_info.json'), 'w') as doc_json:
            json.dump(str(self.output_info), doc_json)

        with open(str(directory + '/denoised_d.json'), 'w') as doc_json:
            json.dump(str(self.denoised_d_output), doc_json)

        with open(str(directory + '/denoised_ratio.json'), 'w') as doc_json:
            json.dump(str(self.denoised_ratio_output), doc_json)

        with open(str(directory + '/denoised_ratio_d_output.json'), 'w') as doc_json:
            json.dump(str(self.denoised_ratio_d_output), doc_json)

    def read_variables(self):
        directory = str(self.MOTUoutfile + '_dir')

        with open(str(directory + '/variables.json'), 'r') as doc_json:
            variables = literal_eval(json.load(doc_json))

        with open(str(directory + '/data_initial.json'), 'r') as doc_json:
            self.data_initial = pd.read_json(doc_json, orient='table')

        with open(str(directory + '/runned_list.json'), 'r') as doc_json:
            self.runned_list = pd.read_json(doc_json, orient='table')

        with open(str(directory + '/good_seq.json'), 'r') as doc_json:
            self.good_seq = literal_eval(json.load(doc_json))

        with open(str(directory + '/output_info.json'), 'r') as doc_json:
            self.output_info = literal_eval(json.load(doc_json))

        with open(str(directory + '/denoised_d.json'), 'r') as doc_json:
            self.denoised_d_output = literal_eval(json.load(doc_json))

        with open(str(directory + '/denoised_ratio.json'), 'r') as doc_json:
            self.denoised_ratio_output = literal_eval(json.load(doc_json))

        with open(str(directory + '/denoised_ratio_d_output.json'), 'r') as doc_json:
            self.denoised_ratio_d_output = literal_eval(json.load(doc_json))

        self.entropy = variables['entropy']
        self.cores = variables['cores']
        self.alpha = variables['alpha']
        self.min_mother = variables['min_mother']
        self.possibleMothers = variables['possibleMothers']
        self.not_possibleMothers = variables['not_possibleMothers']
        self.justcount = variables['justcount']
        self.abund_col_names = variables['abund_col_names']
        self.first_col_names = variables['first_col_names']
        self.seq = variables['seq']
        self.count = variables['count']
        self.fasta_output = variables['fasta_output']

        if variables['entropy']:
            self.Ad1 = variables['Ad1']
            self.Ad2 = variables['Ad2']
            self.Ad3 = variables['Ad3']

        # shutil.rmtree(directory, ignore_errors=True)

    def read_variables2(self):
        directory = str(self.MOTUoutfile + '_dir')
        # part 1
        with open(str(directory + '/variables.json'), 'r') as doc_json:
            variables = literal_eval(json.load(doc_json))

        with open(str(directory + '/data_initial.json'), 'r') as doc_json:
            self.data_initial = pd.read_json(doc_json, orient='table')

        with open(str(directory + '/runned_list.json'), 'r') as doc_json:
            self.runned_list = pd.read_json(doc_json, orient='table')

        with open(str(directory + '/good_seq.json'), 'r') as doc_json:
            self.good_seq = literal_eval(json.load(doc_json))

        with open(str(directory + '/output_info.json'), 'r') as doc_json:
            self.output_info = literal_eval(json.load(doc_json))

        with open(str(directory + '/denoised_d.json'), 'r') as doc_json:
            self.denoised_d_output = literal_eval(json.load(doc_json))

        with open(str(directory + '/denoised_ratio.json'), 'r') as doc_json:
            self.denoised_ratio_output = literal_eval(json.load(doc_json))

        with open(str(directory + '/denoised_ratio_d_output.json'), 'r') as doc_json:
            self.denoised_ratio_d_output = literal_eval(json.load(doc_json))
        # part 2
        with open(str(directory + '/good_seq_daughters.json'), 'r') as doc_json:
            self.good_seq_2 = literal_eval(json.load(doc_json))

        with open(str(directory + '/output_info_daughters.json'), 'r') as doc_json:
            self.output_info_2 = literal_eval(json.load(doc_json))

        with open(str(directory + '/denoised_d_daughters.json'), 'r') as doc_json:
            self.denoised_d_output_2 = literal_eval(json.load(doc_json))

        with open(str(directory + '/denoised_ratio_daughters.json'), 'r') as doc_json:
            self.denoised_ratio_output_2 = literal_eval(json.load(doc_json))

        with open(str(directory + '/denoised_ratio_d_output_daughters.json'), 'r') as doc_json:
            self.denoised_ratio_d_output_2 = literal_eval(json.load(doc_json))

        self.entropy = variables['entropy']
        self.cores = variables['cores']
        self.alpha = variables['alpha']
        self.min_mother = variables['min_mother']
        self.possibleMothers = variables['possibleMothers']
        self.not_possibleMothers = variables['not_possibleMothers']
        self.justcount = variables['justcount']
        self.abund_col_names = variables['abund_col_names']
        self.first_col_names = variables['first_col_names']
        self.seq = variables['seq']
        self.count = variables['count']
        self.fasta_output = variables['fasta_output']

        if variables['entropy']:
            self.Ad1 = variables['Ad1']
            self.Ad2 = variables['Ad2']
            self.Ad3 = variables['Ad3']

        # shutil.rmtree(directory, ignore_errors=True)

    def write_outputs(self):

        directory = str(self.MOTUoutfile + '_dir')

        print(str('outputs in ' + directory))

        if not os.path.isdir(directory):
            os.mkdir(directory)

        if not os.path.isfile(str(directory + '/data_initial.json')):
            with open(str(directory + '/data_initial.json'), 'w') as doc_json:
                self.data_initial.to_json(doc_json, orient='table')
        if not os.path.isfile(str(directory + '/runned_list.json')):
            with open(str(directory + '/runned_list.json'), 'w') as doc_json:
                self.runned_list.to_json(doc_json, orient='table')
        if not os.path.isfile(str(directory + '/good_seq.json')):
            with open(str(directory + '/good_seq.json'), 'w') as doc_json:
                json.dump(str(self.good_seq), doc_json)
        if not os.path.isfile(str(directory + '/output_info.json')):
            with open(str(directory + '/output_info.json'), 'w') as doc_json:
                json.dump(str(self.output_info), doc_json)
        if not os.path.isfile(str(directory + '/denoised_d.json')):
            with open(str(directory + '/denoised_d.json'), 'w') as doc_json:
                json.dump(str(self.denoised_d_output), doc_json)
        if not os.path.isfile(str(directory + '/denoised_ratio.json')):
            with open(str(directory + '/denoised_ratio.json'), 'w') as doc_json:
                json.dump(str(self.denoised_ratio_output), doc_json)
        if not os.path.isfile(str(directory + '/denoised_ratio_d_output.json')):
            with open(str(directory + '/denoised_ratio_d_output.json'), 'w') as doc_json:
                json.dump(str(self.denoised_ratio_d_output), doc_json)

        with open(str(directory + '/good_seq_daughters.json'), 'w') as doc_json:
            json.dump(str(self.good_seq_2), doc_json)

        with open(str(directory + '/output_info_daughters.json'), 'w') as doc_json:
            json.dump(str(self.output_info_2), doc_json)

        with open(str(directory + '/denoised_d_daughters.json'), 'w') as doc_json:
            json.dump(str(self.denoised_d_output_2), doc_json)

        with open(str(directory + '/denoised_ratio_daughters.json'), 'w') as doc_json:
            json.dump(str(self.denoised_ratio_output_2), doc_json)

        with open(str(directory + '/denoised_ratio_d_output_daughters.json'), 'w') as doc_json:
            json.dump(str(self.denoised_ratio_d_output_2), doc_json)
