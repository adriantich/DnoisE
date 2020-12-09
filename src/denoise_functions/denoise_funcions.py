
import pandas as pd
import Levenshtein as lv
import json
from ast import literal_eval
import os
# import shutil

# not_possibleMothers = []
# possibleMothers = []
# data_inicial = []
# runned_list = []
# alpha = []
# min_mother = []
# Ad1 = []
# Ad2 = []
# Ad3 = []
# MOTUfile = ""
# good_seq = []
# output_info = []
# denoised_d_output = []
# denoised_ratio_output = []
# denoised_ratio_d_output = []
# entropy = []
# cores = []
# MOTUoutfile = ""

class denoise_functions:
    not_possibleMothers = 0
    possibleMothers = 0
    data_inicial = pd.DataFrame()
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
    inicial_pos = 2
    justcount = True
    abund_col_names = []
    first_col_names = []
    fasta_output = True

    def __init__(self):
        print("starting to denoise")

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
        pD = self.data_inicial.loc[pos, 'id']
        pDseq = self.data_inicial.loc[pos, self.seq]
        pDabund = self.data_inicial.loc[pos, self.count]
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
            pMseq = self.data_inicial.loc[pMpos, self.seq]
            pMabund = self.data_inicial.loc[pMpos, self.count]
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

        pD = self.data_inicial.loc[pos, 'id']
        pDseq = self.data_inicial.loc[pos, self.seq]
        pDabund = self.data_inicial.loc[pos, self.count]
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
            pMseq = self.data_inicial.loc[pMpos, self.seq]
            pMabund = self.data_inicial.loc[pMpos, self.count]
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
            # pDinfo = np.transpose(pd.DataFrame(self.data_inicial.loc[pos, abund_col_names]))
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
        for i, o in enumerate(seq1):
            if o != seq2[i]:
                dif_position = (i + self.inicial_pos) % 3
                if dif_position == 1:
                    difcount += self.Ad1
                elif dif_position == 2:
                    difcount += self.Ad2
                elif dif_position == 0:
                    difcount += self.Ad3
        return difcount

    def denoising_Adcorrected(self, pos):
        pD = self.data_inicial.loc[pos, 'id']
        pDseq = self.data_inicial.loc[pos, self.seq]
        pDabund = self.data_inicial.loc[pos, self.count]

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
            pMseq = self.data_inicial.loc[pMpos, self.seq]
            pMabund = self.data_inicial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > (1 / self.min_mother):
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d = self.difference(seq1=pDseq, seq2=pMseq)
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

        pD = self.data_inicial.loc[pos, 'id']
        pDseq = self.data_inicial.loc[pos, self.seq]
        pDabund = self.data_inicial.loc[pos, self.count]
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
            pMseq = self.data_inicial.loc[pMpos, self.seq]
            pMabund = self.data_inicial.loc[pMpos, self.count]
            # obtain ratio ---> total_reads pD / total_reads pM
            b_ratio = pDabund / pMabund
            # if ratio is less than minimum (1/64)
            # break comparing with more pM
            if b_ratio > (1 / self.min_mother):
                break
                # obtain d ---> external function:
                # input: must be seq of both pM and pD
                # output: d
            d = self.difference(pDseq, pMseq)
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
            # pDinfo = np.transpose(pd.DataFrame(self.data_inicial.loc[pos, abund_col_names]))
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

        with open(str(directory + '/data_inicial.json'), 'w') as doc_json:
            self.data_inicial.to_json(doc_json, orient='table')

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

        with open(str(directory + '/data_inicial.json'), 'r') as doc_json:
            self.data_inicial = pd.read_json(doc_json, orient='table')

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

        with open(str(directory + '/data_inicial.json'), 'r') as doc_json:
            self.data_inicial = pd.read_json(doc_json, orient='table')

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

        if not os.path.isfile(str(directory + '/data_inicial.json')):
            with open(str(directory + '/data_inicial.json'), 'w') as doc_json:
                self.data_inicial.to_json(doc_json, orient='table')
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
