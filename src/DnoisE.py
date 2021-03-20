#import all the modules that will be used
import pandas as pd
import numpy as np
import multiprocessing as mp
import csv
import sys
from tqdm import tqdm
from needed_modules.Bio import SeqIO
# from Bio import SeqIO
import re
import stats
import itertools
from denoise_functions import denoise_functions

de = denoise_functions()

# Get full command-line arguments
full_cmd_arguments = sys.argv

# Keep all but the first
# argument_list = ['-i', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/PHY1bis_final_subset.csv', '-o', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/PHY1bis_final_subset.csv_Adcorr_nou',
#                  '-P', '3', '-f', 'F', '-F', 'F', '-c', '2', '-n', 'reads', '-a', '5', '-q', 'seq', '-p', '2', '-e', '0.4727,0.2266,1.0212', '-y', 'T']
argument_list = full_cmd_arguments[1:]

print(argument_list)
de.read_parameters(argument_list)

if de.part != 3:
    # preparing data

    if de.part == 1:
        print('reading input file')
        if de.fasta:
            # data_initial = pd.DataFrame()
            for fastaseq in SeqIO.parse(de.MOTUfile, "fasta"):
                seq_seq = fastaseq.seq._data
                print(fastaseq.description)
                seq_id = re.findall('(.+);', fastaseq.description)[0]
                seq_count = re.findall('size=(.+);', fastaseq.description)[0]
                de.data_initial = pd.concat(
                    [de.data_initial, pd.DataFrame({"id": seq_id, "count": int(seq_count), "sequence": [seq_seq]})])
        else:
            de.data_initial = pd.read_csv(de.MOTUfile, sep=de.sep)

        print('input file read')

        # define alpha and minimum min_mother reads

        # de.min_mother = min(de.data_initial.loc[:, de.count]) / (1 / 2) ** (de.alpha * 1 + 1 + (min(de.Ad1, de.Ad2, de.Ad3)-1) * 1) # for sum version
        # de.min_mother = min(de.data_initial.loc[:, de.count]) / (1 / 2) ** (de.alpha * 1 * min(de.Ad1, de.Ad2, de.Ad3) + 1)
        # print('min_mother equals to %s' % de.min_mother)
        print('and Ad corr:')
        print(de.Ad1, de.Ad2, de.Ad3)
        # the program should work with different seq_length, if not, filter and get de mode

        # obtain a column with total reads per seq.
        if not de.justcount:
            de.abund_col_names = list(de.data_initial.columns)[(de.start - 1):de.end]
            de.first_col_names = list(de.data_initial.columns)[0:(de.start - 2)]
            de.data_initial.loc[:, de.count] = de.data_initial.loc[:, de.abund_col_names].sum(axis=1)
        else:
            de.first_col_names = ['id']

        # sort by total reads
        de.data_initial = de.data_initial.sort_values([de.count], axis=0, ascending=False)

        # delete seqs with 0 reads
        de.data_initial = de.data_initial.loc[(de.data_initial.loc[:, de.count] != 0)]
        if de.entropy:
            # remove seq out of mode length
            de.data_initial.index = list(range(de.data_initial.shape[0]))
            seq_length = []
            seq_length_per_read = []
            for i in list(range(de.data_initial.shape[0])):
                i_seq = de.data_initial.loc[i, de.seq]
                i_count = de.data_initial.loc[i, de.count]
                seq_length.append(len(i_seq))
                seq_length_per_read.append([len(i_seq)]*i_count)
            seq_length_per_read = list(itertools.chain.from_iterable(seq_length_per_read))
            try:
                de.data_initial = de.data_initial.loc[(np.asarray(seq_length) == stats.mode(seq_length_per_read))]
            except:
                print('MOTU no available to run with Entropy. Equal number of seqs with different seq length')
                exit()
            del seq_length, seq_length_per_read

        # reorder index
        de.data_initial.index = list(range(de.data_initial.shape[0]))

        # analising each seq. two main groups.
        # A. possible mothers
        # B. less than min_mother reads. cannot be mothers



        ###
        ###
        # here starts the main function
        ###
        ###
        #first position is always the same
        de.denoised_d_output = [de.data_initial.loc[0, 'id']]
        de.denoised_ratio_output = [de.data_initial.loc[0, 'id']]
        de.denoised_ratio_d_output = [de.data_initial.loc[0, 'id']]
        if de.entropy:
            de.output_info = [{'daughter': de.data_initial.loc[0, 'id'], 'mother_d': None, 'd': None,
                               'mother_ratio': None, 'ratio': None,
                               'mother_xavier_criteria': None, 'xavier_criteria': None,
                               'difpos1': None, 'difpos2': None, 'difpos3': None}]
        else:
            de.output_info = [{'daughter': de.data_initial.loc[0, 'id'], 'mother_d': None, 'd': None,
                               'mother_ratio': None, 'ratio': None,
                               'mother_xavier_criteria': None, 'xavier_criteria': None}]
        de.good_seq = [True]
        de.abund_col_names.insert(0, de.count)
        de.runned_list = [{'id': de.data_initial.loc[0, 'id'], de.count: de.data_initial.loc[0, de.count],
                           'runned': True, 'daughter': False}]

        ##############################################################################################
        # for each seq testing (possible Daughter, pD)
        # here is presented 2 different situations of similar algorithm but changing the final testing
        # there are two main ways to do it. test first whether a seq is A. or B. and run all the following
        # programm or test whether a seq is A. or B. to run the last part of the code
    if de.entropy:
        if de.cores > 1:
            while len(de.runned_list) < de.data_initial.shape[0]:
                de.min_mother = de.runned_list[-1].get(de.count) * (1 / 2) ** (
                            de.alpha * 1 * min(de.Ad1, de.Ad2, de.Ad3) + 1)
                run_to = sum(de.data_initial.loc[:, de.count] > de.min_mother)
                if run_to == len(de.runned_list):
                    run_to += 1
                de.quartiles_runned()
                pool = mp.Pool(de.cores)
                print('running until %s reads' % de.min_mother)
                print(len(de.runned_list) / de.data_initial.shape[0] * 100, '%')
                [de.good_seq_2,
                 de.output_info_2,
                 de.denoised_d_output_2,
                 de.denoised_ratio_output_2,
                 de.denoised_ratio_d_output_2,
                 de.runned_list_2] = zip(*pool.map(de.denoising_Adcorrected_parallel, [pos for pos in
                                                                                       range(len(de.runned_list),
                                                                                             run_to)]))
                de.good_seq.extend(list(de.good_seq_2))
                de.output_info.extend(list(de.output_info_2))
                de.denoised_d_output.extend(list(de.denoised_d_output_2))
                de.denoised_ratio_output.extend(list(de.denoised_ratio_output_2))
                de.denoised_ratio_d_output.extend(list(de.denoised_ratio_d_output_2))
                de.runned_list.extend(list(de.runned_list_2))

                pool.close()
            del (pool, de.cores, de.alpha, de.min_mother)
        else:
            for pos in tqdm(range(1, de.data_initial.shape[0])):
                [de.good_seq[len(de.good_seq):],
                 de.output_info[len(de.output_info):],
                 de.denoised_d_output[len(de.denoised_d_output):],
                 de.denoised_ratio_output[len(de.denoised_ratio_output):],
                 de.denoised_ratio_d_output[len(de.denoised_ratio_d_output):],
                 de.runned_list[len(de.runned_list):]] = de.denoising_Adcorrected(pos)

    else:
        if de.cores > 1:
            while len(de.runned_list) < de.data_initial.shape[0]:
                de.min_mother = de.runned_list[-1].get(de.count) * (1 / 2) ** (de.alpha * 1 * min(de.Ad1, de.Ad2, de.Ad3) + 1)
                run_to = sum(de.data_initial.loc[:, de.count] > de.min_mother)
                if run_to == len(de.runned_list):
                    run_to += 1
                de.quartiles_runned()
                pool = mp.Pool(de.cores)
                print('running until %s reads' % de.min_mother)
                print(len(de.runned_list)/de.data_initial.shape[0] *100, '%')
                [de.good_seq_2,
                 de.output_info_2,
                 de.denoised_d_output_2,
                 de.denoised_ratio_output_2,
                 de.denoised_ratio_d_output_2,
                 de.runned_list_2] = zip(*pool.map(de.denoising_parallel, [pos for pos in
                                                                           range(len(de.runned_list),
                                                                                 run_to)]))
                de.good_seq.extend(list(de.good_seq_2))
                de.output_info.extend(list(de.output_info_2))
                de.denoised_d_output.extend(list(de.denoised_d_output_2))
                de.denoised_ratio_output.extend(list(de.denoised_ratio_output_2))
                de.denoised_ratio_d_output.extend(list(de.denoised_ratio_d_output_2))
                de.runned_list.extend(list(de.runned_list_2))

                pool.close()
            del (pool, de.cores, de.alpha, de.min_mother)
        else:
            for pos in tqdm(range(1,de.data_initial.shape[0])):
                [de.good_seq[len(de.good_seq):],
                 de.output_info[len(de.output_info):],
                 de.denoised_d_output[len(de.denoised_d_output):],
                 de.denoised_ratio_output[len(de.denoised_ratio_output):],
                 de.denoised_ratio_d_output[len(de.denoised_ratio_d_output):],
                 de.runned_list[len(de.runned_list):]] = de.denoising(pos)

    de.write_variables()

if de.part == 3:
    de.read_variables2()

if de.entropy:
    de.MOTUoutfile = str(de.MOTUoutfile + '_Adcorr')

print('writing output_info')
fieldnames = de.output_info[0].keys()


with open(str(de.MOTUoutfile + '_denoising_info.csv'), 'w') as output_file:
    dict_writer = csv.DictWriter(output_file, fieldnames=fieldnames)
    dict_writer.writeheader()
    dict_writer.writerows(list(de.output_info))
del(dict_writer, de.output_info, fieldnames)


de.data_initial = de.data_initial.set_index(de.data_initial.loc[:, 'id'])

good_mothers = de.data_initial.loc[de.good_seq][de.first_col_names + de.abund_col_names + [de.seq]]

print('writing output_d')
# writing d
denoised_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])

for mother in good_mothers.loc[:, 'id']:
    row = [good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
           list(de.data_initial.loc[[de.mother_id(x, mother) for x in list(
               de.denoised_d_output)], de.abund_col_names].sum(0)) +
           good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
    row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
    denoised_d = denoised_d.append(row, ignore_index=True)


if de.fasta_output:
    denoised_d = denoised_d.to_dict(orient='index')
    ofile = open(str(de.MOTUoutfile + '_denoised_d.fasta'), "w")
    for i in range(len(denoised_d)):
        ofile.write(">" + denoised_d[i]['id'] + ';size=' + str(denoised_d[i]['count']) + ";\n" +
                    denoised_d[i]['sequence'].upper() + "\n")
    # do not forget to close it
    ofile.close()
else:
    denoised_d.to_csv(str(de.MOTUoutfile + '_denoised_d.csv'), index=False)

del(denoised_d, de.denoised_d_output)


print('writing output_ratio')
# writing ratio
denoised_ratio = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])

for mother in good_mothers.loc[:, 'id']:
    row = [good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
           list(de.data_initial.loc[[de.mother_id(x, mother) for x in list(
               de.denoised_ratio_output)], de.abund_col_names].sum(0)) +
           good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
    row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
    denoised_ratio = denoised_ratio.append(row, ignore_index=True)

if de.fasta_output:
    denoised_ratio = denoised_ratio.to_dict(orient='index')
    ofile = open(str(de.MOTUoutfile + '_denoised_ratio.fasta'), "w")
    for i in range(len(denoised_ratio)):
        ofile.write(">" + denoised_ratio[i]['id'] + ';size=' + str(denoised_ratio[i]['count']) +
                    ";\n" + denoised_ratio[i]['sequence'].upper() + "\n")
    # do not forget to close it
    ofile.close()
else:
    denoised_ratio.to_csv(str(de.MOTUoutfile + '_denoised_ratio.csv'), index=False)

del(denoised_ratio, de.denoised_ratio_output)


print('writing output_ratio_d')
# writing ratio_d
denoised_ratio_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])

for mother in good_mothers.loc[:, 'id']:
    row = [good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
           list(de.data_initial.loc[[de.mother_id(x, mother) for x in list(
               de.denoised_ratio_d_output)], de.abund_col_names].sum(0)) +
           good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
    row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
    denoised_ratio_d = denoised_ratio_d.append(row, ignore_index=True)

if de.fasta_output:
    denoised_ratio_d = denoised_ratio_d.to_dict(orient='index')
    ofile = open(str(de.MOTUoutfile + '_denoised_ratio_d.fasta'), "w")
    for i in range(len(denoised_ratio_d)):
        ofile.write(">" + denoised_ratio_d[i]['id'] + ';size=' + str(denoised_ratio_d[i]['count']) +
                    ";\n" + denoised_ratio_d[i]['sequence'].upper() + "\n")
    # do not forget to close it
    ofile.close()
else:
    denoised_ratio_d.to_csv(str(de.MOTUoutfile + '_denoised_ratio_d.csv'), index=False)

del(denoised_ratio_d, de.denoised_ratio_d_output)

print('done')
