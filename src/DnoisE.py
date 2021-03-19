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
        de.min_mother = min(de.data_initial.loc[:, de.count]) / (1 / 2) ** (de.alpha * 1 * min(de.Ad1, de.Ad2, de.Ad3) + 1)
        print('min_mother equals to %s' % de.min_mother)
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

        # create a runned_list for A seqs with two columns
        de.runned_list = de.data_initial.loc[(de.data_initial.loc[:, de.count] >= de.min_mother), ['id', de.count]]

        # runned ---> FALSE
        de.runned_list['runned'] = np.repeat(False, de.runned_list.shape[0])

        # daughter --> FALSE
        de.runned_list['daughter'] = np.repeat(False, de.runned_list.shape[0])

        de.possibleMothers = de.runned_list.shape[0]
        de.not_possibleMothers = de.data_initial.shape[0]

        ###
        ###
        # here starts the main function
        ###
        ###

        de.denoised_d_output = []
        de.denoised_ratio_output = []
        de.denoised_ratio_d_output = []
        de.output_info = []
        de.good_seq = []
        de.abund_col_names.insert(0, de.count)

        ##############################################################################################
        # for each seq testing (possible Daughter, pD)
        # here is presented 2 different situations of similar algorithm but changing the final testing
        # there are two main ways to do it. test first whether a seq is A. or B. and run all the following
        # programm or test whether a seq is A. or B. to run the last part of the code

    if not de.parted:
        if de.entropy:
            print("running possible mothers")
            for pos in tqdm(range(de.possibleMothers)):
                [de.good_seq[len(de.good_seq):],
                 de.output_info[len(de.output_info):],
                 de.denoised_d_output[len(de.denoised_d_output):],
                 de.denoised_ratio_output[len(de.denoised_ratio_output):],
                 de.denoised_ratio_d_output[len(de.denoised_ratio_d_output):]] = de.denoising_Adcorrected(pos)

            print("running all the other seqs")

            if de.cores > 1:
                de.quartiles_runned()
                pool = mp.Pool(de.cores)

                if de.possibleMothers != de.not_possibleMothers:
                    [de.good_seq_2,
                     de.output_info_2,
                     de.denoised_d_output_2,
                     de.denoised_ratio_output_2,
                     de.denoised_ratio_d_output_2] = zip(
                        *pool.map(de.denoising_Adcorrected_parallel, [pos for pos in range(
                            de.possibleMothers, de.not_possibleMothers)]))

                pool.close()
                del (pool, de.cores, de.alpha, de.min_mother, de.possibleMothers, de.not_possibleMothers)
            else:
                de.denoised_d_output_2 = []
                de.denoised_ratio_output_2 = []
                de.denoised_ratio_d_output_2 = []
                de.output_info_2 = []
                de.good_seq_2 = []
                for pos in tqdm(range(de.possibleMothers, de.not_possibleMothers)):
                    [de.good_seq_2[len(de.good_seq_2):],
                     de.output_info_2[len(de.output_info_2):],
                     de.denoised_d_output_2[len(de.denoised_d_output_2):],
                     de.denoised_ratio_output_2[len(de.denoised_ratio_output_2):],
                     de.denoised_ratio_d_output_2[len(de.denoised_ratio_d_output_2):]] = de.denoising_Adcorrected(pos)

                del (de.cores, de.alpha, de.min_mother, de.possibleMothers, de.not_possibleMothers)
        else:
            print("running possible mothers")
            for pos in tqdm(range(de.possibleMothers)):
                [de.good_seq[len(de.good_seq):],
                 de.output_info[len(de.output_info):],
                 de.denoised_d_output[len(de.denoised_d_output):],
                 de.denoised_ratio_output[len(de.denoised_ratio_output):],
                 de.denoised_ratio_d_output[len(de.denoised_ratio_d_output):]] = de.denoising(pos)

            print("running all the other seqs")

            if de.cores > 1:
                de.quartiles_runned()
                pool = mp.Pool(de.cores)

                if de.possibleMothers != de.not_possibleMothers:
                    [de.good_seq_2,
                     de.output_info_2,
                     de.denoised_d_output_2,
                     de.denoised_ratio_output_2,
                     de.denoised_ratio_d_output_2] = zip(*pool.map(de.denoising_parallel, [pos for pos in
                                                                                           range(de.possibleMothers,
                                                                                                 de.not_possibleMothers)]))

                pool.close()
                del (pool, de.cores, de.alpha, de.min_mother, de.possibleMothers, de.not_possibleMothers)
            else:
                de.denoised_d_output_2 = []
                de.denoised_ratio_output_2 = []
                de.denoised_ratio_d_output_2 = []
                de.output_info_2 = []
                de.good_seq_2 = []
                for pos in tqdm(range(de.possibleMothers, de.not_possibleMothers)):
                    [de.good_seq_2[len(de.good_seq_2):],
                     de.output_info_2[len(de.output_info_2):],
                     de.denoised_d_output_2[len(de.denoised_d_output_2):],
                     de.denoised_ratio_output_2[len(de.denoised_ratio_output_2):],
                     de.denoised_ratio_d_output_2[len(de.denoised_ratio_d_output_2):]] = de.denoising(pos)

                del (de.cores, de.alpha, de.min_mother, de.possibleMothers, de.not_possibleMothers)
    if de.parted:
        if de.part == 1:
            print("running possible mothers")
            if de.entropy:
                for pos in tqdm(range(de.possibleMothers)):
                    [de.good_seq[len(de.good_seq):],
                     de.output_info[len(de.output_info):],
                     de.denoised_d_output[len(de.denoised_d_output):],
                     de.denoised_ratio_output[len(de.denoised_ratio_output):],
                     de.denoised_ratio_d_output[len(de.denoised_ratio_d_output):]] = de.denoising_Adcorrected(pos)
            else:
                for pos in tqdm(range(de.possibleMothers)):
                    [de.good_seq[len(de.good_seq):],
                     de.output_info[len(de.output_info):],
                     de.denoised_d_output[len(de.denoised_d_output):],
                     de.denoised_ratio_output[len(de.denoised_ratio_output):],
                     de.denoised_ratio_d_output[len(de.denoised_ratio_d_output):]] = de.denoising(pos)

            de.write_variables()
            exit()
        else:
            print("running all the other seqs")
            print('running part 2')
            de.read_variables()
            if de.entropy:
                if de.cores > 1:
                    de.quartiles_runned()
                    pool = mp.Pool(de.cores)

                    if de.possibleMothers != de.not_possibleMothers:
                        [de.good_seq_2,
                         de.output_info_2,
                         de.denoised_d_output_2,
                         de.denoised_ratio_output_2,
                         de.denoised_ratio_d_output_2] = zip(
                            *pool.map(de.denoising_Adcorrected_parallel, [pos for pos in
                                                                          range(
                                                                              de.possibleMothers,
                                                                              de.not_possibleMothers)]))

                    pool.close()
                    del (pool, de.cores, de.alpha, de.min_mother, de.possibleMothers, de.not_possibleMothers)
                else:
                    de.denoised_d_output_2 = []
                    de.denoised_ratio_output_2 = []
                    de.denoised_ratio_d_output_2 = []
                    de.output_info_2 = []
                    de.good_seq_2 = []
                    for pos in tqdm(range(de.possibleMothers, de.not_possibleMothers)):
                        [de.good_seq_2[len(de.good_seq_2):],
                         de.output_info_2[len(de.output_info_2):],
                         de.denoised_d_output_2[len(de.denoised_d_output_2):],
                         de.denoised_ratio_output_2[len(de.denoised_ratio_output_2):],
                         de.denoised_ratio_d_output_2[len(de.denoised_ratio_d_output_2):]] = de.denoising_Adcorrected(
                            pos)

                    del (de.cores, de.alpha, de.min_mother, de.possibleMothers, de.not_possibleMothers)
            else:
                if de.cores > 1:
                    de.quartiles_runned()
                    pool = mp.Pool(de.cores)

                    # pool.map(denoising, [pos for pos in range(data_initial.shape[0])])
                    if de.possibleMothers != de.not_possibleMothers:
                        [de.good_seq_2,
                         de.output_info_2,
                         de.denoised_d_output_2,
                         de.denoised_ratio_output_2,
                         de.denoised_ratio_d_output_2] = zip(*pool.map(de.denoising_parallel, [pos for pos in
                                                                                               range(de.possibleMothers,
                                                                                                     de.not_possibleMothers)]))

                    pool.close()
                    del (pool, de.cores, de.alpha, de.min_mother, de.possibleMothers, de.not_possibleMothers)
                else:
                    de.denoised_d_output_2 = []
                    de.denoised_ratio_output_2 = []
                    de.denoised_ratio_d_output_2 = []
                    de.output_info_2 = []
                    de.good_seq_2 = []
                    for pos in tqdm(range(de.possibleMothers, de.not_possibleMothers)):
                        [de.good_seq_2[len(de.good_seq_2):],
                         de.output_info_2[len(de.output_info_2):],
                         de.denoised_d_output_2[len(de.denoised_d_output_2):],
                         de.denoised_ratio_output_2[len(de.denoised_ratio_output_2):],
                         de.denoised_ratio_d_output_2[len(de.denoised_ratio_d_output_2):]] = de.denoising(pos)

                    del (de.cores, de.alpha, de.min_mother, de.possibleMothers, de.not_possibleMothers)

    print("writing files")

    de.write_outputs()

if de.part == 3:
    de.read_variables2()

if de.entropy:
    de.MOTUoutfile = str(de.MOTUoutfile + '_Adcorr')

print('writing output_info')
if de.possibleMothers > 0:
    fieldnames = de.output_info[0].keys()
else:
    fieldnames = de.output_info_2[0].keys()

with open(str(de.MOTUoutfile + '_denoising_info.csv'), 'w') as output_file:
    dict_writer = csv.DictWriter(output_file, fieldnames=fieldnames)
    dict_writer.writeheader()
    dict_writer.writerows(list(de.output_info + list(de.output_info_2)))
if de.output_info_2 == []:
    del(dict_writer, de.output_info, fieldnames)
else:
    del (dict_writer, de.output_info, de.output_info_2, fieldnames)


de.data_initial = de.data_initial.set_index(de.data_initial.loc[:, 'id'])

good_mothers = de.data_initial.loc[de.good_seq + [False] * len(de.good_seq_2)][de.first_col_names + de.abund_col_names + [de.seq]]
good_daughters = de.data_initial.loc[[False] * len(de.good_seq) + list(de.good_seq_2)][de.first_col_names + de.abund_col_names + [de.seq]]


print('writing output_d')
# writing d
denoised_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])

for mother in good_mothers.loc[:, 'id']:
    row = [good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
           list(de.data_initial.loc[[de.mother_id(x, mother) for x in list(
               de.denoised_d_output + list(de.denoised_d_output_2))], de.abund_col_names].sum(0)) +
           good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
    row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
    denoised_d = denoised_d.append(row, ignore_index=True)

denoised_d = pd.concat([denoised_d, good_daughters], ignore_index=True)
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

if de.denoised_d_output_2 == []:
    del(denoised_d, de.denoised_d_output)
else:
    del(denoised_d, de.denoised_d_output, de.denoised_d_output_2)


print('writing output_ratio')
# writing ratio
denoised_ratio = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])

for mother in good_mothers.loc[:, 'id']:
    row = [good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
           list(de.data_initial.loc[[de.mother_id(x, mother) for x in list(
               de.denoised_ratio_output + list(de.denoised_ratio_output_2))], de.abund_col_names].sum(0)) +
           good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
    row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
    denoised_ratio = denoised_ratio.append(row, ignore_index=True)

denoised_ratio = pd.concat([denoised_ratio, good_daughters], ignore_index=True)

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

if de.denoised_ratio_output_2 == []:
    del(denoised_ratio, de.denoised_ratio_output)
else:
    del(denoised_ratio, de.denoised_ratio_output, de.denoised_ratio_output_2)


print('writing output_ratio_d')
# writing ratio_d
denoised_ratio_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])

for mother in good_mothers.loc[:, 'id']:
    row = [good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
           list(de.data_initial.loc[[de.mother_id(x, mother) for x in list(
               de.denoised_ratio_d_output + list(de.denoised_ratio_d_output_2))], de.abund_col_names].sum(0)) +
           good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
    row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
    denoised_ratio_d = denoised_ratio_d.append(row, ignore_index=True)

denoised_ratio_d = pd.concat([denoised_ratio_d, good_daughters], ignore_index=True)
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

if de.denoised_ratio_d_output_2 == []:
    del(denoised_ratio_d, de.denoised_ratio_d_output)
else:
    del(denoised_ratio_d, de.denoised_ratio_d_output, de.denoised_ratio_d_output_2)

print('done')