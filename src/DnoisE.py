# import all the modules that will be used
import pandas as pd
import numpy as np
import multiprocessing as mp
import sys
from tqdm import tqdm
import itertools
from denoise_functions import denoise_functions
from transform_data import *

de = denoise_functions()

# Get full command-line arguments
argument_list = sys.argv

# Keep all but the first
argument_list = argument_list[1:]

print(argument_list)
de.read_parameters(argument_list)
import_data(de)
transform_data(de)

if de.part != 2:
    # preparing data

    if de.part == 1:

        ###
        ###
        # here starts the main function
        ###
        ###
        #first position is always the same
        if de.output_type == 'ratio':
            de.denoised_ratio_output = [de.data_initial.loc[0, 'id']]
        else:
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
            if de.output_type == 'ratio':
                while len(de.runned_list) < de.data_initial.shape[0]:
                    de.min_mother = de.runned_list[-1].get(de.count) * de.max_ratio
                    run_to = sum(de.data_initial.loc[:, de.count] > de.min_mother)
                    if run_to == len(de.runned_list):
                        run_to += 1
                    de.quartiles_runned()
                    print('running until %s reads' % de.min_mother)
                    print(len(de.runned_list) / de.data_initial.shape[0] * 100, '%')
                    pool = mp.Pool(de.cores)
                    [de.good_seq_2,
                     de.output_info_2,
                     de.denoised_ratio_output_2,
                     de.runned_list_2] = zip(*pool.map(de.denoising_Adcorrected_parallel_ratio, [pos for pos in
                                                                                                 range(len(de.runned_list),
                                                                                                       run_to)]))
                    pool.close()
                    del pool
                    de.good_seq.extend(list(de.good_seq_2))
                    de.output_info.extend(list(de.output_info_2))
                    de.denoised_ratio_output.extend(list(de.denoised_ratio_output_2))
                    de.runned_list.extend(list(de.runned_list_2))
                    print('runned until %s reads' % de.min_mother)
                    print(len(de.runned_list) / de.data_initial.shape[0] * 100, '%')
                    del (de.good_seq_2, de.output_info_2, de.denoised_ratio_output_2, de.runned_list_2)

            else:
                while len(de.runned_list) < de.data_initial.shape[0]:
                    de.min_mother = de.runned_list[-1].get(de.count) * de.max_ratio
                    run_to = sum(de.data_initial.loc[:, de.count] > de.min_mother)
                    if run_to == len(de.runned_list):
                        run_to += 1
                    de.quartiles_runned()
                    print('running until %s reads' % de.min_mother)
                    print(len(de.runned_list) / de.data_initial.shape[0] * 100, '%')
                    pool = mp.Pool(de.cores)
                    [de.good_seq_2,
                     de.output_info_2,
                     de.denoised_d_output_2,
                     de.denoised_ratio_output_2,
                     de.denoised_ratio_d_output_2,
                     de.runned_list_2] = zip(*pool.map(de.denoising_Adcorrected_parallel, [pos for pos in
                                                                                           range(len(de.runned_list),
                                                                                                 run_to)]))
                    pool.close()
                    del pool
                    de.good_seq.extend(list(de.good_seq_2))
                    de.output_info.extend(list(de.output_info_2))
                    de.denoised_d_output.extend(list(de.denoised_d_output_2))
                    de.denoised_ratio_output.extend(list(de.denoised_ratio_output_2))
                    de.denoised_ratio_d_output.extend(list(de.denoised_ratio_d_output_2))
                    de.runned_list.extend(list(de.runned_list_2))
                    print('runned until %s reads' % de.min_mother)
                    print(len(de.runned_list) / de.data_initial.shape[0] * 100, '%')
                    del (de.good_seq_2, de.output_info_2, de.denoised_d_output_2, de.denoised_ratio_output_2,
                         de.denoised_ratio_d_output_2, de.runned_list_2)
            if 'de.min_mother' in locals():
                del de.min_mother
            del (de.cores, de.alpha)
        else:
            if de.output_type == 'ratio':
                for pos in tqdm(range(1, de.data_initial.shape[0])):
                    [de.good_seq[len(de.good_seq):],
                     de.output_info[len(de.output_info):],
                     de.denoised_ratio_output[len(de.denoised_ratio_output):],
                     de.runned_list[len(de.runned_list):]] = de.denoising_Adcorrected_ratio(pos)

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
            if de.output_type == 'ratio':
                while len(de.runned_list) < de.data_initial.shape[0]:
                    de.min_mother = de.runned_list[-1].get(de.count) * de.max_ratio
                    run_to = sum(de.data_initial.loc[:, de.count] > de.min_mother)
                    if run_to == len(de.runned_list):
                        run_to += 1
                    de.quartiles_runned()
                    print('running until %s reads' % de.min_mother)
                    print(len(de.runned_list) / de.data_initial.shape[0] * 100, '%')
                    pool = mp.Pool(de.cores)
                    [de.good_seq_2,
                     de.output_info_2,
                     de.denoised_ratio_output_2,
                     de.runned_list_2] = zip(*pool.map(de.denoising_parallel_ratio, [pos for pos in
                                                                                     range(len(de.runned_list),
                                                                                           run_to)]))
                    pool.close()
                    del pool
                    de.good_seq.extend(list(de.good_seq_2))
                    de.output_info.extend(list(de.output_info_2))
                    de.denoised_d_output.extend(list(de.denoised_d_output_2))
                    de.denoised_ratio_output.extend(list(de.denoised_ratio_output_2))
                    de.denoised_ratio_d_output.extend(list(de.denoised_ratio_d_output_2))
                    de.runned_list.extend(list(de.runned_list_2))
                    print('runned until %s reads' % de.min_mother)
                    print(len(de.runned_list) / de.data_initial.shape[0] * 100, '%')
                    del (de.good_seq_2, de.output_info_2, de.denoised_ratio_output_2, de.runned_list_2)
            else:
                while len(de.runned_list) < de.data_initial.shape[0]:
                    de.min_mother = de.runned_list[-1].get(de.count) * de.max_ratio
                    run_to = sum(de.data_initial.loc[:, de.count] > de.min_mother)
                    if run_to == len(de.runned_list):
                        run_to += 1
                    de.quartiles_runned()
                    print('running until %s reads' % de.min_mother)
                    print(len(de.runned_list) / de.data_initial.shape[0] * 100, '%')
                    pool = mp.Pool(de.cores)
                    [de.good_seq_2,
                     de.output_info_2,
                     de.denoised_d_output_2,
                     de.denoised_ratio_output_2,
                     de.denoised_ratio_d_output_2,
                     de.runned_list_2] = zip(*pool.map(de.denoising_parallel, [pos for pos in
                                                                               range(len(de.runned_list),
                                                                                     run_to)]))
                    pool.close()
                    del pool
                    de.good_seq.extend(list(de.good_seq_2))
                    de.output_info.extend(list(de.output_info_2))
                    de.denoised_d_output.extend(list(de.denoised_d_output_2))
                    de.denoised_ratio_output.extend(list(de.denoised_ratio_output_2))
                    de.denoised_ratio_d_output.extend(list(de.denoised_ratio_d_output_2))
                    de.runned_list.extend(list(de.runned_list_2))
                    print('runned until %s reads' % de.min_mother)
                    print(len(de.runned_list) / de.data_initial.shape[0] * 100, '%')
                    del (de.good_seq_2, de.output_info_2, de.denoised_d_output_2, de.denoised_ratio_output_2,
                         de.denoised_ratio_d_output_2, de.runned_list_2)
            if 'de.min_mother' in locals():
                del de.min_mother
            del (de.cores, de.alpha)
        else:
            if de.output_type == 'ratio':
                for pos in tqdm(range(1, de.data_initial.shape[0])):
                    [de.good_seq[len(de.good_seq):],
                     de.output_info[len(de.output_info):],
                     de.denoised_ratio_output[len(de.denoised_ratio_output):],
                     de.runned_list[len(de.runned_list):]] = de.denoising_ratio(pos)
            else:
                for pos in tqdm(range(1, de.data_initial.shape[0])):
                    [de.good_seq[len(de.good_seq):],
                     de.output_info[len(de.output_info):],
                     de.denoised_d_output[len(de.denoised_d_output):],
                     de.denoised_ratio_output[len(de.denoised_ratio_output):],
                     de.denoised_ratio_d_output[len(de.denoised_ratio_d_output):],
                     de.runned_list[len(de.runned_list):]] = de.denoising(pos)

    de.write_variables()

if de.part == 2:
    de.read_variables2()

if de.entropy:
    de.MOTUoutfile = str(de.MOTUoutfile + '_Adcorr')

print('writing output_info')
# fieldnames = de.output_info[0].keys()
#
#
# with open(str(de.MOTUoutfile + '_denoising_info.csv'), 'w') as output_file:
#     dict_writer = csv.DictWriter(output_file, fieldnames=fieldnames)
#     dict_writer.writeheader()
#     dict_writer.writerows(list(de.output_info))
# del(dict_writer, de.output_info, fieldnames)

de.output_info = pd.DataFrame.from_dict(de.output_info)
de.output_info.to_csv(str(de.MOTUoutfile + '_denoising_info.csv'), index=False)

if (de.output_type == 'ratio') or (de.output_type == 'all'):
    mothers_ratio = de.output_info.mother_ratio.unique()[1:]
if (de.output_type == 'd') or (de.output_type == 'all'):
    mothers_d = de.output_info.mother_d.unique()[1:]
if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
    mothers_ratio_d = de.output_info.mother_xavier_criteria.unique()[1:]

del de.output_info


de.data_initial = de.data_initial.set_index(de.data_initial.loc[:, 'id'])


if (de.output_type == 'ratio') or (de.output_type == 'all'):
    de.good_mothers = de.data_initial.loc[de.good_seq][de.first_col_names + de.abund_col_names + [de.seq]]
    print('writing output_ratio')
    # writing ratio
    if de.cores > 1:
        pool = mp.Pool(de.cores)
        [row] = zip(*pool.map(de.write_output_ratio, [mother for mother in mothers_ratio]))
        pool.close()
        del pool
        denoised_ratio = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
        de.good_mothers = de.good_mothers.drop(index=mothers_ratio)
        denoised_ratio = denoised_ratio.append(de.good_mothers, ignore_index=True)
        denoised_ratio = denoised_ratio.sort_values([de.count], axis=0, ascending=False)
    else:
        # denoised_ratio = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
        # for mother in tqdm(de.good_mothers.loc[:, 'id']):
        #     row = [de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
        #            list(de.data_initial.loc[[de.mother_id(x, mother) for x in list(
        #                de.denoised_ratio_output)], de.abund_col_names].sum(0)) +
        #            de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
        #     row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
        #     denoised_ratio = denoised_ratio.append(row, ignore_index=True)

        denoised_ratio = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
        for mother in tqdm(mothers_ratio):
            row = [
                de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
                list(de.data_initial.loc[list(pd.Series(de.denoised_ratio_output) == mother), de.abund_col_names].sum(0)) +
                de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
            row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            denoised_ratio = denoised_ratio.append(row, ignore_index=True)
            de.good_mothers = de.good_mothers.drop(index=mother)
        denoised_ratio = denoised_ratio.append(de.good_mothers, ignore_index=True)
        denoised_ratio = denoised_ratio.sort_values([de.count], axis=0, ascending=False)
    if 'row' in locals():
        del row
    del de.good_mothers, mothers_ratio, de.denoised_ratio_output

    if de.fasta_output:
        denoised_ratio = denoised_ratio.to_dict(orient='index')
        ofile = open(str(de.MOTUoutfile + '_denoised_ratio.fasta'), "w")
        for i in tqdm(range(len(denoised_ratio))):
            ofile.write(">" + denoised_ratio[i]['id'] + ';size=' + str(denoised_ratio[i][de.count]) +
                        ";\n" + denoised_ratio[i][de.seq].upper() + "\n")
        # do not forget to close it
        ofile.close()
    else:
        denoised_ratio.to_csv(str(de.MOTUoutfile + '_denoised_ratio.csv'), index=False)

    del denoised_ratio

if (de.output_type == 'd') or (de.output_type == 'all'):
    de.good_mothers = de.data_initial.loc[de.good_seq][de.first_col_names + de.abund_col_names + [de.seq]]
    print('writing output_d')
    # writing d
    if de.cores > 1:
        pool = mp.Pool(de.cores)
        [row] = zip(*pool.map(de.write_output_d, [mother for mother in mothers_d]))
        pool.close()
        del pool
        denoised_d = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
        de.good_mothers = de.good_mothers.drop(index=mothers_d)
        denoised_d = denoised_d.append(de.good_mothers, ignore_index=True)
        denoised_d = denoised_d.sort_values([de.count], axis=0, ascending=False)
    else:
        denoised_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
        for mother in tqdm(mothers_d):
            row = [
                de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
                list(de.data_initial.loc[list(pd.Series(de.denoised_d_output) == mother), de.abund_col_names].sum(
                    0)) +
                de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
            row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            denoised_d = denoised_d.append(row, ignore_index=True)
            de.good_mothers = de.good_mothers.drop(index=mother)
        denoised_d = denoised_d.append(de.good_mothers, ignore_index=True)
        denoised_d = denoised_d.sort_values([de.count], axis=0, ascending=False)
    if 'row' in locals():
        del row
    del de.good_mothers, mothers_d, de.denoised_d_output

    if de.fasta_output:
        denoised_d = denoised_d.to_dict(orient='index')
        ofile = open(str(de.MOTUoutfile + '_denoised_d.fasta'), "w")
        for i in tqdm(range(len(denoised_d))):
            ofile.write(">" + denoised_d[i]['id'] + ';size=' + str(denoised_d[i][de.count]) + ";\n" +
                        denoised_d[i][de.seq].upper() + "\n")
        # do not forget to close it
        ofile.close()
    else:
        denoised_d.to_csv(str(de.MOTUoutfile + '_denoised_d.csv'), index=False)

    del denoised_d




if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
    de.good_mothers = de.data_initial.loc[de.good_seq][de.first_col_names + de.abund_col_names + [de.seq]]
    print('writing output_ratio_d')
    # writing ratio_d
    if de.cores > 1:
        pool = mp.Pool(de.cores)
        [row] = zip(*pool.map(de.write_output_ratio_d, [mother for mother in mothers_ratio_d]))
        pool.close()
        del pool
        denoised_ratio_d = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
        de.good_mothers = de.good_mothers.drop(index=mothers_ratio_d)
        denoised_ratio_d = denoised_ratio_d.append(de.good_mothers, ignore_index=True)
        denoised_ratio_d = denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
    else:
        denoised_ratio_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
        for mother in tqdm(mothers_ratio_d):
            row = [
                de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
                list(de.data_initial.loc[list(pd.Series(de.denoised_ratio_d_output) == mother), de.abund_col_names].sum(
                    0)) +
                de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
            row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            denoised_ratio_d = denoised_ratio_d.append(row, ignore_index=True)
            de.good_mothers = de.good_mothers.drop(index=mother)
        denoised_ratio_d = denoised_ratio_d.append(de.good_mothers, ignore_index=True)
        denoised_ratio_d = denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
    if 'row' in locals():
        del row
    del de.good_mothers, mothers_ratio_d, de.denoised_ratio_d_output

    if de.fasta_output:
        denoised_ratio_d = denoised_ratio_d.to_dict(orient='index')
        ofile = open(str(de.MOTUoutfile + '_denoised_ratio_d.fasta'), "w")
        for i in tqdm(range(len(denoised_ratio_d))):
            ofile.write(">" + denoised_ratio_d[i]['id'] + ';size=' + str(denoised_ratio_d[i][de.count]) + ";\n" +
                        denoised_ratio_d[i][de.seq].upper() + "\n")
        # do not forget to close it
        ofile.close()
    else:
        denoised_ratio_d.to_csv(str(de.MOTUoutfile + '_denoised_ratio_d.csv'), index=False)

    del denoised_ratio_d

print('done')
