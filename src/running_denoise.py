#!/usr/bin/env python3

"""
.. codeauthor:: Adrià Antich <adriantich@gmail.com>

This programme is called by the DnoisE.

runnin_denoise.py runs the algorithm of denoising sequences using Levenshtein distance or entropy correction.

"""

import multiprocessing as mp
from tqdm import tqdm
import itertools
import numpy as np
import entropy as en
import pandas as pd
from denoise_functions import *


def run_denoise(de):
    if de.output_type == 'ratio':
        de.denoised_ratio_output = [de.data_initial.loc[0, 'id']]
    else:
        de.denoised_d_output = [de.data_initial.loc[0, 'id']]
        de.denoised_ratio_output = [de.data_initial.loc[0, 'id']]
        de.denoised_ratio_d_output = [de.data_initial.loc[0, 'id']]

    de.output_info = [{'daughter': de.data_initial.loc[0, 'id'], 'mother_d': None, 'd': None,
                       'mother_ratio': None, 'ratio': None,
                       'mother_ratio_d': None, 'xavier_criteria': None}]

    de.good_seq = [True]
    de.abund_col_names.insert(0, de.count)
    de.run_list = [{'id': de.data_initial.loc[0, 'id'], de.count: de.data_initial.loc[0, de.count],
                    'run': True, 'daughter': False}]

    run_dnoise_testing(de)  # function in denoise_funtions.py

    de.output_info = pd.DataFrame.from_dict(de.output_info)
    de.output_info.to_csv(str(de.MOTUoutfile + '_denoising_info.csv'), index=False)

    if (de.output_type == 'ratio') or (de.output_type == 'all'):
        mothers_ratio = de.output_info.mother_ratio.unique()[1:]
    if (de.output_type == 'd') or (de.output_type == 'all'):
        mothers_d = de.output_info.mother_d.unique()[1:]
    if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
        mothers_ratio_d = de.output_info.mother_ratio_d.unique()[1:]

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
            de.denoised_ratio = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            de.good_mothers = de.good_mothers.drop(index=mothers_ratio)
            de.denoised_ratio = de.denoised_ratio.append(de.good_mothers, ignore_index=True)
            de.denoised_ratio = de.denoised_ratio.sort_values([de.count], axis=0, ascending=False)
        else:
            de.denoised_ratio = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            for mother in tqdm(mothers_ratio):
                row = [
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[
                        0] +
                    list(de.data_initial.loc[
                             list(pd.Series(de.denoised_ratio_output) == mother), de.abund_col_names].sum(0)) +
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
                row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
                de.denoised_ratio = de.denoised_ratio.append(row, ignore_index=True)
                de.good_mothers = de.good_mothers.drop(index=mother)
            de.denoised_ratio = de.denoised_ratio.append(de.good_mothers, ignore_index=True)
            de.denoised_ratio = de.denoised_ratio.sort_values([de.count], axis=0, ascending=False)
        if 'row' in locals():
            del row
        del de.good_mothers, mothers_ratio, de.denoised_ratio_output

    if (de.output_type == 'd') or (de.output_type == 'all'):
        de.good_mothers = de.data_initial.loc[de.good_seq][de.first_col_names + de.abund_col_names + [de.seq]]
        print('writing output_d')
        # writing d
        if de.cores > 1:
            pool = mp.Pool(de.cores)
            [row] = zip(*pool.map(de.write_output_d, [mother for mother in mothers_d]))
            pool.close()
            del pool
            de.denoised_d = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            de.good_mothers = de.good_mothers.drop(index=mothers_d)
            de.denoised_d = de.denoised_d.append(de.good_mothers, ignore_index=True)
            de.denoised_d = de.denoised_d.sort_values([de.count], axis=0, ascending=False)
        else:
            de.denoised_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            for mother in tqdm(mothers_d):
                row = [
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[
                        0] +
                    list(de.data_initial.loc[list(pd.Series(de.denoised_d_output) == mother), de.abund_col_names].sum(
                        0)) +
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
                row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
                de.denoised_d = de.denoised_d.append(row, ignore_index=True)
                de.good_mothers = de.good_mothers.drop(index=mother)
            de.denoised_d = de.denoised_d.append(de.good_mothers, ignore_index=True)
            de.denoised_d = de.denoised_d.sort_values([de.count], axis=0, ascending=False)
        if 'row' in locals():
            del row
        del de.good_mothers, mothers_d, de.denoised_d_output

    if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
        de.good_mothers = de.data_initial.loc[de.good_seq][de.first_col_names + de.abund_col_names + [de.seq]]
        print('writing output_ratio_d')
        # writing ratio_d
        if de.cores > 1:
            pool = mp.Pool(de.cores)
            [row] = zip(*pool.map(de.write_output_ratio_d, [mother for mother in mothers_ratio_d]))
            pool.close()
            del pool
            de.denoised_ratio_d = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            de.good_mothers = de.good_mothers.drop(index=mothers_ratio_d)
            de.denoised_ratio_d = de.denoised_ratio_d.append(de.good_mothers, ignore_index=True)
            de.denoised_ratio_d = de.denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
        else:
            de.denoised_ratio_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            for mother in tqdm(mothers_ratio_d):
                row = [
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[
                        0] +
                    list(de.data_initial.loc[
                        list(pd.Series(de.denoised_ratio_d_output) == mother), de.abund_col_names].sum(
                        0)) +
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
                row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
                de.denoised_ratio_d = de.denoised_ratio_d.append(row, ignore_index=True)
                de.good_mothers = de.good_mothers.drop(index=mother)
            de.denoised_ratio_d = de.denoised_ratio_d.append(de.good_mothers, ignore_index=True)
            de.denoised_ratio_d = de.denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
        if 'row' in locals():
            del row
        del de.good_mothers, mothers_ratio_d, de.denoised_ratio_d_output


def run_denoise_entropy(de):
    seq_length = []
    seq_length_per_read = []
    for i in list(range(de.data_initial.shape[0])):
        i_seq = de.data_initial.loc[i, de.seq]
        i_count = de.data_initial.loc[i, de.count]
        seq_length.append(len(i_seq))
        seq_length_per_read.append([len(i_seq)] * i_count)
    seq_length_per_read = list(itertools.chain.from_iterable(seq_length_per_read))

    uniq_seq_lengths = set()
    uniq_seq_lengths.update(seq_length)
    uniq_seq_lengths = list(uniq_seq_lengths)

    # separate data in different DataFrames by sequence length
    if len(de.modal_length_value) == 0:
        de.modal_length_value = modal_length(seq_length_per_read)

    if len(de.modal_length_value) != 1:
        for e in range(0, len(de.modal_length_value)):
            if ((de.modal_length_value[e] - 1) % 3) == 0:
                good_modal_length_value = de.modal_length_value[e]
                break
        if 'good_modal_length_value' not in locals():
            good_modal_length_value = de.modal_length_value[0]

        print('WARNING!! %s not available to run with Entropy. '
              'Equal number of seqs with different seq length' % de.MOTUfile)
        print('set -m as one value of the following: %s ' % de.modal_length_value)
        print('DnoisE will run with sequence length %s and its multiples' % good_modal_length_value)
    else:
        good_modal_length_value = de.modal_length_value

    allowed_lengths = np.array(uniq_seq_lengths) - good_modal_length_value
    allowed_lengths = list(allowed_lengths % 3 == 0)
    allowed_lengths = list(itertools.compress(uniq_seq_lengths, allowed_lengths))

    del good_modal_length_value, seq_length_per_read, de.modal_length_value

    de.output_info = pd.DataFrame()
    if (de.output_type == 'ratio') or (de.output_type == 'all'):
        de.denoised_ratio = pd.DataFrame()
    if (de.output_type == 'd') or (de.output_type == 'all'):
        de.denoised_d = pd.DataFrame()
    if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
        de.denoised_ratio_d = pd.DataFrame()

    for len_seq in allowed_lengths:

        desub = DnoisEFunctions()
        copy_to_subset(declass=de, desub=desub, seq_length=seq_length, len_seq=len_seq)

        if de.compute_entropy:
            if desub.initial_pos == 1:
                e1, e2, e3 = en.mean_entropy(desub.data_initial)
            if desub.initial_pos == 2:
                e2, e3, e1 = en.mean_entropy(desub.data_initial)
            if desub.initial_pos == 3:
                e3, e1, e2 = en.mean_entropy(desub.data_initial)
            desub.Ad1 = e1 / (e1 + e2 + e3)
            desub.Ad2 = e2 / (e1 + e2 + e3)
            desub.Ad3 = e3 / (e1 + e2 + e3)
            print('entropy values for length {:.0f} (first nt is a position {:.0f}):\n'
                  '\t {:.3f} for first position of codon\n'
                  '\t {:.3f} for second position of codon\n'
                  '\t {:.3f} for third position of codon'.format(len_seq, desub.initial_pos, e1, e2, e3))

        # maximum ratio allowed
        desub.max_ratio = (1 / 2) ** (desub.alpha * 1 * min(desub.Ad1, desub.Ad2, desub.Ad3) + 1)

        if desub.output_type == 'ratio':
            desub.denoised_ratio_output = [desub.data_initial.loc[0, 'id']]
        else:
            desub.denoised_d_output = [desub.data_initial.loc[0, 'id']]
            desub.denoised_ratio_output = [desub.data_initial.loc[0, 'id']]
            desub.denoised_ratio_d_output = [desub.data_initial.loc[0, 'id']]

        desub.output_info = [{'daughter': desub.data_initial.loc[0, 'id'], 'mother_d': None, 'd': None,
                              'mother_ratio': None, 'ratio': None,
                              'mother_ratio_d': None, 'xavier_criteria': None,
                              'difpos1': None, 'difpos2': None, 'difpos3': None}]
        desub.good_seq = [True]
        desub.abund_col_names.insert(0, de.count)
        desub.run_list = [{'id': desub.data_initial.loc[0, 'id'], de.count: desub.data_initial.loc[0, de.count],
                           'run': True, 'daughter': False}]

        run_dnoise_testing(desub)

        desub.output_info = pd.DataFrame.from_dict(desub.output_info)

        if (desub.output_type == 'ratio') or (desub.output_type == 'all'):
            mothers_ratio = desub.output_info.mother_ratio.unique()[1:]
        if (desub.output_type == 'd') or (desub.output_type == 'all'):
            mothers_d = desub.output_info.mother_d.unique()[1:]
        if (desub.output_type == 'ratio_d') or (desub.output_type == 'all'):
            mothers_ratio_d = desub.output_info.mother_ratio_d.unique()[1:]

        de.output_info = pd.concat([de.output_info, desub.output_info], ignore_index=True)

        del desub.output_info

        desub.data_initial = desub.data_initial.set_index(desub.data_initial.loc[:, 'id'])

        if (desub.output_type == 'ratio') or (desub.output_type == 'all'):
            desub.good_mothers = desub.data_initial.loc[desub.good_seq][de.first_col_names +
                                                                        desub.abund_col_names + [de.seq]]
            # writing ratio
            if desub.cores > 1:
                pool = mp.Pool(desub.cores)
                [row] = zip(*pool.map(desub.write_output_ratio, [mother for mother in mothers_ratio]))
                pool.close()
                del pool
                desub.denoised_ratio = pd.DataFrame(row,
                                                    columns=[de.first_col_names + desub.abund_col_names +
                                                             [de.seq]][0])
                desub.good_mothers = desub.good_mothers.drop(index=mothers_ratio)
                desub.denoised_ratio = desub.denoised_ratio.append(desub.good_mothers, ignore_index=True)
                desub.denoised_ratio = desub.denoised_ratio.sort_values([desub.count], axis=0, ascending=False)
            else:
                desub.denoised_ratio = pd.DataFrame(columns=[de.first_col_names + desub.abund_col_names + [de.seq]][0])
                for mother in tqdm(mothers_ratio):
                    row = [
                        desub.good_mothers[list(desub.good_mothers.loc[:, 'id'] == mother)][de.first_col_names
                        ].values.tolist()[0] +
                        list(desub.data_initial.loc[
                            list(pd.Series(desub.denoised_ratio_output) == mother), desub.abund_col_names].sum(
                            0)) +
                        desub.good_mothers[list(desub.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
                    row = pd.Series(row[0], index=[de.first_col_names + desub.abund_col_names + [de.seq]][0])
                    desub.denoised_ratio = desub.denoised_ratio.append(row, ignore_index=True)
                    desub.good_mothers = desub.good_mothers.drop(index=mother)
                desub.denoised_ratio = desub.denoised_ratio.append(desub.good_mothers, ignore_index=True)
                desub.denoised_ratio = desub.denoised_ratio.sort_values([de.count], axis=0, ascending=False)
            if 'row' in locals():
                del row
            de.denoised_ratio = pd.concat([de.denoised_ratio, desub.denoised_ratio], ignore_index=True)
            del desub.denoised_ratio, desub.good_mothers, mothers_ratio, desub.denoised_ratio_output

        if (desub.output_type == 'd') or (desub.output_type == 'all'):
            desub.good_mothers = desub.data_initial.loc[desub.good_seq][de.first_col_names +
                                                                        desub.abund_col_names + [de.seq]]
            # writing d
            if desub.cores > 1:
                pool = mp.Pool(desub.cores)
                [row] = zip(*pool.map(desub.write_output_d, [mother for mother in mothers_d]))
                pool.close()
                del pool
                desub.denoised_d = pd.DataFrame(row, columns=[de.first_col_names + desub.abund_col_names + [de.seq]][0])
                desub.good_mothers = desub.good_mothers.drop(index=mothers_d)
                desub.denoised_d = desub.denoised_d.append(desub.good_mothers, ignore_index=True)
                desub.denoised_d = desub.denoised_d.sort_values([de.count], axis=0, ascending=False)
            else:
                desub.denoised_d = pd.DataFrame(columns=[de.first_col_names + desub.abund_col_names + [de.seq]][0])
                for mother in tqdm(mothers_d):
                    row = [
                        desub.good_mothers[list(desub.good_mothers.loc[:, 'id'] == mother)][
                            de.first_col_names].values.tolist()[
                            0] +
                        list(desub.data_initial.loc[
                            list(pd.Series(desub.denoised_d_output) == mother), desub.abund_col_names].sum(
                            0)) +
                        desub.good_mothers[list(desub.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
                    row = pd.Series(row[0], index=[de.first_col_names + desub.abund_col_names + [de.seq]][0])
                    desub.denoised_d = desub.denoised_d.append(row, ignore_index=True)
                    desub.good_mothers = desub.good_mothers.drop(index=mother)
                desub.denoised_d = desub.denoised_d.append(desub.good_mothers, ignore_index=True)
                desub.denoised_d = desub.denoised_d.sort_values([de.count], axis=0, ascending=False)
            if 'row' in locals():
                del row
            de.denoised_d = pd.concat([de.denoised_d, desub.denoised_d], ignore_index=True)
            del desub.denoised_d, desub.good_mothers, mothers_d, desub.denoised_d_output

        if (desub.output_type == 'ratio_d') or (desub.output_type == 'all'):
            desub.good_mothers = desub.data_initial.loc[desub.good_seq][
                de.first_col_names + desub.abund_col_names + [de.seq]]
            # writing ratio_d
            if desub.cores > 1:
                pool = mp.Pool(desub.cores)
                [row] = zip(*pool.map(desub.write_output_ratio_d, [mother for mother in mothers_ratio_d]))
                pool.close()
                del pool
                desub.denoised_ratio_d = pd.DataFrame(row,
                                                      columns=[de.first_col_names + desub.abund_col_names + [de.seq]][
                                                          0])
                desub.good_mothers = desub.good_mothers.drop(index=mothers_ratio_d)
                desub.denoised_ratio_d = desub.denoised_ratio_d.append(desub.good_mothers, ignore_index=True)
                desub.denoised_ratio_d = desub.denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
            else:
                desub.denoised_ratio_d = pd.DataFrame(
                    columns=[de.first_col_names + desub.abund_col_names + [de.seq]][0])
                for mother in tqdm(mothers_ratio_d):
                    row = [
                        desub.good_mothers[list(desub.good_mothers.loc[:, 'id'] == mother)][
                            de.first_col_names].values.tolist()[
                            0] +
                        list(desub.data_initial.loc[
                            list(pd.Series(desub.denoised_ratio_d_output) == mother), desub.abund_col_names].sum(
                            0)) +
                        desub.good_mothers[list(desub.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
                    row = pd.Series(row[0], index=[de.first_col_names + desub.abund_col_names + [de.seq]][0])
                    desub.denoised_ratio_d = desub.denoised_ratio_d.append(row, ignore_index=True)
                    desub.good_mothers = desub.good_mothers.drop(index=mother)
                desub.denoised_ratio_d = desub.denoised_ratio_d.append(desub.good_mothers, ignore_index=True)
                desub.denoised_ratio_d = desub.denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
            if 'row' in locals():
                del row

            de.denoised_ratio_d = pd.concat([de.denoised_ratio_d, desub.denoised_ratio_d], ignore_index=True)
            del desub.denoised_ratio_d, desub.good_mothers, mothers_ratio_d, desub.denoised_ratio_d_output

        del desub

    de.output_info.to_csv(str(de.MOTUoutfile + '_denoising_info.csv'), index=False)
    del de.output_info


def run_from_info(de):
    de.data_initial.index = de.data_initial.id
    de.abund_col_names.insert(0, de.count)
    if (de.output_type == 'ratio') or (de.output_type == 'all'):
        mothers_ratio = de.merge_data.mother_ratio.unique()[1:]
        de.good_mothers = de.data_initial.loc[de.merge_data.daughter[de.merge_data.mother_d.isna()]][de.first_col_names +
                                                                                                     de.abund_col_names +
                                                                                                     [de.seq]]
        print('writing output_ratio')
        # writing ratio
        if de.cores > 1:
            pool = mp.Pool(de.cores)
            [row] = zip(*pool.map(de.write_ratio_from_info, [mother for mother in mothers_ratio]))
            pool.close()
            del pool
            de.denoised_ratio = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            de.good_mothers = de.good_mothers.drop(index=mothers_ratio)
            de.denoised_ratio = de.denoised_ratio.append(de.good_mothers, ignore_index=True)
            de.denoised_ratio = de.denoised_ratio.sort_values([de.count], axis=0, ascending=False)
        else:
            de.denoised_ratio = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            for mother in tqdm(mothers_ratio):
                row = [
                    de.good_mothers[list(de.good_mothers.id == mother)][de.first_col_names].values.tolist()[
                        0] +
                    list(de.data_initial.loc[[[mother] +
                                              list(de.merge_data.daughter[de.merge_data.mother_ratio == mother])][0],
                                             de.abund_col_names].sum(0)) +
                    de.good_mothers[list(de.good_mothers.id == mother)][de.seq].values.tolist()]
                row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
                de.denoised_ratio = de.denoised_ratio.append(row, ignore_index=True)
                de.good_mothers = de.good_mothers.drop(index=mother)
            de.denoised_ratio = de.denoised_ratio.append(de.good_mothers, ignore_index=True)
            de.denoised_ratio = de.denoised_ratio.sort_values([de.count], axis=0, ascending=False)
        if 'row' in locals():
            del row
        del de.good_mothers, mothers_ratio

    if (de.output_type == 'd') or (de.output_type == 'all'):
        mothers_d = de.merge_data.mother_d.unique()[1:]
        de.good_mothers = de.data_initial.loc[de.merge_data.daughter[de.merge_data.mother_d.isna()]
        ][de.first_col_names + de.abund_col_names + [de.seq]]
        print('writing output_d')
        # writing ratio
        if de.cores > 1:
            pool = mp.Pool(de.cores)
            [row] = zip(*pool.map(de.write_d_from_info, [mother for mother in mothers_d]))
            pool.close()
            del pool
            de.denoised_d = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            de.good_mothers = de.good_mothers.drop(index=mothers_d)
            de.denoised_d = de.denoised_d.append(de.good_mothers, ignore_index=True)
            de.denoised_d = de.denoised_d.sort_values([de.count], axis=0, ascending=False)
        else:
            de.denoised_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            for mother in tqdm(mothers_d):
                row = [
                    de.good_mothers[list(de.good_mothers.id == mother)][de.first_col_names].values.tolist()[
                        0] +
                    list(de.data_initial.loc[[[mother] +
                                              list(de.merge_data.daughter[de.merge_data.mother_d == mother])][0],
                                             de.abund_col_names].sum(0)) +
                    de.good_mothers[list(de.good_mothers.id == mother)][de.seq].values.tolist()]
                row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
                de.denoised_d = de.denoised_d.append(row, ignore_index=True)
                de.good_mothers = de.good_mothers.drop(index=mother)
            de.denoised_d = de.denoised_d.append(de.good_mothers, ignore_index=True)
            de.denoised_d = de.denoised_d.sort_values([de.count], axis=0, ascending=False)
        if 'row' in locals():
            del row
        del de.good_mothers, mothers_d
    if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
        mothers_ratio_d = de.merge_data.mother_ratio_d.unique()[1:]
        de.good_mothers = de.data_initial.loc[de.merge_data.daughter[de.merge_data.mother_d.isna()]
        ][de.first_col_names + de.abund_col_names + [de.seq]]
        print('writing output_ratio_d')
        # writing ratio
        if de.cores > 1:
            pool = mp.Pool(de.cores)
            [row] = zip(*pool.map(de.write_ratio_d_from_info, [mother for mother in mothers_ratio_d]))
            pool.close()
            del pool
            de.denoised_ratio_d = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            de.good_mothers = de.good_mothers.drop(index=mothers_ratio_d)
            de.denoised_ratio_d = de.denoised_ratio_d.append(de.good_mothers, ignore_index=True)
            de.denoised_ratio_d = de.denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
        else:
            de.denoised_ratio_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            for mother in tqdm(mothers_ratio_d):
                row = [
                    de.good_mothers[list(de.good_mothers.id == mother)][de.first_col_names].values.tolist()[
                        0] +
                    list(de.data_initial.loc[[[mother] +
                                              list(de.merge_data.daughter[de.merge_data.mother_ratio_d == mother])][0],
                                             de.abund_col_names].sum(0)) +
                    de.good_mothers[list(de.good_mothers.id == mother)][de.seq].values.tolist()]
                row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
                de.denoised_ratio_d = de.denoised_ratio_d.append(row, ignore_index=True)
                de.good_mothers = de.good_mothers.drop(index=mother)
            de.denoised_ratio_d = de.denoised_ratio_d.append(de.good_mothers, ignore_index=True)
            de.denoised_ratio_d = de.denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
        if 'row' in locals():
            del row
        del de.good_mothers, mothers_ratio_d

    del de.merge_data
