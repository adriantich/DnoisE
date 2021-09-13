#!/usr/bin/env python3

"""
.. codeauthor:: Adrià Antich <adriantich@gmail.com>

This programme is called by the DnoisE.

transform_data.py transforms data imported by import_data.py script as pandas.DataFrame. Data is filtered and sorted by
read abundance. Also other variables necessary for following steps od DnoisE are defined.
This script also computes the entropy if specified.

"""

import itertools
import numpy as np
import entropy as en


def transform_data(de):
    print('transforming data')
    if de.part == 1:

        de.data_initial[de.seq] = de.data_initial[de.seq].str.upper()

        # obtain a column with total reads per seq.
        if not de.justcount:
            de.abund_col_names = list(de.data_initial.columns)[(de.start - 1):de.end]
            de.first_col_names = list(de.data_initial.columns)[0:(de.start - 2)]
            if de.count in de.first_col_names:
                de.first_col_names.remove(de.count)
            if de.seq in de.first_col_names:
                de.first_col_names.remove(de.seq)
            de.data_initial.loc[:, de.count] = de.data_initial.loc[:, de.abund_col_names].sum(axis=1)
        else:
            de.first_col_names = ['id']

        # sort by total reads
        de.data_initial = de.data_initial.sort_values([de.count], axis=0, ascending=False)

        # delete seqs with 0 reads
        de.data_initial = de.data_initial.loc[(de.data_initial.loc[:, de.count] != 0)]

        # maximum ratio allowed
        de.max_ratio = (1 / 2) ** (de.alpha * 1 * min(de.Ad1, de.Ad2, de.Ad3) + 1)
        # reorder index
        de.data_initial.index = list(range(de.data_initial.shape[0]))

    else:
        print('no data transformation')


def transform_data_entropy_correction(de):
    print('transforming data')
    if de.part == 1:

        de.data_initial[de.seq] = de.data_initial[de.seq].str.upper()

        # obtain a column with total reads per seq.
        if not de.justcount:
            de.abund_col_names = list(de.data_initial.columns)[(de.start - 1):de.end]
            de.first_col_names = list(de.data_initial.columns)[0:(de.start - 2)]
            if de.count in de.first_col_names:
                de.first_col_names.remove(de.count)
            if de.seq in de.first_col_names:
                de.first_col_names.remove(de.seq)
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
                seq_length_per_read.append([len(i_seq)] * i_count)
            seq_length_per_read = list(itertools.chain.from_iterable(seq_length_per_read))

            uniq_seq_lengths = set()
            for x in seq_length:
                uniq_seq_lengths.add(x)
            print(uniq_seq_lengths)

            # separate data in different DataFrames by sequence length
            if len(de.modal_length_value) == 0:
                de.modal_length_value = de.modal_length(seq_length_per_read)

            if len(de.modal_length_value) == 1:
                de.data_initial = de.data_initial.loc[(np.asarray(seq_length) == de.modal_length_value)]
            else:
                for e in range(0, len(de.modal_length_value)):
                    if ((de.modal_length_value[e] - 1) % 3) == 0:
                        good_modal_length_value = de.modal_length_value[e]
                        break
                if 'good_modal_length_value' not in locals():
                    good_modal_length_value = de.modal_length_value[0]

                print('WARNING!! %s not available to run with Entropy. '
                      'Equal number of seqs with different seq length' % de.MOTUfile)
                print('set -m as one value of the following: %s ' % de.modal_length_value)
                print('DnoisE will run with sequence length %s' % good_modal_length_value)
                de.data_initial = de.data_initial.loc[(np.asarray(seq_length) == good_modal_length_value)]

                del good_modal_length_value

            del seq_length, seq_length_per_read, de.modal_length_value

            if de.compute_entropy:
                if de.initial_pos == 1:
                    e1, e2, e3 = en.mean_entropy(de.data_initial)
                if de.initial_pos == 2:
                    e2, e3, e1 = en.mean_entropy(de.data_initial)
                if de.initial_pos == 3:
                    e3, e1, e2 = en.mean_entropy(de.data_initial)
                de.Ad1 = e1 / (e1 + e2 + e3)
                de.Ad2 = e2 / (e1 + e2 + e3)
                de.Ad3 = e3 / (e1 + e2 + e3)
            print('entropy values (first nt is a position {:.0f}:\n'
                  '\t {:.3f} for first position of codon\n'
                  '\t {:.3f} for second position of codon\n'
                  '\t {:.3f} for third position of codon'.format(de.initial_pos, e1, e2, e3))

            # maximum ratio allowed
            de.max_ratio = (1 / 2) ** (de.alpha * 1 * min(de.Ad1, de.Ad2, de.Ad3) + 1)
            # reorder index
            de.data_initial.index = list(range(de.data_initial.shape[0]))

    else:
        print('no data transformation')
