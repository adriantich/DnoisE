import pandas as pd
import itertools
import numpy as np
import entropy as en


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

def transform_data(de):
    print('start')
    if de.part == 1:

        de.data_initial[de.seq] = de.data_initial[de.seq].str.upper()

        print('input file read')

        # the program should work with different seq_length, if not, filter and get de mode

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

                print(
                    'WARNING!! %s no available to run with Entropy. Equal number of seqs with different seq length' % de.MOTUfile)
                print('set -m as one value of the following: %s ' % de.modal_length_value)
                print('DnoisE will run with sequence length %s' % good_modal_length_value)
                de.data_initial = de.data_initial.loc[(np.asarray(seq_length) == good_modal_length_value)]

                del good_modal_length_value

            del seq_length, seq_length_per_read, de.modal_length_value

        # reorder index
        de.data_initial.index = list(range(de.data_initial.shape[0]))
        if de.compute_entropy:
            if de.initial_pos == 1:
                de.Ad1, de.Ad2, de.Ad3 = en.mean_entropy(de.data_initial)
            if de.initial_pos == 2:
                de.Ad2, de.Ad3, de.Ad1 = en.mean_entropy(de.data_initial)
            if de.initial_pos == 3:
                de.Ad3, de.Ad1, de.Ad2 = en.mean_entropy(de.data_initial)
            print('entropy values:\n'
                  '\t {:.3f} for first position of codon\n'
                  '\t {:.3f} for second position of codon\n'
                  '\t {:.3f} for third position of codon'.format(de.Ad1, de.Ad2, de.Ad3))
            # maximum ratio allowed
        de.max_ratio = (1 / 2) ** (de.alpha * 1 * min(de.Ad1, de.Ad2, de.Ad3) + 1)


    else:
        print('no data transformation')

