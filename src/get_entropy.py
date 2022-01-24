import pandas as pd

from denoise_functions import *
import entropy as en


def get_entropy_func(de):
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
                good_modal_length_value = [de.modal_length_value[e]]
                break
        if 'good_modal_length_value' not in locals():
            good_modal_length_value = de.modal_length_value[0]

        print('WARNING!! %s not available to run with entropy correction. '
              'Equal number of seqs with different seq length' % de.MOTUfile)
        print('set -m as one value of the following: %s ' % de.modal_length_value)
        print('DnoisE will run with sequence length %s and its accepted variations (multiples of 3 '
              'nucleotides)' % good_modal_length_value)
    else:
        good_modal_length_value = de.modal_length_value

    if de.unique_length:
        allowed_lengths = good_modal_length_value
    else:
        allowed_lengths = np.array(uniq_seq_lengths) - good_modal_length_value
        allowed_lengths = list(allowed_lengths % 3 == 0)
        allowed_lengths = list(itertools.compress(uniq_seq_lengths, allowed_lengths))
        allowed_lengths.remove(good_modal_length_value[0])
        allowed_lengths.insert(0, good_modal_length_value[0])

    del seq_length_per_read, de.modal_length_value
    entropy_df = pd.DataFrame(columns=['seq_length', 'total_count', 'e1', 'e2', 'e3'])

    for i in list(range(len(allowed_lengths))):
        len_seq = allowed_lengths[i]

        desub = DnoisEFunctions()
        copy_to_subset(declass=de, desub=desub, seq_length=seq_length, len_seq=len_seq)

        if desub.initial_pos == 1:
            e1, e2, e3 = en.mean_entropy(desub.data_initial, de.seq, de.count)
        if desub.initial_pos == 2:
            e2, e3, e1 = en.mean_entropy(desub.data_initial, de.seq, de.count)
        if desub.initial_pos == 3:
            e3, e1, e2 = en.mean_entropy(desub.data_initial, de.seq, de.count)
        entr = [{'seq_length': len_seq, 'total_count': sum(desub.data_initial[de.count]), 'e1': e1, 'e2': e2, 'e3': e3}]
        entropy_df = entropy_df.append(pd.DataFrame(entr), ignore_index=True)
        del desub

    entropy_df.to_csv(de.MOTUoutfile, index=False)
    sys.exit()

