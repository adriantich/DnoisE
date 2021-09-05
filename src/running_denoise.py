
import multiprocessing as mp
from tqdm import tqdm

def run_denoise(de):
    # first position is always the same
    if de.part == 1:
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
                                                                                                     range(
                                                                                                         len(de.runned_list),
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
                                                                                               range(
                                                                                                   len(de.runned_list),
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
        if de.json_database:
            de.write_variables()
    if de.part == 2:
        de.read_variables2()

