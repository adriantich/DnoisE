#import all the modules that will be used
import pandas as pd
import numpy as np
import multiprocessing as mp
import csv
import getopt, sys
from tqdm import tqdm
import Bio
import re
import stats
from denoise_functions import denoise_functions


de = denoise_functions()

# Get full command-line arguments
full_cmd_arguments = sys.argv

# Keep all but the first
argument_list = full_cmd_arguments[1:]
# per part 1
# argument_list = ['-i', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/altres/versions_millorant/PHY1_001014928', '-o', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/altres/versions_millorant/PHY1_001014928', '-P', '1', '-f', 'False', '-c', '2', '-s', '4', '-z', '79', '-n', 'count', '-a', '5', '-q', 'sequence', '-p', '1', '-e', '0.4298,0.1833,0.9256', '-y', 'T']
# per part 2
# argument_list = ['-i', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/altres/versions_millorant/PHY1_001014928', '-o', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/altres/versions_millorant/PHY1_001014928', '-P', '2']
# sense parts
argument_list = ['-i', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/altres/versions_millorant/PHY1_000372023', '-o', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/altres/versions_millorant/PHY1_000372023', '-f', 'False', '-c', '2', '-n', 'count', '-a', '5', '-q', 'sequence', '-p', '1', '-e', '0.4298,0.1833,0.9256', '-y', 'T']


print(argument_list)
short_options = "hi:o:P:f:F:c:s:z:n:a:q:p:e:y:"
long_options = ["help", "input=", "output=", "part=", "fasta_input=", "fasta_output=", "cores=", "start_sample_cols=",
                "end_sample_cols=", "count_name=", "alpha=", "sequence=", "separation=", "entropy=",  "entropy_influence="]
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
              " -f --fasta_input logical, if T, fasta file as input, if F (default) .csv as input\n"
              " -F --fasta_output logical, if T, fasta file as input, if F (default) .csv as input\n"
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
        de.MOTUfile = current_value
        arg_i = True
    elif current_argument in ("-o", "--output"):
        print("Output files will be named %s*" % current_value)
        de.MOTUoutfile = current_value
        arg_o = True
    elif current_argument in ("-P", "--Part"):
        print("Part %s running" % current_value)
        part = int(current_value)
        parted = True
        arg_P = True
    elif current_argument in ("-c", "--cores"):
        print("Running with %s cores" % current_value)
        de.cores = int(current_value)
        arg_c = True
    elif current_argument in ("-f", "--fasta_input"):
        if current_value == "T":
            fasta = True
        else:
            fasta = False
        arg_f = True
        print("Fasta file: %s" % fasta)
    elif current_argument in ("-F", "--fasta_output"):
        if current_value == "T":
            fasta_output = True
        else:
            fasta_output = False
        arg_F = True
        print("Fasta output file: %s" % fasta_output)
    elif current_argument in ("-s", "--start_sample_cols"):
        start = int(current_value)
        arg_s = True
        print("Abundant cols starts in: %s" % current_value)
    elif current_argument in ("-z", "--end_sample_cols"):
        end = int(current_value)
        arg_z = True
        print("Abundant cols ends in: %s" % current_value)
    elif current_argument in ("-n", "--count_name"):
        de.count = current_value
        arg_n = True
        print("Count col name: %s" % current_value)
    elif current_argument in ("-a", "--alpha"):
        de.alpha = int(current_value)
        arg_a = True
        print("Alpha: %s" % current_value)
    elif current_argument in ("-q", "--sequence"):
        de.seq = current_value
        arg_q = True
        print("Seq: %s" % current_value)
    elif current_argument in ("-p", "--sep"):
        print(current_value)
        if current_value == "1":
            sep = '\t'
        if current_value == "2":
            sep = ','
        if current_value == "3":
            sep = ';'
        arg_p = True
        print("Sep: %s" % sep)
    elif current_argument in ("-e", "--entropy"):
        E1, E2, E3 = current_value.split(",")
        print(str("E1: " + E1 + " \nE2: " + E2 + " \nE3: " + E3))
        E1 = float(E1)
        E2 = float(E2)
        E3 = float(E3)
        arg_e = True
    elif current_argument in ("-y", "--entropy_influence"):
        if current_value == "T":
            de.entropy = True
        else:
            de.entropy = False
        arg_y = True
        print("Is entropy taken into account: %s" % de.entropy)

if 'arg_i' not in locals():
    print("Err: input file needed")
    exit()
if 'arg_o' not in locals():
    print("Err: output path needed")
    exit()
if 'arg_P' not in locals():
    print("Running all at once")
    parted = False
    part = 1
if part == 1:
    if 'arg_c' not in locals():
        print("cores not given, 1 core by default")
        de.cores = 1
    if 'arg_f' not in locals():
        # de moment per a la opcio amb fasta no hi ha la opció de posar samples
        print("by default, fasta file expected")
        fasta = True
        de.count = "count"
        de.justcount = True
        if 'arg_F' not in locals():
            de.fasta_output = True
            print("by default, fasta output")
    else:
        if "arg_n" not in locals():
            print("count_name not given, 'count' by default")
            de.count = "count"
        if "arg_s" not in locals():
            if "arg_z" not in locals():
                print("start and end not given, no samples in file")
                abund_col_names = [de.count]
                de.justcount = True
            else:
                print("Err: end given but no start")
                exit()
        elif "arg_z" not in locals():
            print("Err: start given but no end")
            exit()
        else:
            de.justcount = False
            if fasta_output:
                print('WARNING!!!! .csv output when samples information is given')
            de.fasta_output = False
            print("output as .csv")
        if "arg_q" not in locals():
            print("sequence not given, 'sequence' by default")
            de.seq = 'sequence'
        if "arg_p" not in locals():
            print("Separation not given, '\t' by default")
            sep = '\t'
    if "arg_a" not in locals():
        print("alpha not given, 5 by default")
        de.alpha = 5
    if "arg_y" not in locals():
        de.entropy = False
        print("Ad correction not applied")
    else:
        if "arg_e" not in locals():
            E1, E2, E3 = (0.4298, 0.1833, 0.9256)
            print("Entropy set as 0.4298, 0.1833 and 0.9256 by default")







#import motu_tab file
#
# import pandas as pd
# import numpy as np
# import Levenshtein as lv
# import multiprocessing as mp
# import csv
# import getopt, sys
# from tqdm import tqdm
# import Bio
# import re
# import stats
# from denoise_functions import denoise_functions
#
# de=denoise_functions()
#
# # objects from de
# de.MOTUfile = '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/versions_millorant/PHY1_000161710'
# de.MOTUoutfile = '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/versions_millorant/PHY1_000161710'
#
# de.cores = 1
# de.count = "count"
# de.alpha = 5
# de.seq = 'sequence'
# de.entropy = True
#
# fasta = False
# start = 4
# end = 67
# de.justcount = False
# sep = '\t'
# E1, E2, E3 = (0.4298, 0.1833, 0.9256)
# part = 1
# parted = True


if part != 3:
    # prepearing data

    if part == 1:
        print('reading input file')
        if fasta:
            # data_inicial = pd.DataFrame()
            for fastaseq in Bio.SeqIO.parse("fasta.fasta", "fasta"):
                seq_seq = fastaseq.seq._data
                seq_id = re.findall('(.+); ', fastaseq.description)[0]
                seq_count = re.findall('size=(.+);', fastaseq.description)[0]
                de.data_inicial = pd.concat(
                    [de.data_inicial, pd.DataFrame({"id": seq_id, "count": int(seq_count), "sequence": [seq_seq]})])
        else:
            de.data_inicial = pd.read_csv(de.MOTUfile, sep=sep)

        print('input file read')

        # Afegim el factor corrector Ad. Aquest valor corregeix la d segons cada posició tenint en compte que cada posició té una
        # entropia diferent. De manera natural la posició 2 té menys entropia que la 1 i aquesta menys que la 3.
        # Hem optat per calcular aquest factor de manera que Adi=Ei/sum(Ej) * 3. Adi es la Ad a la posició i, Ei es la entropia en
        # aquesta posició i sum(Ej) és la suma de les entropies a les 3 posicions. Es multiplica per 3 de manera que així les Ad no tenen
        # un efecte conjunt en la funció de Edgar

        if de.entropy:
            de.Ad1 = E1 / (E1 + E2 + E3) * 3
            de.Ad2 = E2 / (E1 + E2 + E3) * 3
            de.Ad3 = E3 / (E1 + E2 + E3) * 3
        else:
            de.Ad1, de.Ad2, de.Ad3 = (1, 1, 1)

        # define alpha and minimum min_mother reads

        de.min_mother = 1 / (1 / 2) ** (de.alpha * min(de.Ad1, de.Ad2, de.Ad3) + 1)
        print('min_mother equals to %s' % de.min_mother)
        print('and Ad corr:')
        print(de.Ad1, de.Ad2, de.Ad3)
        # the program should work with different seq_length, if not, filter and get de mode

        # obtain a column with total reads per seq.
        if not de.justcount:
            de.abund_col_names = list(de.data_inicial.columns)[(start - 1):end]
            de.first_col_names = list(de.data_inicial.columns)[0:(start - 2)]
            de.data_inicial.loc[:, de.count] = de.data_inicial.loc[:, de.abund_col_names].sum(axis=1)
        else:
            de.first_col_names = ['id']

            # sort by total reads
        de.data_inicial = de.data_inicial.sort_values([de.count], axis=0, ascending=False)

        # delete seqs with 0 reads
        de.data_inicial = de.data_inicial.loc[(de.data_inicial.loc[:, de.count] != 0)]
        if de.entropy:
            # remove seq out of mode length
            seq_length = []
            for i in de.data_inicial.loc[:, de.seq]:
                seq_length.append(len(i))
            try:
                de.data_inicial = de.data_inicial.loc[(np.asarray(seq_length) == stats.mode(seq_length))]
            except:
                print('MOTU no available to run with Entropy. Equal number of seqs with different seq length')
                exit()
            del seq_length

            # reorder index
        de.data_inicial.index = list(range(de.data_inicial.shape[0]))

        # analising each seq. two main groups.
        # A. possible mothers (more than 64 reads i think)
        # B. less than 64 reads. cannot be mothers

        # create a runnedlist (rl) for A seqs with two columns

        de.runned_list = de.data_inicial.loc[(de.data_inicial.loc[:, de.count] >= de.min_mother), ['id', de.count]]

        # runned ---> FALSE
        de.runned_list['runned'] = np.repeat(False, de.runned_list.shape[0])

        # daughter --> FALSE
        de.runned_list['daughter'] = np.repeat(False, de.runned_list.shape[0])

        de.possibleMothers = de.runned_list.shape[0]
        de.not_possibleMothers = de.data_inicial.shape[0]

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

        # mothers_d = data_inicial.loc[(data_inicial.loc[:, count] >= min_mother), :]
        # mothers_ratio = data_inicial.loc[(data_inicial.loc[:, count] >= min_mother), :]
        # mothers_ratio_d = data_inicial.loc[(data_inicial.loc[:, count] >= min_mother), :]

        de.abund_col_names.insert(0, de.count)

        ##############################################################################################

        # for each seq testing (possible Daughter, pD)

        # here is presented 2 different situations of similar algorithm but changing the final testing
        # there are two main ways to do it. test firsth whether a seq is A. or B. and run all the following
        # programm or test whether a seq is A. or B. to run the last part of the code

    if not parted:
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

                # pool.map(denoising, [pos for pos in range(data_inicial.shape[0])])
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

                # pool.map(denoising, [pos for pos in range(data_inicial.shape[0])])
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
    if parted:
        if part == 1:
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

                    # pool.map(denoising, [pos for pos in range(data_inicial.shape[0])])
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

                    # pool.map(denoising, [pos for pos in range(data_inicial.shape[0])])
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

if part == 3:
    de.read_variables2()

if de.entropy:
    de.MOTUoutfile = str(de.MOTUoutfile + '_Adcorr')

print('writing output_info')
fieldnames = de.output_info[0].keys()
with open(str(de.MOTUoutfile + '_denoising_info.csv'), 'w') as output_file:
    dict_writer = csv.DictWriter(output_file, fieldnames=fieldnames)
    dict_writer.writeheader()
    dict_writer.writerows(list(de.output_info + list(de.output_info_2)))
if de.output_info_2 == []:
    del(dict_writer, de.output_info, fieldnames)
else:
    del (dict_writer, de.output_info, de.output_info_2, fieldnames)


de.data_inicial = de.data_inicial.set_index(de.data_inicial.loc[:, 'id'])

good_mothers = de.data_inicial.loc[de.good_seq + [False] * len(de.good_seq_2)][de.first_col_names + de.abund_col_names + [de.seq]]
good_daughters = de.data_inicial.loc[[False] * len(de.good_seq) + list(de.good_seq_2)][de.first_col_names + de.abund_col_names + [de.seq]]


print('writing output_d')
# writing d
denoised_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])

for mother in good_mothers.loc[:, 'id']:
    row = [good_mothers[list(good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[0] +
           list(de.data_inicial.loc[[de.mother_id(x, mother) for x in list(
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
           list(de.data_inicial.loc[[de.mother_id(x, mother) for x in list(
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
           list(de.data_inicial.loc[[de.mother_id(x, mother) for x in list(
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