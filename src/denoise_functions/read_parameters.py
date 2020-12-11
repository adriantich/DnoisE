
"import parameters and return variables"

argument_list = ['-i', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/altres/versions_millorant/PHY1_000372023', '-o', '/home/adriantich/Nextcloud/1_tesi_Adrià/Denoise/altres/versions_millorant/PHY1_000372023', '-f', 'False', '-c', '2', '-n', 'count', '-a', '5', '-q', 'sequence', '-p', '1', '-e', '0.4298,0.1833,0.9256', '-y', 'T']


import getopt
import sys

class read_parameters:



    def __init__(self):
        print("reading variables")
        print(argument_list)

short_options = "hi:o:P:f:F:c:s:z:n:a:q:p:e:y:"
long_options = ["help", "input=", "output=", "part=", "fasta_input=", "fasta_output=", "cores=",
                    "start_sample_cols=",
                    "end_sample_cols=", "count_name=", "alpha=", "sequence=", "separation=", "entropy=",
                    "entropy_influence="]
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
            de.fasta_output = True
        else:
            de.fasta_output = False
        arg_F = True
        print("Fasta output file: %s" % de.fasta_output)
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
            if de.fasta_output:
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
