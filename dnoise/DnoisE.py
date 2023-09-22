#!/usr/bin/env python3

"""
.. codeauthor:: Adrià Antich <adriantich@gmail.com>

DnoisE is designed to denoise data sets from illumina sequencing using (optionally) a correction of the distance
measure (parameter d) according to the entropy of each codon position.
"""

import sys
import multiprocessing as mp
from dnoise.denoise_functions import *
from dnoise.import_data import *
from dnoise.running_denoise import *
from dnoise.transform_data import *
from dnoise.write_output import *
from dnoise.get_entropy import *
from dnoise.within_MOTU import *


def main():
    # if platform.system() == 'Linux':
    #     mp.set_start_method('fork')
    # else:
    #     print('not Linux system detected')
    #     mp.set_start_method('spawn')

    de = DnoisEFunctions()

    # Get full command-line arguments
    argument_list = sys.argv

    # Keep all but the first
    argument_list = argument_list[1:]

    print(argument_list)
    de.read_parameters(argument_list)
    import_data(de)
    if de.within_motu:
        run_within_MOTU(de)
    else:
        transform_data(de)
        if de.get_entropy:
            get_entropy_func(de)
            # this will break to return a file with entropy values
        if not de.merge_from_info:
            if de.entropy:
                run_denoise_entropy(de)
            else:
                run_denoise(de)
        else:
            run_from_info(de)

    write_output(de)

    print('done')


if __name__ == '__main__':
    main()
