#!/usr/bin/env python3

"""
.. codeauthor:: Adrià Antich <adriantich@gmail.com>

DnoisE is designed to denoise data sets from illumina using parameter d (distance) corrected (optionally)
according to the entropy of each codon position.
"""

import sys
from denoise_functions import *
from import_data import *
from transform_data import *
from running_denoise import *
from write_output import *


if __name__ == '__main__':
    de = DnoisEFunctions()

    # Get full command-line arguments
    argument_list = sys.argv

    # Keep all but the first
    argument_list = argument_list[1:]

    print(argument_list)
    de.read_parameters(argument_list)
    import_data(de)
    transform_data(de)
    if not de.merge_from_info:
        if de.entropy:
            run_denoise_entropy(de)
        else:
            run_denoise(de)

    write_output(de)

    print('done')



