#!/usr/bin/env python3

"""
.. codeauthor:: Adrià Antich <adriantich@gmail.com>

This programme is called by the DnoisE.

within_MOTU.py contains the code to run DnoisE within each MOTU.

"""

import sys
import multiprocessing as mp
from dnoise.denoise_functions import *
from dnoise.running_denoise import *
from dnoise.transform_data import *
from dnoise.write_output import *
import copy
import pandas as pd


def run_within_MOTU(de):
    sort_MOTU(de)
    for motu in list(de.data_initial[de.motu_column].unique()):
        demotu = copy.deepcopy(de)
        demotu.data_initial = demotu.data_initial.loc[(np.asarray(list(de.data_initial[de.motu_column])) == motu)]
        demotu.data_initial.index = list(range(demotu.data_initial.shape[0]))
        try:
            del demotu.output_info, demotu.denoised_ratio, demotu.denoised_d, demotu.denoised_ratio_d
        except:
            pass

        transform_data(demotu)
        if de.entropy:
            run_denoise_entropy(demotu)
        else:
            run_denoise(demotu)
        if list(de.data_initial[de.motu_column].unique())[0] == motu:
            de.output_info = demotu.output_info
            if (de.output_type == 'ratio') or (de.output_type == 'all'):
                de.denoised_ratio = demotu.denoised_ratio
            if (de.output_type == 'd') or (de.output_type == 'all'):
                de.denoised_d = demotu.denoised_d
            if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
                de.denoised_ratio_d = demotu.denoised_ratio_d
        else:
            de.output_info = pd.concat([de.output_info, demotu.output_info], ignore_index=True)
            if (de.output_type == 'ratio') or (de.output_type == 'all'):
                de.denoised_ratio = pd.concat([de.denoised_ratio, demotu.denoised_ratio], ignore_index=True)
            if (de.output_type == 'd') or (de.output_type == 'all'):
                de.denoised_d = pd.concat([de.denoised_d, demotu.denoised_d], ignore_index=True)
            if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
                de.denoised_ratio_d = pd.concat([de.denoised_ratio_d, demotu.denoised_ratio_d], ignore_index=True)


