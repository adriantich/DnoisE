# import all the modules that will be used
import pandas as pd
import numpy as np
import multiprocessing as mp
import sys
from tqdm import tqdm
import itertools
from denoise_functions import denoise_functions
from transform_data import *
from running_denoise import run_denoise
from write_output import write_ouput

de = denoise_functions()

# Get full command-line arguments
argument_list = sys.argv

# Keep all but the first
argument_list = argument_list[1:]

print(argument_list)
de.read_parameters(argument_list)
import_data(de)
transform_data(de)
run_denoise(de)
write_ouput(de)


print('done')
