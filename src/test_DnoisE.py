
import unittest
from denoise_functions import *
from running_denoise import *
from entropy import *
import json
from copy import deepcopy


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)


def create_de_for_test(desub):
    with open('../test-DnoisE/test.json', 'r') as doc_json:
        dictionary = json.load(doc_json)
    de = Struct(**dictionary)
    desub.alpha = de.alpha
    desub.count = de.count
    desub.output_type = de.output_type
    desub.cores = de.cores
    desub.initial_pos = de.initial_pos
    desub.Ad1 = de.Ad1
    desub.Ad2 = de.Ad2
    desub.Ad3 = de.Ad3
    desub.entropy = de.entropy
    desub.first_col_names = de.first_col_names
    desub.output_type = 'ratio_d'
    desub.merge_from_info = False


class MyTestCase(unittest.TestCase):

    with open('../test-DnoisE/output_entropy.json', 'r') as doc_json:
        output_entropy = pd.read_json(doc_json, orient='table')

    with open('../test-DnoisE/output_no_entropy.json', 'r') as doc_json:
        output_no_entropy = pd.read_json(doc_json, orient='table')

    # entropy correction
    def test_entropy_correction(self):
        de = DnoisEFunctions()
        with open('../test-DnoisE/data.json', 'r') as doc_json:
            de.data_initial = pd.read_json(doc_json, orient='table')

        create_de_for_test(de)
        de.entropy = True
        de.compute_entropy = False
        de.cores = 1
        run_denoise_entropy(de)
        a = self.output_entropy != de.denoised_ratio_d
        del de
        self.assertEqual(a.any(axis=None), False) # add assertion here

    # entropy correction parallel
    def test_entropy_correction_par(self):
        de = DnoisEFunctions()
        with open('../test-DnoisE/data.json', 'r') as doc_json:
            de.data_initial = pd.read_json(doc_json, orient='table')

        create_de_for_test(de)
        de.entropy = True
        de.compute_entropy = False
        de.cores = 2
        run_denoise_entropy(de)
        a = self.output_entropy != de.denoised_ratio_d
        del de
        self.assertEqual(a.any(axis=None), False)  # add assertion here

    # without entropy correction
    def test_no_entropy_correction(self):
        de = DnoisEFunctions()
        with open('../test-DnoisE/data.json', 'r') as doc_json:
            de.data_initial = pd.read_json(doc_json, orient='table')

        create_de_for_test(de)
        de.entropy = False
        de.cores = 1
        run_denoise(de, test=True)
        a = self.output_no_entropy != de.denoised_ratio_d
        self.assertEqual(a.any(axis=None), False)  # add assertion here

    # without entropy correction parallel
    def test_no_entropy_correction_par(self):
        de = DnoisEFunctions()
        with open('../test-DnoisE/data.json', 'r') as doc_json:
            de.data_initial = pd.read_json(doc_json, orient='table')

        create_de_for_test(de)
        de.entropy = False
        de.cores = 2
        de.abund_col_names = []
        run_denoise(de, test=True)
        a = self.output_no_entropy != de.denoised_ratio_d
        del de
        self.assertEqual(a.any(axis=None), False)  # add assertion here

    # entropy correction calculation
    def test_entropy(self):
        de = DnoisEFunctions()
        with open('../test-DnoisE/data.json', 'r') as doc_json:
            de.data_initial = pd.read_json(doc_json, orient='table')

        create_de_for_test(de)
        de.entropy = False
        de.cores = 2
        de.abund_col_names = []
        run_denoise(de, test=True)
        a = self.output_no_entropy != de.denoised_ratio_d
        del de
        self.assertEqual(a.any(axis=None), False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
