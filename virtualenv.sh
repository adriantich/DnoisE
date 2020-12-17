#!/bin/bash

# create a virtualenv with python=>3.6 in DnoisE directory
virtualenv 3.6-dev ./venv

source ./venv/bin/activate

pip3 install pandas
pip3 install tqdm
pip3 install biopython
pip3 install python-Levenshtein

deactivate
