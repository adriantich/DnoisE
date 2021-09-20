#!/bin/bash

# create required modules with python=>3.6
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
pyenv shell 3.6-dev
pyenv virtualenv 3.6-dev DnoisE
pyenv shell DnoisE3.6

pip3 install pandas
pip3 install tqdm
pip3 install python-Levenshtein
pip3 install pyinstaller

cd ./src

pyinstaller DnoisE.py --onefile --distpath ../bin



