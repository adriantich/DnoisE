#!/bin/bash

# create required modules with python=>3.6
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
pyenv shell 3.8.2
pyenv virtualenv 3.8.2 DnoisE
pyenv shell DnoisE

pip3 install pandas
pip3 install pyinstaller
pip3 install python-Levenshtein
pip3 install tqdm

python3 setup.py install

cd ./src

pyinstaller DnoisE.py --onefile --distpath ../bin
