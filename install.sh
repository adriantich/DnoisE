#!/bin/bash

python3 setup.py install

cd ./src

python -m nuitka --enable-plugin=multiprocessing --follow-imports DnoisE.py

