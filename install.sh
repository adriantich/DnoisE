#!/bin/bash

python3 setup.py install

cd ./src

python3 -m nuitka --enable-plugin=multiprocessing --follow-imports DnoisE.py

