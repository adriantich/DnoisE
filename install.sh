#!/bin/bash

python3 setup.py install

cd ./src

pyinstaller DnoisE.py --onefile --distpath ../bin
