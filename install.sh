#!/bin/bash

# python -m pip install .

python -m nuitka --enable-plugin=multiprocessing --static-libpython=no --standalone --follow-imports dnoise/DnoisE.py
