#!/bin/bash

python -m pip install .

python -m nuitka --static-libpython=no --standalone dnoise/DnoisE.py
