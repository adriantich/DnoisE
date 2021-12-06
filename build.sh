$PYTHON setup.py install --single-version-externally-managed --record=record.txt

cd ./src

python -m nuitka --enable-plugin=multiprocessing --follow-imports DnoisE.py


