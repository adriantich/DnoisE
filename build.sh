$PYTHON setup.py install --single-version-externally-managed --record=record.txt

cd ./src

python -m nuitka --enable-plugin=multiprocessing --static-libpython=no --follow-imports DnoisE.py


