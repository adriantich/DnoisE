$PYTHON setup.py install --single-version-externally-managed --record=record.txt

cd ./src

pyinstaller DnoisE.py --onefile --distpath ../bin

