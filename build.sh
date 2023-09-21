${PYTHON} -m pip install . -vv --no-deps --ignore-installed

python -m nuitka --enable-plugin=multiprocessing --standalone --static-libpython=no --follow-imports ${SRC_DIR}/dnoise/DnoisE.py

mkdir -p ${PREFIX}/opt
mkdir -p ${PREFIX}/bin

mv DnoisE.dist ${PREFIX}/opt

ln -s -r ${PREFIX}/opt/DnoisE.dist/DnoisE.bin ${PREFIX}/bin/DnoisE.bin
