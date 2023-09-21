${PYTHON} -m pip install . -vv --no-deps --ignore-installed

python -m nuitka --standalone --static-libpython=no ${SRC_DIR}/dnoise/DnoisE.py

mkdir -p ${PREFIX}/opt
mkdir -p ${PREFIX}/bin

mv DnoisE.dist ${PREFIX}/opt

ln -s -r ${PREFIX}/opt/DnoisE.dist/DnoisE.bin ${PREFIX}/bin/DnoisE.bin
