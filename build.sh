${PYTHON} -m pip install . -vv --no-deps --ignore-installed

python -m nuitka --standalone --static-libpython=no ${SRC_DIR}/dnoise/DnoisE.py

mkdir -p ${PREFIX}/opt
mkdir -p ${PREFIX}/bin

mv DnoisE.dist ${PREFIX}/opt

# Set environment variable to skip RPATH checking
export CONDA_BUILD_SKIP_RPATH_CHECK=1

ln -s -r ${PREFIX}/opt/DnoisE.dist/DnoisE.bin ${PREFIX}/bin/DnoisE.bin
