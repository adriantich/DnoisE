from setuptools import setup

setup(
    name="DnoisE",
    version="1.2",
    description="Denoise sequence data sets from Illumina using distance corrected according to the entropy"
                " of each codon position",
    author="Adrià Antich",
    author_email="adriantich@gmail.com",
    url="https://github.com/adriantich/DnoisE",
    packages=["src"],
    install_requires=["pandas", "nuitka", "levenshtein", "tqdm"]
)

