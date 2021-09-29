# __DnoisE__: Distance denoise by Entropy
## __An open source parallelizable alternative to Unoise__
by Adrià Antich (CEAB-CSIC, Center of Advanced Studies of Blanes)

Here we present a new program to denoise sequence data sets from Illumina using a distance parameter d corrected 
(optionally) according to the entropy of each codon position. DnoisE is a denoising software that uses the Unoise 
algorithm (Edgar 2016) to detect incorrect sequences from PCR and sequencing errors. The incorrect (“daughter”) 
sequences are merged with the correct “mother” sequence. For coding sequences where the entropy of each codon 
position is highly variable, a correction is advisable to avoid merging correct sequences that have changes in 
position 3 of the codons (highly variable in nature). DnoisE has been tested with the Leray fragment of the COI 
barcode region in Antich et al. (2021).


Pros of DnoisE versus Unoise:

1 - DnoisE can weight distances depending on the codon position of nucleotides where changes occur, based on entropy 
values of each codon position or any other user-settable measure of variability.

2 - DnoisE algorithm is parallelizable leading to high computational speed depending on computational hardware. 
It is a very good option if a multicore computer is available.

3 - DnoisE is written in Python3 and open access, which makes it user customizable.

4 - It accepts both .csv, .fasta and .fastq input files and can return .csv and .fasta file-types too.

5 - DnoisE allows the user to choose among three joining methods. Following Edgar's equation (beta(d)=.5^(alpha\*d+1)), 
the Unoise algorithm joins incorrect “daughter” sequences to the most abundant “mother” sequence with which they have 
an abundance ratio below beta(d). From our point of view this can lead to over-joining of sequences to the most 
abundant ones. Our algorithm returns two extra types of joining criteria outputs. For a given sequence, all potential 
“mothers” that satisfy the condition abundance skew ratio<beta(d) are stored. We then choose the correct “mother” as 
(1) the one having the lowest skew ratio with the potential “daughter” (ratio criterion, corresponding to the original 
Unoise formulation); (2) the “mother” with which it has the lowest d (distance criterion), or (3) the “mother” for which
the skew abundance ratio divided by beta(d) is the lowest (ratio_distance criterion). These criteria are short-named r, 
d, and r_d criteria.

### __INSTALLING DnoisE__

DnoisE depends on the following software:

* Python>=3.6
* C (for Levenshtein module)
* R
* bash

DnoisE is easy to install from terminal. Presently it runs on UNIX/Linux systems. An executable file can be optionally created (see below), otherwise the program can be run directly from the python script DnoisE.py.

#### __A. INSTALLATION FROM GITHUB__
#### __1. First clone the github repository__

```console
git clone https://github.com/adriantich/DnoisE.git
```

#### __2a. INSTALL WITH install.sh__
install.sh is a bash script that installs DnoisE and the required modules and creates a "./bin" directory with an executable file (DnoisE). It is recommended
to use a pyenv environment to avoid potential issues with python versions (see pyenv [documentation](https://github.com/pyenv/pyenv)).

```console
cd DnoisE/
bash install.sh
```

#### __2b. INSTALL MANUALLY__
The manual installation requires the following modules

   + [Pandas](https://pandas.pydata.org/)
   + [pyinstaller](https://www.pyinstaller.org/)
   + [python-Levenshtein](https://pypi.org/project/python-Levenshtein/)
   + [tqdm](https://pypi.org/project/tqdm/)

 To make a stand-alone application (optional), pyinstaller creates an executable file (DnoisE) in a "./bin" directory as follows:

```console
cd ./DnoisE/

pip3 install pandas
pip3 install pyinstaller
pip3 install python-Levenshtein
pip3 install tqdm

cd ./src

pyinstaller DnoisE.py --onefile --distpath ../bin
```
#### __B. INSTALLATION FROM CONDA__
#### __1. Install package__

```console
conda install -c adriantich dnoise
```
Different python versions (3.6, 3.7 and 3.8) are available at the [conda repository](https://anaconda.org/adriantich/dnoise/files)

#### __2. Create an executable file (optional)__
To create an executable file (in a "bin" directory) of the program run 

```console
pyinstaller DnoisE.py --onefile --distpath ../bin
```
from the directory where the conda package has been installed (typically ~/anaconda3/pkgs/dnoise-<dnoise_version>/lib/<python_version>/site-packages/src).


### __WORKFLOW__

#### __Running DnoisE__

DnoisE has been created to run both after and before clustering as described in Antich et al. (2021).
Therefore, it accepts and returns both .fasta and .csv files. However, writing output is faster if .csv. 
It also accepts .fastq	input files.

DnoisE can be called from the executable created by *pyinstaller* or directly from
the python script. See example below:

```console
> ./bin/DnoisE -h
> python3 ./src/DnoisE.py -h
```

Parameters of DnoisE are described in help but some are explained in more detail below.

```console
> ./bin/DnoisE -h

*HELP*
Displaying help
		-h --help display help
	Input file options:
		--csv_input [path] input file path in csv format
		--fasta_input [path] input file path in fasta format
		--fastq_input [path] input file path in fastq format
		--joining_file [path] file path of an info output from DnoisE. This option allows to use the information of previous runs of DnoisE to return different joining criteriaoutputs without running all the program again
		-n --count_name [size/reads/count...] count name column 'size' by default
		-p --sep [1/2/3] separator in case of csv input file
				1='	' (tab)
				2=','
				3=';'
		-q --sequence [sequence/seq...] sequence column name, 'sequence' by default
		-s --start_sample_cols [number] first sample column (1 == 1st col) if not given, just one column with total read counts expected (see README.md)
		-z --end_sample_cols [number] last sample column (n == nst col) if not given, just one column with total read counts expected (see README.md)
	Output file options:
		--csv_output [path] common path for csv format
		--fasta_output [path] common path for fasta format
		-j --joining_criteria [1/2/3]
				1-> will join by the lesser [abundance ratio / beta(d)] (default)
				2-> will join by the lesser abundance ratio (r criterion)
				3-> will join by the lesser distance (d) value (d criterion)
				4-> will provide all joining criteria in three different outputs (all)
	Other options:
		-a --alpha [number] alpha value, 5 by default
		-c --cores [number] number of cores, 1 by default
		-e --entropy [number,number,number] entropy values (or any user-settable measure of variability) of the different codon positions [0.47,0.23,1.02] by default
		-m --modal_length [number] when running DnoisE with entropy correction, sequence length expected can be set, if not, modal_length is used and only sequences with modal_length + or - 3*n are accepted
		-u --unique_length only modal length is accepted as sequence length when running with entropy correction
		-x --first_nt_codon_position [number] as DnoisE has been developed for COI sequences amplified with Leray-XT primers, default value is 3 (i.e., the reading frame starts in the third nucleotide).
		-y --entropy_correction a distance correction based on entropy is performed (see ENTROPY CORRECTION below). If set to F, no correction for entropy is performed (corresponding to the standard Unoise formulation)
```


__*INPUT FILES (--csv_input|--fasta_input|--fastq_input|-n|-q|-p)*__

Input files can be in either .csv, .fasta and .fastq format. This can be specified using *-f* parameter set as T 
as default meaning that input file is in .fasta format. All sequences are supposed to be aligned.

Different pipelines use different names for the number of reads (size/count/reads...). This can be specified using 
parameter *-n* followed by a string weith the name (for instance: -n size, default).  Sequence name can also be 
specified using *-q* parameter (sequence/seq...)(*-q* sequence, default)

If input is a fasta file both id (first qualifier) and size must end by ";". Any other qualifiers will be ignored.

If input file is a .csv, the separator between columns can be specified using the *-p* parameter (see help).

Examples of *how to run* and input files are available in the test-DnoisE subdirectory.


__*OUTPUT FILES (--csv_output|--fasta_output|-j)*__

Outoput files can be both in .csv and .fasta files but the path and a common pattern identifier 
(i.e. "/PATH/TO/DIR/common_pattern") must be specified. However, writing output is faster if .csv is chosen, especially 
in case of large files.

DnoisE can return three different types of output files. When a "daughter" sequence is found, different joining 
criteria can be applied if more than one possible "mother" meets Edgar’s equation requirement. As comparisons are 
done sequentially from higher to lower abundances, when a sequence meets a "mother", comparisons will stop if r 
criteria (*-j* 2) is chosen (i.e., the lesser ratio value between "daughter" abundance and "mother" abundance). 
However, this doesn't guarantee that the "mother" found is the best. 

If d criteria (*-j* 3 , lesser d value), r_d criteria (*-j* 1, lesser value of the skew abundance ratio divided by beta(d)), or all criteria (*-j* 4), are chosen, comparisons continue and all potential “mothers” are stored. The program will choose afterwards the best “mother” according to the preferred criterion.
Therefore, when r criteria is chosen, computation time is lower because less comparisons are performed. However, we recommend the r_d criteria which is set as default value.


__*MERGING FROM INFO FILE (--joining_file)*__

In order to return different types of output criteria (*-j*), even after DnoisE is finished, DnoisE creates an info file 
which contains information on whether a sequence is correct or if not, and how each sequence is merged to a mother with 
the proper distance value. This output can be re-analysed again by running DnoisE if the info file is specified with 
*--joining_file* and *-j* as desired. Note that if *-j* 2 has been initially chosen, it is not possible to obtain the 
other joining criteria without re-running DnoisE, because comparisons are halted when the first “mother” is encountered.


__*I/O PATHS AND SAMPLE SPECIFICATIONS (-s|-z)*__

Path before file name is required but can be in the form of ./ to avoid large strings.
Output path is also required to specify output name.
For instance, the output path './file_to_denoise will' will return a file './file_to_denoise_denoising_info.csv' 
among others.

Sample information cannot be processed if the input is a fasta or a fastq file. When this information is relevant, 
the input should be a .csv file, and the first and last sample columns should be indicated with parameters 
*-s* and *-z*.


__*ENTROPY CORRECTION (-e|-y|-x|-m|-u)*__

As described in Antich et al. (2021) a correction of the distance value (d) in Edgar's algorithm (2016) can be 
performed using the entropy values of each codon position in coding barcodes. We performed DnoisE for COI 
Leray/Leray-XT primers (Leray et al. 2013; Wangensteen et al. 2018) and consequently sequences start with a 
codon position 2 and the first initial codon position is the third nucleotide as follows:

!!!!!!!!!!!PER REVISAAR!!!!!
```console
seq       --> T-T-T-G-A-G-T-T-C-A-A-T-...
position  --> 2-3-1-2-3-1-2-3-1-2-3-1-...
```
*-x/--first_nt_position* is set as 3 by default.

Note that, in Edgar’s formula, the d used is the Levenshtein distance. This is the one used by DnoisE if no correction is applied. However, when entropy correction is selected, the Levenshtein distance is not applicable as the codon position needs to be considered, and a d based simply on the number of nucleotide differences is used instead. With sequences of equal length and aligned, both distances are practically equivalent.

The use of Levenshtein distance allowed us to compare sequences of inequal length, both in the complete dataset or within 
MOTUs (depending on whether DnoisE is performed before or after clustering). However, with the entropy correction length 
should be constant when comparing two sequences. Dataset is thus analysed by separate sequence length sets. These sets 
must differ from the modal length of all dataset by *n* number of codons (groups of 3 nucleotides). The modal 
sequence length can be set using *-m* parameter (i.e. in case of sequence lengths of 304, 310, 311, 313, 314 and 316 
nucleotides, if the modal length is 313, lengths of 311 and 314 will be considered as incorrect and will be eliminated). If *-u*,
only modal length is used and all sequences with different depth are not included.

Entropy values are given as E_1, E_2, E_3, where 1, 2, and 3 are the codon positions (default as *-e* 0.47,0.23,1.02). 
By default, the program computes the entropies from the complete dataset. Any user-derived value of variability of 
each codon position can be used instead of entropy.

The correction is applied as follows:

The original Edgar’s formula is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\beta&space;(d)=(1/2)^{\alpha&space;*&space;d&plus;1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta&space;(d)=(1/2)^{\alpha&space;*&space;d&plus;1}" title="\beta (d)=(1/2)^{\alpha * d+1}" /></a>

We correct the d value as:

<a href="https://www.codecogs.com/eqnedit.php?latex=d&space;=&space;\sum\limits_{i=1}^{3}&space;d_i&space;*&space;\frac{E_i&space;*&space;3}{E_1&space;&plus;&space;E_2&space;&plus;&space;E_3}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?d&space;=&space;\sum\limits_{i=1}^{3}&space;d_i&space;*&space;\frac{E_i&space;*&space;3}{E_1&space;&plus;&space;E_2&space;&plus;&space;E_3}" title="d = \sum\limits_{i=1}^{3} d_i * \frac{E_i * 3}{E_1 + E_2 + E_3}" /></a>

Entropy is computed by DnoisE if no entropy values are given when *-y*. In practice, when entropy is computed, 
*-x* is not mandatory so the program will compute three independent entropy values associated with codon 
positions of the current dataset.


#### __Running DnoisE after SWARM within MOTU__

If DnoisE is run after SWARM (see [Torognes](https://github.com/torognes/swarm)) a separate .csv file for each MOTU is needed.
MOTUs_from_Swarm.sh will return a directory were all MOTUs will be stored as separate .csv files

```console
> bash DnoisE/src/MOTUs_from_Swarm.sh -h

Generating a .csv file of each MOTU sequences using output of SWARM

Syntax: bash MOTUs_from_SWARM.sh [-h] [-i] [motu_list_file] [-o] [output_file_from_SWARM] [-r] [-t] [sample_abundances] [-d] [output_directory] [-l] [output_file_from_lulu]
options:
h     Print this Help.
i     .txt containing MOTU ids for which to create .csv files
      example: cat motus_to_denoise.txt
                seq1_id
                seq2_id
                seq3_id
o     swarm output file
r     remove databases that will be created during the process in the output directory
t     .tab file containing sample information of original sequences (from obitab)
d     output directory
l     lulu corrected_sequences file, an output file from lulu (opcional)

```

__*MANDATORY INPUTS*__

* a file containing a list of MOTUs do denoise
* output file from SWARM
* a .tab file containing sample information of original sequences

__*OUTPUTS*__

* directory with a .csv file of each MOTU
* a database with all the files that have been created during the process if not '-r'

__*LULU FILE*__

* output file from lulu with deleted/corrected sequences

Some pipelines (e.g. [MJOLNIR](https://github.com/uit-metabarcoding/MJOLNIR)) use [LULU](https://github.com/tobiasgf/lulu) to merge incorrect MOTUs to correct ones. This script retrieves MOTU sequence composition from original MOTUs generated by SWARM. Therefore, a *-l* [output_file_from_lulu] file us mandatory if LULU is used, as it gives information on which MOTUs have to be merged according to LULU results.
