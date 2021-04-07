# __DnoisE__: Distance denoise by Entropy
## __An open source parallelizable alternative to Unoise__
by Adrià Antich (CEAB-CSIC, Center of Advanced Studies of Blanes)

Here we present a new program to denoise sequence data sets from Illumina using parameter d (distance) corrected (optionally) according to the entropy of each codon position. DnoisE is a denoising software that uses the Unoise algorithm (Edgar 2016) to detect incorrect sequences from PCR and sequencing errors. The incorrect (“daughter”) sequences are merged with the correct “mother” sequence. For coding sequences where the entropy of each codon position is highly variable, a correction is advisable to avoid merging correct sequences that have changes in position 3 of the codons (highly variable in nature). DnoisE has been tested with the Leray fragment of the COI barcode region in Antich et al. (2021).


Pros of DnoisE versus Unoise:

1 - DnoisE can weight distances depending on the codon position of nucleotides where changes occur, based on entropy values of each codon position or any other user-settable measure of variability.

2 - DnoisE algorithm is parallelizable leading to high computational speed depending on computational hardware. It is a very good option if a multicore computer is available.

3 - DnoisE is written in Python3 and open access which makes it user customizable.

4 - It accepts both .csv and .fasta input files and can return both file-types too.

5 - DnoisE allows the user to choose among three joining method. Following Edgar's equation (beta(d)=.5^(alpha\*d+1)), the Unoise algorithm joins incorrect “daughter” sequences to the most abundant “mother” sequence with which they have an abundance ratio below beta(d). From our point of view this can lead to over-joining of sequences to the most abundant ones. Our algorithm returns two extra types of joining criteria outputs. For a given sequence, all potential “mothers” that satisfy the condition abundance skew ratio<beta(d) are stored. We then choose the correct “mother” as (1) the one having the lowest skew ratio with the potential “daughter” (ratio criterion, corresponding to the original Unoise formulation); (2) the “mother” with which it has the lowest d (distance criterion), or (3) the “mother” for which the skew abundance ratio divided by beta(d) is the lowest (ratio_distance criterion). These criteria are short-named r, d, and r_d criteria.

### __INSTALLING DnoisE__

DnoisE depends on the following software:

* Python3
* C (for Levenshtein module)
* R
* bash

DnoisE is easy to install from terminal. 

1. First clone the github repository

```console
git clone https://github.com/adriantich/DnoisE.git
```

2. Install the required modules (python3)

   + [Pandas](https://pandas.pydata.org/)
   + [tqdm](https://pypi.org/project/tqdm/)
   + [python-Levenshtein](https://pypi.org/project/python-Levenshtein/)
   + [stats](https://pypi.org/project/stats/)

```console
cd DnoisE/
bash required_modules.sh
```

3. Install required packages (R)

   + [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
   + [entropy](https://cran.r-project.org/web/packages/entropy/) 
   
```console
# R session:
R
> install.packages("optparse")
> install.packages("entropy")
> install.packages("stringr")
```

We also recomend to use pyenv to create an environment to run DnoisE (see pyenv [documentation](https://github.com/pyenv/pyenv))

### __WORKFLOW__

#### __Running DnoisE__

DnoisE has been created to run both after and before clustering as described in Antich et al. (2021).
Therefore, it accepts and returns both .fasta and .csv files. However, writing output is faster if .csv.

Parameters of DnoisE are described in help but some are explained in more detail below.

```console
> python3 DnoisE/src/DnoisE.py -h

*HELP*
Displaying help
 -h --help display help
 -i --input input file path
 -o --output common output files path
 -P --part DnoisE can be run by parts, part 1 runs the main program and returns the specified output and a database where results are stored
                 part 2 can re-analyse this database and return further outputs without running again the program (see README.md)
                     - If part = 1 (default) runs normally and a directory and database is named as --output
                     - If part = 2 returns outputs from database
                         Part 2 requires --input, --output and --cores if necessary
 -f --fasta_input logical, if T (default), fasta file as input, if F, .csv as input
 -F --fasta_output logical, if T (default), fasta file as output, if F, .csv as output
 -j --joining_criteria   1-> will join by the lesser [abundance ratio / beta(d)] (r_d criterion) (default)
                         2-> will join by the lesser abundance ratio (r criterion)
                         3-> will join by the lesser d value (d criterion)
                         4-> will provide all joining criteria in three different outputs
 -c --cores number of cores, 1 by default
 -s --start_sample_cols first sample column (1 == 1st col) if not given, just one column with total read count expected (see README.md)
 -z --end_sample_cols last sample column (n == nst col) if not given, just one column with total read count expected (see README.md)
 -n --count_name count name column (size/reads/count...) 'size' by default
 -a --alpha alpha value, 5 by default
 -q --sequence sequence column name (sequence/seq...), 'sequence' by default
 -p --sep separation 1='        ' (tab)
                     2=','
                     3=';'
 -e --entropy entropy (or any user-settable measure of variability) of the different codon positions [0.47,0.23,1.02] by default
 -y --entropy_correction logical, if T, a distance correction based on entropy is performed (see ENTROPY CORRECTION below). If set to F, no correction for entropy is performed (corresponding to the standard Unoise formulation)
 -x --first_nt_codon_position as DnoisE has been developed for COI sequences amplified with Leray-XT primers, default value is 3
 -m --modal_length when running DnoisE with entropy correction, sequence length accepted can be set, if not, modal_length is used
```

__*INPUT FILES (-i|-f|-n|-q|-p)*__

Input files can be in both .csv and .fasta format. This can be specified using *-f* parameter set as T as default meaning that input file is in .fasta format.

Different pipelines use different names to number of reads (size/count/reads...). This can be specified using parameter *-n* followed by string name (for instance: -n size, default).  Sequence name can also be specified using *-q* parameter (sequence/seq...)(*-q* sequence, default)

If input is a fasta file, the sequence must be in a single line and both id (first qualifier) and size must end by ";". Any other qualifier will be ignored.


      >Seq_000000012;size=433081;
      TTTGAGTTCAATACAAAGTCATTCAGGAGCTGCTATTGACTTAGCTATCTTCAGTTTACATCTTTCAGGAGCTTCTTCGATTCTAGGAGCAATTAATTTTATTTCTACCATTATAAATATGCGAAATCCTGGACAAACATTTTATCGCATTCCTTTATTTGTTTGATCGATTTTCGTAACTGCTTTACTACTATTATTAGCAGTACCAGTTTTAGCAGGAGCTATTACCATGTTACTAACTGATCGTAATTTTAATACAGCCTTTTTTGACCCTTCTGGAGGTGGTGATCCTGTACTTTATCAACATTTATTT

If input file is a .csv (*-f* F), the separator between columns can be specified using the *-p* parameter (see help).


__*OUTPUT FILES (-j|-F|-P)*__

DnoisE can return three different types of output files. When a "daughter" sequence is found, different joining criteria can be applied if more than one possible "mother" meets Edgar’s equation requirement. As comparisons are done sequentially from higher to lower abundances, when a sequence meets a "mother", comparisons will stop if r criteria (*-j* 2) is chosen (i.e., the lesser ratio value between "daughter" abundance and "mother" abundance). However, this doesn't guarantee that the "mother" found is the best. 

If d criteria (*-j* 3 , lesser d value), r_d criteria (*-j* 1, lesser value of the skew abundance ratio divided by beta(d)), or all criteria (*-j* 4), are chosen, comparisons continue and all potential “mothers” are stored. The program will choose afterwards the best “mother” according to the preferred criterion.
Therefore, when r criteria is chosen, computation time is lower because less comparisons are performed. However, we recommend the r_d criteria which is set as default value.

In order to return different types of output for both criteria (*-j*) and format (*-F*; T set as default for .fasta files and F for .csv files) even after DnoisE is finished, DnoisE creates a database in output directory which contains information of how to join all incorrect sequences. Output can be re-analysed again by running DnoisE with the parameter *-P* set as 2 and resetting both *-j* and *-F* parameters as desired. Note that if *-j* 2 has been initially chosen, it is not possible to obtain the other joining criteria without re-running DnoisE, because comparisons are halted when the first “mother” is encountered.



__*INPUT AND OUTPUT PATHS (-i|-o)*__

Path before file name is required but can be in form of ./ to avoid large strings.
Output path is also required as far as it can be used to specify output name.
For instance, the output path './file_to_denoise will' will return a file './file_to_denoise_denoising_info.csv' among others.


__*ENTROPY CORRECTION (-e|-y|-x)*__

As described in Antich et al. (2021) a correction of the distance value (d) in Edgar's algorithm (2016) can be performed using the entropy values of each codon position in coding barcodes.
We performed DnoisE for COI Leray/Leray-XT primers (Leray et al. 2013; Wangensteen et al. 2018) and consequently sequences start with a codon position 3 and the first codon position is in the second sequence position as follow

```console
seq       --> T-T-T-G-A-G-T-T-C-A-A-T-...
position  --> 3-1-2-3-1-2-3-1-2-3-1-2-...
```
*-x/--first_nt_position* is set as 3 by default.

Note that, in Edgar’s formula, the d used is the Levenshtein distance. This is the one used by DnoisE if no correction is applied. However, when entropy correction is selected, the Levenshtein distance is not applicable as the codon position needs to be considered, and a d based simply on the number of nucleotide differences is used instead. With sequences of equal length and aligned, both distances are practically equivalent.

The use of Levenshtein distance allowed us to compare sequences of inequal length, both in the complete dataset or within MOTUs (depending on whether DnoisE is performed before or after clustering, see below). However, with the entropy correction length should be constant. If denoise is performed before clustering with the Leray fragment, only 313 bp-long sequences will be compared. If done after clustering, only sequences of the modal length within each MOTU will be compared.

Entropy values are given as E_1, E_2, E_3, where 1, 2, and 3 are the codon positions (default as *-e* 0.47,0.23,0.1.02). Any user-derived value of variability of each codon position can be used instead of entropy.

The correction is applied as follows:

The original Edgar’s formula is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\beta&space;(d)=(1/2)^{\alpha&space;*&space;d&plus;1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta&space;(d)=(1/2)^{\alpha&space;*&space;d&plus;1}" title="\beta (d)=(1/2)^{\alpha * d+1}" /></a>

We correct the d value as:

<a href="https://www.codecogs.com/eqnedit.php?latex=d&space;=&space;\sum\limits_{i=1}^{3}&space;d_i&space;*&space;\frac{E_i&space;*&space;3}{E_1&space;&plus;&space;E_2&space;&plus;&space;E_3}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?d&space;=&space;\sum\limits_{i=1}^{3}&space;d_i&space;*&space;\frac{E_i&space;*&space;3}{E_1&space;&plus;&space;E_2&space;&plus;&space;E_3}" title="d = \sum\limits_{i=1}^{3} d_i * \frac{E_i * 3}{E_1 + E_2 + E_3}" /></a>

Entropy can be calculated using the entrpy.R script as follows:

```console
Rscript --vanilla entrpy.R -i [input_file] -o [output_file] -f [TRUE|FALSE] 

# display help

Rscript --vanilla entrpy.R --help

Usage: entrpy.R [options]


Options:
        -i CHARACTER, --input=CHARACTER
                dataset file name of format .fa/.fasta/.csv with just id, size and sequence. If .csv, only ',' accepted

        -x NUMERIC, --first_nt_position=NUMERIC
                first nucleotide position, 3 by default

        -o CHARACTER, --output_name=CHARACTER
                output file name

        -h, --help
                Show this help message and exit

```

__*SAMPLE INFORMATION (-s|-z)*__

Sample information cannot be processed if the input is a fasta file. When this information is relevant, the input should be a .csv file, and the first and last sample columns should be indicated with parameters *-s* and *-z*.

#### __Running DnoisE after SWARM within MOTU__

If DnoisE is run after SWARM (see [Torognes](https://github.com/torognes/swarm)) a separate .csv file for each MOTU is needed.
MOTUs_from_Swarm.sh will return a directory were all MOTUs will be stored as separate .csv files

```console
> bash DnoisE/src/MOTUs_from_Swarm.sh -h

Generating a .csv file of each MOTU sequences using output of SWARM

Syntax: bash MOTUs_from_SWARM.sh [-h] [-i] [motu_list_file] [-o] [output_file_from_SWARM] [-r] [-t] [sample_abundances] [-d] [output_directory] [-l] [output_file_from_lulu]
options:
h     Print this Help.
i     .txt containing MOTU ids for which to create a .csv file will be created
      example: cat motus_to_denoise.txt
                seq1_id
                seq2_id
                seq3_id
o     output swarm
r     remove databases that will be created during process on the output directory
t     .tab file containing sample information of original sequences
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
