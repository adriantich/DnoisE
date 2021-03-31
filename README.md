# __DnoisE__: Distance denoise by Entropy
## __An open source parallelizable alternative to Unoise__
by Adrià Antich (CEAB-CSIC, Center of Advanced Studies of Blanes)

Here we present a new program to denoise sequence data sets from Illumina using parameter d (distance) corrected (optionally) according to the entropy of each codon position. DnoisE is a denoising software that uses the Unoise algorithm (Edgar 2016) to detect incorrect sequences from PCR and sequencing errors. The incorrect (“daughter”) sequences are merged with the correct “mother” sequence. For coding sequences where the entropy of each codon position is highly variable, a correction is advisable to avoid merging correct sequences that have changes in position 3 of the codons (highly variable in nature). DnoisE has been tested with the Leray fragment of the COI barcode region in Antich et al. (2021).


Pros of DnoisE versus Unoise:

1 - DnoisE can weight distances depending on the codon position of nucleotides where changes occur, based on entropy values of each codon position.

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
 -h --help display help
 -i --input input file path
 -o --output common output files path
 -P --part DnoisE can be run by parts, part 1 runs the main program, but if crushes, 
               part 3 can return outputs from database
                     - If part = 1 runs normally and a directory as database is named as --output (default)
                     - If part = 3 returns outputs from database
                         Part 3 requires --input, --output and --cores if necessary
 -f --fasta_input logical, if T (default), fasta file as input, if F, .csv as input
 -F --fasta_output logical, if T (default), fasta file as output, if F, .csv as output
 -j --joining_type 1-> will join by the lesser [abundance ratio / beta(d)] (ratio_d) (default)
                   2-> will join by the lesser abundance ratio (ratio)
                   3-> will join by the lesser d value (d)
                   4-> will give all joining types in three different outputs
 -c --cores number of cores, 1 by default
 -s --start_sample_cols first sample column (1 == 1st col) if not given, just total read count expected
 -z --end_sample_cols first sample column (1 == 1st col) if not given, just total read count expected
 -n --count_name count name column (count/reads/size..) 'count' by default
 -a --alpha alpha value, 5 by default
 -q --sequence sequence column name, 'sequence' by default
 -p --sep separation 1='        ' (tab)
                     2=','
                     3=';'
 -e --entropy entropy of the different codon positions [0.4298,0.1833,0.9256] by default
 -y --entropy_influence logical, if T, Ad correction parameter is used. In this case, only sequences with mode length are processed
 -x --first_nt_position as far as DnoisE has been performed for COI barcode amplified with Leray-XT primers, default value is 3

```

__*INPUT FILES*__

Input files can be in both .csv and .fasta format

If input is a fasta file, the sequence must be in a single line and both id and size must end by ";".


      >Seq_000000012;size=433081;
      TTTGAGTTCAATACAAAGTCATTCAGGAGCTGCTATTGACTTAGCTATCTTCAGTTTACATCTTTCAGGAGCTTCTTCGATTCTAGGAGCAATTAATTTTATTTCTACCATTATAAATATGCGAAATCCTGGACAAACATTTTATCGCATTCCTTTATTTGTTTGATCGATTTTCGTAACTGCTTTACTACTATTATTAGCAGTACCAGTTTTAGCAGGAGCTATTACCATGTTACTAACTGATCGTAATTTTAATACAGCCTTTTTTGACCCTTCTGGAGGTGGTGATCCTGTACTTTATCAACATTTATTT

If input file is a .csv, separation between columns can be specified using the -p parameter (see help)

__*INPUT AND OUTPUT PATHS*__

Path before file name is required but can be in form of ./ to avoid larger strings.
Output path is also required as far as it can be used to specify output name.
For instance, the output path './file_to_denoise will' will return a file './file_to_denoise_denoising_info.csv' among others.

__*ENTROPY CORRECTION*__

As described in Antich et al. (2021) a correction to the importance of d in Edgar's algorith (2016) from entropy values of each nucleotide position in codon for coding barcodes.
We performed DnoisE for COI Leray/Leray-XT primers (Leray et al. 2013; Wangensteen et al. 2018) and therefore sequences starts with a position 3 and de first first codon position is in the second sequence position as follow
```console
seq       --> T-T-T-G-A-G-T-T-C-A-A-T-...
position  --> 3-1-2-3-1-2-3-1-2-3-1-2-...
```
--first_nt_position is set as 3 by default.

Entropy values are given as E1,E2,E3 independent from later parameter.
The correction is aplied as following.

If the original Edgar formula is:

$\beta(d) = (1/2)^{\alpha * d + 1}$

We have performed a correction to d value as following:

$d = \sum\limits_{i=1}^{3} d_i * \frac{E_i * 3}{E_1 + E_2 + E_3}$

Entropy values can be calculated using the entrpy.R script as following:

```console
Rscript --vanilla entrpy.R -i [input_file] -o [output_file] -f [TRUE|FALSE] 

# display help

Rscript --vanilla entrpy.R --help

Usage: entrpy.R [options]


Options:
        -i CHARACTER, --input=CHARACTER
                dataset file name

        -x NUMERIC, --first_nt_position=NUMERIC
                first nucleotide position, 3 by default

        -f LOGICAL, --fasta_input=LOGICAL
                input file is a .fasta file (TRUE) or .csv (FALSE, default)

        -o CHARACTER, --output_name=CHARACTER
                output file name

        -h, --help
                Show this help message and exit

```


#### __Running DnoisE after SWARM within MOTU__

If DnoisE is runned after SWARM (see [Torognes](https://github.com/torognes/swarm)) creating an independent .csv for each MOTU is needed.
MOTUs_from_Swarm.sh will return a directory were all MOTUs composition will be stored as independent .csv

```console
> bash DnoisE/src/MOTUs_from_Swarm.sh -h

Generating a .csv file of each MOTU sequences using output of SWARM

Syntax: bash MOTUs_from_SWARM.sh [-h] [-i] [motu_list_file] [-o] [output_file_from_SWARM] [-r] [-t] [sample_abundances] [-d] [output_directory] [-l] [lulu_deleted_seqs]
options:
h     Print this Help.
i     .txt containing MOTU ids for which to create a .csv file
      example: cat motus_to_denoise.txt
      		seq1_id
      		seq2_id
      		seq3_id
o     output swarm
r     remove databases that will be created during process on the output directory
t     .tab file containing sample information of original sequences
d     output directory
l     lulu delated sequences (opcional) 

```
__*MANDATORY INPUTS*__

* a file containing a list of MOTUs do denoise
* output file from SWARM
* a .tab file containing sample information of original sequences

__*OUTPUTS*__

* directory with .csv of each MOTU
* a database with all the files that have been created during the process if not '-r'
