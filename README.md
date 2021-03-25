# DnoisE
## _Distance denoise by Entropy_
### _An open source parallelizable alternative to Unoise_
by Adrià Antich (CEAB-CSIC, Center of Advanced Studies of Blanes)

Here we present a new program to denoise sequence data sets from Illumina using d (distance) corrected by Entropy. DnoisE is a denoising software that uses the Unoise algorithm (Edgar 2016) to detect incorrect sequences from PCR and sequencing errors. For coding sequences where the entropy of nucleotides depends on codon position, a correction is needed to avoid merging of correct sequences that have changes in position 3 of the codons (highly variable in nature). DnoisE has been tested with the Leray fragment of the COI barcode region in Antich et al. (2021).

Pros versus Unoise:

1 - DnoisE can weight differences depending on the position of nucleotides in the codon based on Entropy values of each position.

2 - DnoisE algorithm is parallelizable leading to high computational speed depending on computational hardware. It is a very good option if a multicore computer is available.

3 - DnoisE is written in Python3 and open access which makes it user customizable.

4 - It accepts both .csv and .fasta input files and can return both file-types too.

5 - DnoisE allows the user to choose among three joining method. Following Edgar's equation (beta(d)=.5^(alpha\*d+1)), the Unoise algorithm joins incorrect “daughter” sequences to the most abundant “mother” sequence with which they have an abundance ratio below beta(d). From our point of view this can lead to over-joining of sequences to the most abundant ones. We have developed an algorithm that returns two extra types of joining criteria outputs. For a given sequence, all potential “mothers” that satisfy the condition abundance skew ratio<beta(d) are stored. We then choose the correct “mother” as (1) the one having the lowest skew ratio with the potential “daughter” (ratio criterion, corresponding to the original Unoise formulation); (2) the “mother” with which it has the lowest d (distance criterion), or (3) the “mother” for which the skew abundance ratio divided by beta(d) is the lowest (ratio_distance criterion). These criteria are short-named r, d, and r_d criteria.



      *HELP*
      -h --help Display help
       -i --input input file path
       -o --output common output files path
       -P --part Denoise can be run by parts, part 1 runs only possible mothers and part 2 runs the others
                           In part = 1 runs normally and a directory as database is named as --output (default)
                           In part = 3 returns outputs from database
                               Part 3 requires --input, --output and --cores if necessary
       -f --fasta_input logical, if T (default), fasta file as input, if F .csv as input
       -F --fasta_output logical, if T (default), fasta file as input, if F .csv as input
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
       -p --sep separation 1='        '
                           2=','
                           3=';'
       -e --entropy entropy of the different codon positions [0.4298,0.1833,0.9256] by default
       -y --entropy_influence logical, if T, Ad correction parameter is used. In this case, only sequences with mode length are processed

*INPUT FILES*

Input files can be in both .csv and .fasta format

If input is a fasta file, the sequence must be in a single line and both id and size must end by ";".


      >Seq_000000012;size=433081;
      TTTGAGTTCAATACAAAGTCATTCAGGAGCTGCTATTGACTTAGCTATCTTCAGTTTACATCTTTCAGGAGCTTCTTCGATTCTAGGAGCAATTAATTTTATTTCTACCATTATAAATATGCGAAATCCTGGACAAACATTTTATCGCATTCCTTTATTTGTTTGATCGATTTTCGTAACTGCTTTACTACTATTATTAGCAGTACCAGTTTTAGCAGGAGCTATTACCATGTTACTAACTGATCGTAATTTTAATACAGCCTTTTTTGACCCTTCTGGAGGTGGTGATCCTGTACTTTATCAACATTTATTT

If input file is a .csv, separation between columns can be specified using the -p parameter (see help)

